# 実装方針(案): BLAST hit強度のスコアリング反映 (feature/sequence-evidence)

Step 1(現状分析・方針提案)のみ。コードはまだ変更していません。
ベースブランチ: `main`(`13d6ef2`, `feature/integrated-scoring` マージ済み・316テストpass確認済み)。
作業ブランチ: `feature/sequence-evidence`。

## 1. 背景(再掲)

統合設計書3章・9章の指摘:「非常に強いBLAST hitと弱いBLAST hitが、同じpositive hit
として同程度に扱われる可能性がある」。実際、`analysis/interaction_scoring.py` の
`_build_evidence_components_v2` 内 `source_classification` コンポーネントは

```python
source_value = min(1.0, CANDIDATE_PRIORITY_BASE.get(candidate_source, 10.0) / 30.0)
```

という「どの candidate_source バケツから来たか」という粗いカテゴリ値のみを使っており、
`ProteinRecord.positive_hits`(実際のBLASTヒットのidentity/coverage/evalue/bitscore)
を一切見ていません。`legacy_additive` 側の `_candidate_priority_score` も同様です。

## 2. 現状分析

### 2.1 `positive_hits` の実体

- [core/models.py](core/models.py) `BlastHit`: `percent_identity`(0-100)、
  `alignment_length`(int)、`evalue`、`bitscore`、`query_coverage`(プロパティ、
  `alignment_length / query_length * 100`。`query_length` が未設定なら `None`)。
- 生成経路: [analysis/blast_pipeline.py:81-108](analysis/blast_pipeline.py:81) が
  `run_blast_pipeline(query_fasta=target_fasta, subject_fasta=positive_fasta, ...)` を
  実行 → [blast/runner.py:86-131](blast/runner.py:86) の `run_blastp` が
  `-evalue {config.blast.evalue}`(`config.yaml` 既定値 `1e-5`)を **BLAST自体の
  カットオフとして** 渡している。つまり `positive_hits` に入る evalue は
  **既に必ず `<= evalue閾値`**(既定 1e-5)であり、そこから 0 に近い側まで
  指数的に分布する(0 が返るケースもBLAST仕様上あり得る)。
- `-max_target_seqs 10`(既定)なので、1候補につき最大10件のヒットが
  positive_fasta 側の複数配列と対応しうる。
- `config.yaml` の `ortholog_filter`(strong/medium/weak)が今回参考にできる
  「実際にこのプロジェクトで使われているidentity/coverage/evalue閾値」:

  | 強度 | min_identity | min_query_coverage | max_evalue |
  |---|---|---|---|
  | strong | 40.0 | 70.0 | 1e-5 |
  | medium | 30.0 | 70.0 | 1e-5 |
  | weak | 25.0 | 50.0 | 1e-3 |

  (これは negative hit 分類用だが、identity 25-40%・coverage 50-70% が
  「生物学的に意味のある下限」として既にこのプロジェクトで採用されている値。)
- `bitscore` は現状スコアリングに一切使われておらず、既存コード2箇所
  (下記2.2)で「複数ヒットから代表1件を選ぶ」ためのタイブレークにのみ使用。

### 2.2 複数 `positive_hits` の扱い(現状コード)

- [analysis/candidates.py:13-20](analysis/candidates.py:13) `group_hits_by_query`:
  BLASTヒットを `query_id`(=候補タンパク質のprotein_id)でグルーピングし、
  `ProteinRecord.positive_hits` にそのままリストとして格納。**集約はしていない**。
- 代表ヒットを選ぶ既存の2つの実装(**どちらも同じルール**):
  - [analysis/candidates.py:23-28](analysis/candidates.py:23) `get_best_hit`
  - [output/excel.py:663-668](output/excel.py:663) `_best_hit`

  ```python
  return max(hits, key=lambda hit: (hit.bitscore, -hit.evalue))
  ```

  bitscore降順、同点ならevalue昇順。Excel出力の `best_positive_hit` /
  `best_positive_bitscore` / `best_positive_evalue` 列は既にこれを使っている。
- [analysis/ortholog_filter.py:145-159](analysis/ortholog_filter.py:145)
  `_is_better_representative` は negative hit 用に「強度ランク優先、同ランク内は
  (bitscore, -evalue)」という別ルールを使っているが、これは「strong/medium/weak
  という段階分類が先にある」ネガティブ評価特有のロジックであり、positive側には
  段階分類の概念がまだない。

**結論**: 新規コンポーネントでは `get_best_hit(candidate.positive_hits)` を
そのまま再利用し、代表1ヒットのidentity/coverage/evalue/bitscoreを正規化対象と
する。新しい集約ロジックを発明しない(複数ヒットの平均を取る、等はしない)。

### 2.3 正規化の土台

[core/evidence.py](core/evidence.py) に既に `clamp01()` と
`linear_normalize(value, low, high)` がある(genomic_context等で使用中)。
対数変換のヘルパーは現状ないが、`linear_normalize` に `-log10(evalue)` を
通せばそのまま流用できるので、`core/evidence.py` 自体の変更は不要と見込む。

## 3. 正規化方式(案)

`config/scoring_engine.example.yaml` に新セクション `sequence_evidence:` を追加し、
`analysis/scoring_engine_config.py::ScoringEngineConfig` を拡張して読み込む
(既存の `category_caps` 等と同じパターン)。

### 3.1 identity / coverage -- `linear_normalize` を直接適用

```python
identity_score = linear_normalize(hit.percent_identity, identity_floor, identity_ceiling)
coverage_score = linear_normalize(hit.query_coverage, coverage_floor, coverage_ceiling)
```

- 既定値案: `identity_floor=25.0`(weak閾値と同じ)、`identity_ceiling=90.0`。
  `coverage_floor=50.0`、`coverage_ceiling=100.0`。
- `hit.query_coverage` が `None`(`query_length` 未設定)の場合は
  coverage サブシグナルを「このヒットでは計算不能」として3.4の加重平均から
  除外する(0点扱いにしない)。

### 3.2 evalue -- 対数変換してから `linear_normalize`

evalueは指数分布かつ既にBLASTのevalueカットオフで上限が切られているため、
素の値では正規化できない。

```python
safe_evalue = max(hit.evalue, evalue_log_floor_guard)  # 例: 1e-300, log10(0)回避
evalue_score = linear_normalize(
    -math.log10(safe_evalue),
    -math.log10(evalue_reference_ceiling),  # 例 1e-5 -> 5.0
    -math.log10(evalue_reference_floor),    # 例 1e-100 -> 100.0
)
```

- `evalue_reference_ceiling` は「弱いと見なす境界」(既定案 `1e-5`。BLAST自体の
  カットオフと同じ値をデフォルトにする=「カットオフぎりぎりのヒットは
  スコア0」という自然な意味になる)。
- `evalue_reference_floor` は「これ以上強くても頭打ち」の境界(既定案 `1e-100`)。
- `evalue == 0.0`(BLASTが浮動小数点精度限界で0を返すケース)は
  `evalue_log_floor_guard`(既定案 `1e-300`)に丸めてから対数を取る
  (`math.log10(0)` は `ValueError` になるため必須のガード)。
- これらの数式・閾値はすべてYAML設定値であり、コードにハードコードしない
  (設計書50.3章の方針どおり)。

### 3.3 bitscore -- 要ユーザー判断(下記4章で確認したい点)

bitscoreは配列長・スコア行列に依存するため、identity/coverage/evalueと違って
「このプロジェクトで既に妥当性検証済みの絶対閾値」が存在しない
(`ortholog_filter` の閾値はidentity/coverage/evalueのみで、bitscoreを
使っていない)。単純な線形正規化はミスリーディングになりうる。

選択肢:
- **(a) 今回は含めない**(identity+coverage+evalueの3軸のみ)。将来
  「bits per aligned residue」等の長さ非依存な指標を別途検討する余地を残す。
- **(b) `bitscore / alignment_length` を4つ目のサブシグナルとして実装するが、
  既定サブ重みを `0.0` にして無効化しておく**(監査目的で値は出すが、
  デフォルトではスコアに寄与しない。ユーザーが将来値を入れれば有効化される)。

(a)はシンプルだが「後で追加」のための配線コストが発生する。(b)は
今回の実装コストがわずかに増えるが、無根拠な閾値をデフォルトで効かせない
まま拡張性を確保できる。**現時点では(b)を推奨**するが、この判断は
ユーザー確認事項とする(4章)。

### 3.4 サブシグナルの合成

各サブシグナル(identity/coverage/[evalue]/[bitscore])に設定可能なサブ重みを
持たせ、**このヒットで計算できたサブシグナルのみ**の加重平均を
`sequence_evidence` コンポーネントの `normalized_value` とする
(`analysis/scoring_engine.py::_score_categories` が「利用可能な証拠だけで
正規化する」のと同じ考え方を、コンポーネント内のサブシグナルレベルでも
踏襲する)。

## 4. カテゴリ配置案(A/B)

### 案A: `source_classification` カテゴリに2つ目のコンポーネントとして追加(推奨)

`co_occurrence` + `domain_complementarity` が `functional_annotation` の
20点キャップを共有している既存パターンと同一。

- 賛成理由:
  - `candidate_source`(Candidates/Candidates_relaxedなど)自体が
    「positive_hitsの有無・negative_hitsの有無」から機械的に決まる分類
    ([analysis/blast_pipeline.py:132-178](analysis/blast_pipeline.py:132))
    であり、`source_classification` と新しい「ヒット強度」信号は
    **同じBLAST証拠の別の側面**(カテゴリ=どのバケツか、強度=どれだけ強いか)。
    別カテゴリに分けると、実質同じ証拠源を2つの独立カテゴリの分母として
    二重に数えるリスクがある。
  - 設定ファイル・`DEFAULT_CATEGORY_CAPS` の変更が不要(新規カテゴリを
    追加しないため、`category_caps` のキー自体は増えない)。
  - 実装量が少ない(`V2_COMPONENT_WEIGHTS["source_classification"]` に
    2つ目のエントリを足すだけ)。
- 懸念点:
  - `source_classification` の元々の値(候補ソースからの固定値、常に
    `AVAILABLE`)と平均化されるため、「候補ソースは高優先度だがヒットは
    弱い」ケースのスコア低下幅が、独立カテゴリにするより小さくなる
    (2成分の単純平均なら最大で該当カテゴリ30点の半分=15点分の差)。

### 案B: 新規カテゴリ `sequence` として独立させる

- 賛成理由:
  - 設計書3章・9章が名指しした問題(強いhitと弱いhitの区別)を、
    既存の `source_classification` 信号から完全に独立させて可視化できる。
    Excel監査時に「このペアはBLAST hit強度が弱かったせいでスコアが
    下がった」ことが一目でわかる。
  - `analysis/scoring_engine.py` の設計(`total_cap` は
    アクティブなカテゴリのcapの合計を都度計算し、`final_score` は
    その合計に対する比率として正規化される)により、新カテゴリを
    追加してもTier閾値(70/50/25点、0-100スケール)の再校正は
    **不要**(既存のPIHカテゴリ追加時と同じ理由で、スケールが自動的に
    保たれる)。
- 懸念点:
  - `DEFAULT_CATEGORY_CAPS` / `config/scoring_engine.example.yaml` に
    新規キー追加が必要、かつ両者の一致を守るテスト
    (`test_scoring_engine_config.py::test_example_config_matches_defaults`)
    の対象に含める必要がある。
  - 新しいcap値(何点にするか)自体が未校正の新しい判断を要求する
    (例えば15点?20点?)。

### 推奨

**案A**(`source_classification` 内で共有)を推奨します。理由は
「同じBLAST証拠を測る2つの側面」という性質上、`co_occurrence` /
`domain_complementarity` の前例に最も忠実であり、変更範囲も最小になるためです。
ただし案Bの「問題を可視化しやすい」というメリットも実質的な価値があるため、
**この選択はユーザー確認事項とします**(4章参照)。

## 5. `positive_hits` が空の場合の扱い

`No_hit` ソースや、まだBLASTヒットが1件もない候補では `positive_hits == []`。
このとき新コンポーネントのstatusをどうするか:

- `MISSING`: 「必要な入力データが存在しなかった」。
  [analysis/interaction_scoring.py:930-931](analysis/interaction_scoring.py:930)
  `_co_occurrence_status_and_value`(queryレコードなし)や
  [analysis/interaction_scoring.py:958-959](analysis/interaction_scoring.py:958)
  `_domain_complementarity_status_and_value`(説明文もドメインもなし)が
  この文言で使っているステータス。
- `NOT_APPLICABLE`: 「この比較は意味をなさない」。
  [analysis/interaction_scoring.py:868](analysis/interaction_scoring.py:868)
  `_negative_hit_status_and_value` が「negative hitがない」場合に使用。

**この2つの先例は矛盾しているように見えます**(「BLASTを実行したが0件」という
状況は、negative側では `NOT_APPLICABLE`、他の多くの箇所では `MISSING` 的に
扱われている)。今回は「BLAST自体は実行されたが該当データがない」という点で
`MISSING` がより一貫すると考えますが、`negative_hit_strength` との対称性
(「証拠が存在しない」という同種の状況)を重視するなら `NOT_APPLICABLE` も
妥当です。**この判断もユーザー確認事項とします**(4章)。

いずれにせよ、0点ではなくこの2つのいずれかとして扱うことで、「評価したが
弱い」(`AVAILABLE`, 低い`normalized_value`)と「評価できなかった」
(`MISSING`/`NOT_APPLICABLE`, スコア分母から除外)を区別します。

## 6. ユーザー確認事項(実装着手前に回答をお願いします)

1. **bitscoreサブシグナル**: 3.3節の (a) 今回は含めない / (b) 実装するが
   既定重み0で無効化 -- どちらにしますか?
2. **カテゴリ配置**: 4章の 案A(`source_classification`共有・推奨) /
   案B(新規`sequence`カテゴリ) -- どちらにしますか?
3. **`positive_hits`空時のstatus**: 5章の `MISSING`(推奨) /
   `NOT_APPLICABLE` -- どちらにしますか?
4. **identity/coverage/evalueの既定floor/ceiling値**: 3.1/3.2節の提案値
   (identity 25-90、coverage 50-100、evalue参照範囲 1e-5〜1e-100)で
   よいか、調整したい値があれば教えてください。

## 7. Step 2 で touch する見込みのファイル(未実装・参考情報)

- `analysis/scoring_engine_config.py`: `sequence_evidence` 設定ブロック
  (floor/ceiling/サブ重み)を追加。
- `config/scoring_engine.example.yaml`: 上記の既定値を追加。
- `analysis/interaction_scoring.py`: `_build_evidence_components_v2` に
  新コンポーネントを追加する分岐、新関数(仮称)
  `_sequence_evidence_status_and_value(candidate, engine_config)` を追加。
  `V2_COMPONENT_WEIGHTS` のコメント更新。`get_best_hit` を
  `analysis/candidates.py` からimportして再利用。
- `tests/test_interaction_scoring.py`: 強いhit/弱いhit/hitなし
  (MISSINGまたはNOT_APPLICABLE)/複数hit のケースを追加。
- `tests/test_scoring_engine.py`: zero weight のユニットテスト、
  multiple negative(複数negativeコンポーネントの合算後cap)のユニットテスト
  を追加(41章チェックリスト残り2件)。
- `docs/scoring_engine_v2.md`, `CHANGELOG.md`: 変更内容を追記。

`analysis/ortholog_filter.py` と `legacy_additive` 経路
(`_score_pair`, `_candidate_priority_score` 等)は変更しません。
