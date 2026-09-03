# Final Score統合フェーズ:調査と設計提案・実装記録

Status: **実装済み・修正済み**(方式A、cap 30/70。negative_hit_strength
ペナルティは実装後の実データ検証で問題が判明し、Final Scoreからは除外して
NOT_APPLICABLE予約枠に変更済み)。design spec §17-22・§27が求める、
`protein_hunter_score`と`interaction_score`を1つの数値へ統合する
"Final Score"について、調査→承認→実装→実データ検証→設計修正→再検証、の
記録。Phase 6/7/8(Excel/Word再設計、design spec独自番号)より前に着手した。

以下、§1〜§6は実装前の調査内容(変更なし、記録として保持)。実装後の経緯は
文末の「実装後の実データ検証」「negative_hit_strengthペナルティの除外」
「修正後の再検証」を参照。

## 1. 値域・スケール:最重要の発見

**`protein_hunter_score`と`interaction_score`は全く異なるスケールです。**

- **`protein_hunter_score`**: `analysis/scoring.py::build_candidate_score`
  (`DEFAULT_WEIGHTS`: positive_hit=5, no_negative_hit=5, domain_hit=4,
  uniprot_accession=2, alphafold_url=2, annotation_warning=-1)の単純加算。
  理論上の最大値 = 5+5+4+2+2 = **18**(uniprot/alphafoldは実データでは
  ほぼ付かないため、実質的な上限は14前後)。理論上の最小値は-1(通常は0)。
  **正規化・再スケーリングは一切されていません**——`resolve_protein_hunter_scores`
  (`analysis/interaction_scoring.py:656-683`)は`build_candidate_score(record)`
  を無変換で呼ぶだけです。
- **`interaction_score`**: v2は`analysis/scoring_engine.py::score_candidate`
  のカテゴリcap正規化により`output_scale=100`(既定値)に対して常に**0-100**
  スケールで再計算されます。legacy_additiveも`_legacy_interaction_score`
  (`analysis/interaction_scoring.py`)が`raw/max_points*100`で同じく0-100に
  正規化しています。両モデルとも0-100スケールです。

### 実データでの分布(直近の広範囲診断run、76,797ペア。11クエリ×全バケツ有効)

| 指標 | count | mean | median | min | max |
|---|---:|---:|---:|---:|---:|
| `protein_hunter_score` | 76,797 | 8.45 | 9.00 | 4.00 | 14.00 |
| `interaction_score` | 76,797 | 9.01 | 8.93 | 0.00 | 69.68 |
| `interaction_priority_score`(参考) | 76,797 | 15.93 | 18.28 | 0.00 | 65.74 |

`Candidates`シート(151件、`total_score`列=`protein_hunter_score`と同一式)
では**147/151件が値14で完全に一致**、4件が10——理論上の値域[-1,18]の中で
実質2値にほぼ収束しています。`interaction_score`は0〜70弱まで連続的に分布。
**単純に足し算・平均するとスケールの大きい`interaction_score`が支配的になり、
`protein_hunter_score`はほぼ定数のノイズとしてしか働かない**ことが、この
分布だけからも予見できます(§6の実データ検証で実際にそうなることを確認)。

## 2. 両スコアリングモデルで計算されているか

**両方とも`legacy_additive`・`v2_evidence_based`の両モデルで計算されています。**

- `protein_hunter_score`: `resolve_protein_hunter_scores(config, blast_classification)`
  の呼び出し(`analysis/interaction_scoring.py:518`)は`scoring_model`分岐の
  外側にあり、無条件に実行されます。`INTERACTION_PAIR_COLUMNS`のコメントも
  "Present for both scoring models" と明記(同ファイル273行目付近)。
- `interaction_score`: v2は`_interaction_only_breakdown`経由でEvidenceComponent
  ベースに計算、legacyは`_legacy_interaction_score`が独自に同じ0-100正規化を
  再実装(`INTERACTION_SCORE_COMPONENT_NAMES`と同じ対象範囲:
  `gene_neighborhood`+`domain_complementarity`[+`external_ppi`が有効なら]の
  重みのみを分母とする)。ただし**`coexpression_gse77738`除外の変更(PR #11)は
  v2のみに適用済みで、legacyの`_legacy_interaction_score`は元々coexpression
  を一切含んでいない**(STRING由来の`string_ppi_score`のみ)ため、この点は
  自然に整合しています。

## 3. `negative_hit_strength`ペナルティの適用タイミング・範囲

`analysis/scoring_engine.py::score_candidate`内、`final_score = clamp(raw_score
- negative_penalty_points, 0, output_scale)`(104-107行目)で、**カテゴリ正規化
後・最終スケール(0-100)に対して直接減算**されます。ペナルティ自体は
`V2_COMPONENT_WEIGHTS["negative_hit_strength"] = 30.0`を直接output_scale
ポイントとして使い、`negative_penalty_cap`(既定30.0)で上限。

**重要な非対称性**: `negative_hit_strength`は`INTERACTION_SCORE_COMPONENT_NAMES`
に含まれていません。したがって:
- **`interaction_priority_score`(フル合成)には適用される**
- **`interaction_score`(クエリ固有のみ)には適用されない**(そもそも
  `_interaction_only_breakdown`がコンポーネントをフィルタする際に除外される
  ため、ペナルティ計算の土俵に乗らない)
- `protein_hunter_score`にも適用されない(`build_candidate_score`は
  BLAST陽性/陰性ヒットの有無だけを見て、`negative_hit_strength`の強度
  ラベル自体は使わない)

つまり**現状、`negative_hit_strength`ペナルティを実際に反映しているのは
`interaction_priority_score`だけ**です。Final Scoreを`protein_hunter_score`
+`interaction_score`から素朴に合成すると、このペナルティは自動的には
一切引き継がれません(design specが求める3要素目「negative_hitペナルティ」
を明示的に足し戻す設計が必要)。

### ご質問への直接の回答:`negative_hit`系バケツ無効時でもペナルティは効いているか

**「候補として出力されるかどうか」と「ペナルティが計算に反映されるかどうか」
は別の話です。**

- ペナルティの計算自体(`candidate.negative_hit_strength`を読む処理)は、
  その候補が**どれかのcandidate_sourcesバケツで実際にスコアリングされてさえ
  いれば**、`negative_hit`バケツ自体が無効でも適用されます。例:
  `Candidates_relaxed`は既定で有効で、medium/weak負ヒットを許容する候補
  (`negative_hit_strength = medium`または`weak`)を含み得るため、そうした
  候補の`interaction_priority_score`にはペナルティがちゃんと反映されます。
- 一方、`negative_hit_strength = strong`の候補は`Candidates`/`Candidates_relaxed`
  どちらの分類にも入らず、**専用の`Negative_hit`系バケツでしか候補プールに
  現れません**。このバケツが既定で無効(PR #11で確認済み)なので、こうした
  候補は「ペナルティが効いていない」のではなく、**そもそもどの`Interaction_*`
  シートにも一切スコアリングされず出現しません**——ペナルティ云々以前の話です。

## 4. `Evidence_Tier`の計算方法・閾値(`Confidence`列は存在せず)

`analysis/scoring_engine.py::_classify_tier`(186-197行目)、既定値は
`TierThresholds`(`analysis/scoring_engine_config.py`):

| Tier | 条件 |
|---|---|
| `Tier1_VeryStrong` | `final_score >= 70` かつ `evidence_category_count >= 3` |
| `Tier2_Strong` | `final_score >= 50` かつ `evidence_category_count >= 2` |
| `Tier3_Moderate` | `final_score >= 25` かつ `evidence_category_count >= 1` |
| `Tier4_Weak` | 上記いずれにも該当しない(ただしeligible) |
| `Unclassified` | `final_score`が`None`(証拠不足でスコア自体不成立、
  `minimum_evidence.min_categories`/`min_available_weight`未達) |

`interaction_evidence_tier`列も同じ関数・同じ閾値を、`interaction_priority_score`
ではなく`interaction_score`側の`final_score`/`category_count`に対して適用
しています(`_interaction_only_breakdown`)。**リポジトリ内に別途
"Confidence"という名前のフィールド・列は存在しません**——Tierが唯一の
信頼度概念です。

## 5. `scoring_engine.py`の再利用可否

**技術的には再利用可能で、自然な適合です。** `EvidenceComponent`/
`score_candidate`/`CategoryScore`はドメイン非依存の汎用実装(モジュール
docstring自身が「evidence categories, caps」の一般モデルと説明)——
"protein_hunter_score"と"interaction_score"を**それぞれ1コンポーネントだけ
を持つ2つの新しいトップレベルカテゴリ**として扱うことは、既存コードを
一切変更せずに可能です。具体的な形:

```
components = [
    EvidenceComponent.available(
        "protein_hunter_score", "candidate_quality",
        normalized_value=protein_hunter_score / PHS_CEILING,  # 要決定、§7参照
        weight=..., raw_value=protein_hunter_score,
    ),
    EvidenceComponent.available(
        "interaction_score", "query_specific_evidence",
        normalized_value=interaction_score / 100.0,
        weight=..., raw_value=interaction_score,
    ) if interaction_score is not None else EvidenceComponent.unavailable(
        "interaction_score", "query_specific_evidence", EvidenceStatus.MISSING,
    ),
    # 任意: negative_hit_strengthを再度is_negative=Trueとして投入し、
    # トップレベルでもペナルティを利かせる(§3で見つけた抜け漏れの解消)
]
final_breakdown = score_candidate(components, final_score_engine_config)
```

`score_candidate`はMISSING(`interaction_score`がNoneの候補=クエリなし文脈
での`Candidates`シート単体表示など)を自動的に分母から除外するため、
「クエリ文脈が無ければFinal Scoreは`protein_hunter_score`だけで算出」という
design specの"missingはゼロではない"哲学とも自然に整合します。Tier分類
(`_classify_tier`)もそのまま流用でき、新しい"Final_Tier"を無料で得られます。

### 懸念点

1. **正規化係数の恣意性**: `protein_hunter_score / PHS_CEILING`の
   `PHS_CEILING`をいくつにするか(理論値18か、実質的な最大値14か、経験的な
   分位点か)自体が設計判断で、後述§6の通り結果を大きく左右します。
2. **二重のcap正規化**: `interaction_score`自体が既に1回`score_candidate`
   で正規化された出力です。それをもう一度「1コンポーネントとしてcapをかける」
   のは技術的には問題ありませんが、"capの中にcapされた値が入っている"形に
   なるため、監査時にやや分かりにくくなる可能性があります(実害はない)。
3. 上記を許容できない場合の**代替案(単純な正規化後の加重平均)**:
   ```
   final_score = clamp(
       w_phs * (protein_hunter_score / PHS_CEILING * 100)
       + w_int * interaction_score
       - negative_penalty_points,
       0, 100,
   )
   ```
   `scoring_engine.py`を経由しないため実装は数行で済みますが、MISSING処理
   ・Tier分類・監査ログ(`ScoreBreakdown`)を自前で再実装する必要があり、
   このプロジェクトが既に持つ「MISSING≠ゼロ」の哲学を保つには手動での
   注意が要ります。

**結論**: 再利用(前者)を推奨します。既に何度もテストされた分母正規化・
MISSING処理・Tier分類ロジックをそのまま使い回せるメリットが、二重cap構造の
分かりにくさという小さなデメリットを上回ります。

## 6. 実データでの簡易試算(Tier A正例8ペア vs AlphaFold3陰性28件)

同一の直近診断run(`data/output/ProteinHunter_results_calibration_check.xlsx`、
`claude/experimental_interactions_calibration_report_pairs.csv`の
`tier_final=="A"`、8ペア)から`protein_hunter_score`・`interaction_score`を
実際に取得し、3パターンを試算しました(`PHS_CEILING=18`で0-100に正規化):

| 指標 | 正例(n=8) mean(median) | 陰性(n=28) mean(median) | 分離幅(正例-陰性) |
|---|---:|---:|---:|
| `protein_hunter_score`(0-100正規化のみ) | 50.00 (50.00) | 56.94 (50.00) | **-6.94**(逆転) |
| `interaction_score`(単独) | 39.55 (43.85) | 12.95 (12.27) | **+26.61** |
| 50%/50%ブレンド | 44.78 (46.92) | 34.95 (32.33) | +9.83 |
| 30%(phs)/70%(interaction)ブレンド | 42.69 (45.69) | 26.15 (25.26) | +16.54 |

**`protein_hunter_score`を混ぜるほど分離幅が縮みます**——50/50では
`interaction_score`単独の分離幅(+26.61)の1/3強(+9.83)まで劣化、
30/70でも約6割(+16.54)にとどまります。原因は明確: Tier A正例8ペアは
**全件が`Negative_hit`バケツに分類され`protein_hunter_score`が完全に
9(=50.00)で一致**しており(§1の`Candidates`分布とも整合)、
`protein_hunter_score`側に群内分散がゼロ——`interaction_score`が持つ本物の
識別力を、無関係な定数値で薄めているだけの状態です。

**重要な留保**: n=8・正例側の分散ゼロという極端なサンプルなので、
「`protein_hunter_score`は本質的に無価値」という結論は導けません。むしろ、
正例8ペア全件が偶然`Negative_hit`バケツに集中しているという
**バケツ構成の偏り**(PR #11で既出の発見と同根の問題)が主因である可能性が
高いです——真に多様な候補源(Candidates本体も含む)からの正例セットで
再検証する価値があります。詳細CSV:
`.cache/geo_investigation/final_score_sim_positive.csv` /
`final_score_sim_negative.csv`(未コミット、`.cache/`は既にgitignore対象)。

## 設計提案

### 合成方式オプション

| # | 方式 | 概要 | 長所 | 短所 |
|---|---|---|---|---|
| A | **scoring_engine.py再利用(推奨)** | §5の通り2カテゴリ(+任意でnegative_hit_strengthカテゴリ)としてラップし`score_candidate`に通す | 既存のMISSING処理・Tier分類・監査trailをそのまま獲得。`Interaction_Evidence_Detail`と同じ形式で"Final_Score_Detail"も同じコードパスで自然に作れる | 二重cap構造(§5の懸念2) |
| B | 単純加重平均+手動MISSING処理 | §5の代替案の式をそのまま実装 | 実装が単純・追跡しやすい | MISSING/Tier/監査を独自実装する必要あり、将来的にscoring_engine.pyとロジックが乖離するリスク |
| C | `interaction_score`優先・`protein_hunter_score`はタイブレークのみ | Final Score = `interaction_score`本体、`protein_hunter_score`は同点時の第2ソートキーとしてのみ使用(スコアには算入しない) | §6の実データが示す「混ぜると悪化」を最も忠実に反映。coexpression_gse77738を`interaction_score`から除外した先例(PR #11)と一貫した判断 | design specが求める「1つの数値へ統合」という文言との整合性が薄い(統合はするが実質的に一方が支配的) |

A・Cは併存できます(AのFinal Score計算自体で`interaction_score`側の
重み/capを`protein_hunter_score`側より大きく設定すれば、§6の实データに
合わせた挙動をAの枠組みの中で実現可能)。**A(実装方式)+ Cに近い重み付け
(内容)** の組み合わせを暫定推奨としますが、正確な重み/capは§6のn=8という
小標本ゆえ、実装後にキャリブレーションが必須です(design spec §37と同じ
"暫定値、後で校正"のポリシーを踏襲)。

### 段階分割案(M1...、未実装)

- **M1**: Final Score用の新カテゴリ定義(`candidate_quality`
  [=protein_hunter_score]、`query_specific_evidence`[=interaction_score]、
  任意で`negative_hit_strength`の再カテゴリ化)と、`PHS_CEILING`含む正規化
  パラメータを`scoring_engine_config.py`的な別ファイルに切り出し(暫定値、
  §37ポリシー通りコメントで明記)。
- **M2**: `analysis/scoring_engine.py::score_candidate`を呼び出す
  Final Score計算関数の実装(方式A)。`Interaction_*`シートに
  `final_score`・`final_score_tier`列を追加(既存列は一切変更しない、
  追加のみ)。
- **M3**: `Interaction_Evidence_Detail`と同型の"Final_Score_Detail"サブ
  シート、またはEvidence_Detailへの追加行(protein_hunter_score/
  interaction_score/negative_hit_strengthの3コンポーネントの内訳)。
- **M4**: Tier A/AF3陰性データでの本格的な再検証(§6は簡易試算のみ)、
  重み/capの初回キャリブレーション。
- **M5**: `ranking_metric`に`final_score`を追加選択できるようにする
  (既存の`interaction_priority_score`/`interaction_score`との併存)。

### 未決事項(ご判断をお願いします)

1. **`PHS_CEILING`をいくつにするか**:理論値18、実質的な観測最大値14、
   あるいは経験的な分位点(例:観測分布の95パーセンタイル)のいずれか。
   §6の通りこの選択がFinal Scoreの実質的な意味を大きく左右します。
2. **合成方式A/B/Cのどれを採用するか**(併用可、上記参照)。
3. **`negative_hit_strength`ペナルティをFinal Scoreにも明示的に組み込むか**
   (§3の抜け漏れの解消)。組み込む場合、`interaction_priority_score`側の
   既存ペナルティと二重適用にならないよう、Final Scoreは
   `protein_hunter_score`+`interaction_score`の「素の値」から再計算し、
   `interaction_priority_score`(ペナルティ適用済み)は使わない、という
   設計になります——明示的な確認をお願いします。
4. **クエリ文脈が無い場合(`interaction_score`がそもそも存在しない
   `Candidates`シート単体表示など)にFinal Scoreをどう扱うか**:
   `protein_hunter_score`のみで計算するか、Final Score自体を空欄にするか。
5. **§6の結果を受けて、`protein_hunter_score`の重み付けをどの程度小さく
   すべきか**(あるいは方式Cのようにスコアに算入せずタイブレークのみに
   格下げすべきか)。n=8の暫定結果である点を踏まえ、本実装前に正例データを
   増やす価値があるかも含めてご判断ください。

## 実装後の実データ検証:重要な発見

上記の承認方針(cap 30/70、negative_hitペナルティを独立適用)通りに実装し、
同一の診断run条件(Tier A正例8ペア、AlphaFold3陰性28件、同一パイプライン実行)
で`interaction_score`単独 vs `final_score`を比較しました。生データ:
`claude/final_score_verification_positive.csv` /
`final_score_verification_negative.csv`。

| 指標 | 正例(n=8) mean/median | 陰性(n=28) mean/median | 分離幅(正例-陰性) |
|---|---:|---:|---:|
| `interaction_score`単独 | 39.83 / 41.22 | 12.28 / 12.59 | **+27.55** |
| `final_score`(実装後) | 16.05 / 13.86 | 16.00 / 17.94 | **+0.05(ほぼ消失)** |

**`interaction_score`単独では極めて強い分離(+27.55)があったにもかかわらず、
`final_score`ではその分離がほぼ完全に消失しました。** 原因は`protein_hunter_score`
の希釈効果(§6で既に確認済み、小さい)ではなく、**`negative_hit_strength`
ペナルティの独立適用でした**。`Interaction_Evidence_Detail`で直接確認したところ、
**Tier A正例8ペアの候補側は全件`negative_hit_strength = "strong"`**
(`normalized_value=1.0`、`weight=30.0`——`negative_penalty_cap`の上限をちょうど
使い切る)でした。これは偶然ではなく構造的な理由です: Tier A正例は
(先の統合レポートで既に発見した通り)全件が`Negative_hit`/`Negative_strong_hit`
バケツに分類されており、**そのバケツに分類される基準そのものが
`negative_hit_strength = strong`であるため**、正例側は例外なくペナルティを
フルに受けます。一方、AlphaFold3陰性側は`Candidates`/`Candidates_relaxed`
バケツ由来の候補も多く含み、これらは陰性ヒット自体を持たない(ペナルティ
ゼロ)ケースが相当数あります。結果として、ペナルティが**保存性の高い
真の相互作用パートナーを、多くの陰性候補より重く罰する**という、意図とは
逆方向の効果を生んでいます。

この結果は「実装が間違っている」のではなく、**承認済み設計(negative_hit
ペナルティの独立適用)を実データに正確に適用した結果、想定より大きな副作用が
見つかった**ということです。cap/weight自体は今回変更していません
(ご指示通り30/70のまま、ペナルティも既存実装を参考にした値のまま)。
判断が必要な点として報告します——例えば、Final Score用のペナルティ重み/cap
を`interaction_priority_score`側と別に(より小さく)設定する、あるいは
Final Scoreからはペナルティを除外する、といった選択肢が考えられますが、
実装方針の変更はご判断を仰いでから行います。

## negative_hit_strengthペナルティの除外(方針転換の決定)

上記の発見を受けて、以下の方針転換が決定されました:

**`negative_hit_strength`は本来、系統特異性/新規性のシグナル(「このタンパク質
は陰性参照ゲノムにも広く存在するありふれたものか」)であり、design spec §7.7が
定義する「Negative Evidence」(functional contradiction、incompatible
localization、incompatible domain、phylogenetic contradiction——「このペア
自体が相互作用として矛盾している」という反証)とは概念として別物です。**
Tier A正例8ペア(Hdr/Mtp/Nif複合体)が軒並み`negative_hit_strength = strong`
なのは、これらが古くから保存された中心代謝系だからであり、相互作用の妥当性とは
無関係です。この2つの概念を混同してFinal Scoreにペナルティとして組み込んだ
ことが、真の相互作用ペアを一律に減点する結果を招きました
(PR #11で発見した「negative_hitバケツの盲点」と同根の問題が、スコア計算の
場面で再発したもの)。

対応: Final Scoreから`negative_hit_strength`ペナルティを完全に除外。
`final_score_negative_penalty`コンポーネント自体は監査列の枠組みとして残すが、
値は常に`NOT_APPLICABLE`(将来、design spec §7.7が本来意図する真の生物学的
矛盾シグナル——機能矛盾・局在不整合など——が実装された際に使う予約枠であり、
`negative_hit_strength`はその代用にはならない)。`interaction_priority_score`
側の既存の`negative_hit_strength`適用は変更なし(そちらは候補全体の妥当性を
評価するものであり、Final Scoreとは別の設計判断として妥当)。cap配分
(protein_hunter=30/interaction=70)は今回の問題と無関係のため変更していない。

## 修正後の再検証

同一条件(Tier A正例8ペア・AlphaFold3陰性28件・同一パイプライン実行)で
再検証した結果、分離幅が想定通り回復しました:

| 指標 | 正例(n=8) mean/median | 陰性(n=28) mean/median | 分離幅(正例-陰性) |
|---|---:|---:|---:|
| `interaction_score`単独 | 39.83 / 41.22 | 12.28 / 12.59 | +27.55 |
| `final_score`(ペナルティ適用時、修正前) | 16.05 / 13.86 | 16.00 / 17.94 | +0.05 |
| `final_score`(ペナルティ除外後、修正後) | **42.88 / 43.86** | **25.68 / 25.23** | **+17.20** |

`interaction_score`単独の分離幅+27.55に対し、`final_score`は+17.20と、
`protein_hunter_score`混合による希釈分だけやや下がる程度に収まりました
(実装前の簡易試算§6の30/70ブレンド予測値+16.54ともほぼ一致)。ご想定通りの
挙動に回復したことを確認しました。生データ:
`claude/final_score_verification_positive.csv` /
`final_score_verification_negative.csv`(修正後の値に更新済み)。
