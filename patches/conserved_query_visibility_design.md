# クエリ保存性の自動検知と警告 — 設計書

Status: 実装済み(PR #16)。ユニットテストと結合テストに加え、実データ(5クエリ)でも
警告の出力(ログ・Word 5.7)を確認済み。ただし陰性参照ゲノムの取り違えを修正した後は
この5クエリのいずれも strong/medium にならず警告は出ない
(claude/negative_reference_genome_mixup_investigation.md を参照)。実装前にコードベース
(`main`, PR #15まで)を読んでintegration pointを確認済み。

実装時の判断: 警告は `negative_hit` とその3サブバケツ(strong/medium/weak)が
**すべて無効**のときだけ出す(いずれか有効ならユーザーが既に判断済みとみなして黙る)。
マッチする target レコードのない配列のみのクエリは保存性が不明なので黙る。

## 背景・問題

`claude/experimental_interactions_calibration_report.md`(PR #9〜#11)で
判明した「発見2」は、単なるキャリブレーション不足ではなく構造的な死角である。

HdrD1-Mer、MtpA-MtpC、NifD/K-NifI1/I2 といった Tier A の実証済み正例ペアは、
古くから保存された中心代謝系であるため、陰性参照ゲノムにも強くBLASTヒットし、
`candidate_source = Negative_hit`(またはそのサブバケツ)に分類される。
`INTERACTION_CANDIDATE_SOURCE_DEFAULTS`(`analysis/interaction_scoring.py`)は
`negative_hit: False` かつ3つのサブバケツも全て `False` がデフォルトなので、
**これらの真陽性ペアはデフォルト設定では候補プールにすら入らず、
スコアの良し悪し以前にExcel/Wordレポートのどこにも一切現れない。**

現状の対処は `config.yaml` のコメントで「クエリが保存性の高い経路に関わる
場合は明示的に有効化を検討すべき」と注意喚起するだけ(PR #11)。これはクエリ
ごとに手動で気づいて判断する必要があり、見落としやすい。

`consider_cross_species_matches`(`config.py`)という既存のプリセット
スイッチは `False` にすると `negative_hit` を含む3バケツを一括で有効化するが、
これは以下の理由で今回の問題に対する的確な解決策ではない:

- `max_candidates_per_query` が 200→500 に、`ranking_metric` が
  `interaction_priority_score`→`interaction_score` に、
  `annotation_targets(gff)` が全シートtrueに、と**無関係な設定まで
  まとめて変わる**全体スイッチである。
- 複数クエリを1回の実行で扱える現在のパイプライン(Phase 6-8 Stage 2で
  複数クエリ対応済み)において、「あるクエリだけ保存性が高い」という
  クエリ単位の性質に対して、実行全体を一律に broaden させるのは過剰。
- 保存性が低い他のクエリにとっては、このスイッチを入れることで無関係な
  真陰性候補まで大量に候補プールへ入り、ノイズが増える(まさに
  `negative_hit` バケツ本来の設計目的=系統特異的な新規候補探索、を損なう)。

## 採用する方針

**クエリ自身の `negative_hit_strength` を実行時に見て、強い保存性ヒットが
あれば当該クエリについてだけ警告を出す。** 候補プールや既定のスコアリング
挙動は一切変更しない、純粋な診断機能として実装する。

### この方針が低コストで実装できる理由(既存コードの確認結果)

`analysis/blast_pipeline.py` の `run_blast_classification_pipeline` は
`query_fasta=target_fasta`、つまり **ゲノム全体(target.faa)の全タンパク質を
陰性参照ゲノムに対してBLASTし**、`populate_negative_hit_evidence(records, ...)`
(`analysis/ortholog_filter.py`)で `records` 内の**全レコード**に
`negative_hit_strength`/`negative_strong_hit_count` 等を設定する。
クエリタンパク質も `target.faa` 内の1レコードに過ぎないため、
**新たなBLAST実行や新規計算は一切不要で、既に計算済みの値を読むだけで良い。**

`analysis/interaction_scoring.py::run_interaction_scoring` は既に
`_resolve_query(query, index, blast_classification.all_records)` で
各クエリをレコードに解決しており、`warnings: list[str]` という
`InteractionScoringResult` に載る警告リストの仕組みも既に存在する
(`config.redundant_negative_hit_sources()` と同型の「実行時警告、
処理は継続」パターン)。

### 実装内容

1. `analysis/interaction_scoring.py` に
   `_check_conserved_query_visibility(resolved_queries, records, candidate_sources, warnings)`
   のようなヘルパーを追加(新規関数、既存関数は変更しない)。
   - `_resolve_query` の返り値(またはその中の `matched_record`)から
     `negative_hit_strength` を取得。関数の返り値dictに
     `negative_hit_strength` フィールドを1つ追加するだけで済む
     (既存フィールドは変更しない、追加のみ)。
   - `negative_hit_strength in {"strong", "medium"}` かつ
     `candidate_sources.get("negative_hit", False)` が `False`
     (かつ該当サブバケツも無効)の場合にのみ警告を1件生成。
   - `run_interaction_scoring` 内、`resolved_queries` を作った直後
     (line ~577 `query_rows = [...]` の前後)で呼び出し、
     返ってきた警告文字列を既存の `warnings` リストに追加。
2. 警告文言の例(実データの言い回しに寄せる):
   > "Query {query_id} itself has a {strength} negative-reference BLAST
   > hit (negative_hit_strength={strength}). Experimentally confirmed true
   > interaction partners for highly conserved queries (see Tier A in
   > claude/experimental_interactions_calibration_report.md — HdrD1-Mer,
   > MtpA-MtpC, NifD/K-NifI1/NifI2) are known to fall into the Negative_hit
   > bucket for the same reason. With the current
   > interaction_scoring.candidate_sources settings, such candidates are
   > silently excluded from every output sheet regardless of their score.
   > Consider re-running this query with candidate_sources.negative_hit
   > (and/or its strong/medium/weak sub-bucket) enabled."
3. `main.py` 側は変更不要(既存の warnings 集約・ログ出力経路
   — Excel/Wordどちらの生成コードパスでも `interaction_result.warnings`
   は既にログに出す実装があるはず。無ければ、既存の
   `redundant_negative_hit_sources` 警告のログ出力箇所と同じ場所に
   `for w in interaction_result.warnings: logger.warning(w)` を1回追加する
   だけで済む — 既存経路の確認はimplementation時にmain.pyを再読すること)。
4. **Word/Excelレポートへの掲載(ログだけでは見落とすため)**:
   `output/word_report.py::write_word_report` の
   "5.7 Negative Evidence (reserved)" セクション
   (`_write_evidence_architecture` 内、line ~304)に、
   該当クエリがあれば1文追記する。ここは design spec準拠でreservedと
   されている箇所だが、「reservedのまま」と「実際に起きている既知の
   可視性問題を明記する」は矛盾しない — Negative Evidence機能そのものは
   未実装のままで良いが、この警告はその不在が引き起こす副作用の説明として
   同じセクションに収まる。
   `write_word_report` は既に `blast_classification` と
   `interaction_result` を受け取っているので、
   `interaction_result.warnings` からこの種の警告文字列だけを
   フィルタして段落を追加すればよい(接頭辞などで種類を識別できるように、
   警告文字列を組み立てる関数に `kind="conserved_query_visibility"` の
   ような軽いタグ付けをしておくと、ログ用テキストとWord用テキストを
   同じ関数から作れて重複を避けられる)。

### スコープ外(あえてやらないこと)

- **候補プールを自動で広げることはしない。** 警告を出すだけで、
  `candidate_sources` を実行時に自動変更するような「賢い」挙動は入れない。
  ユーザーが意図的に `negative_hit` を有効にするかどうかを判断できるよう、
  常に人間の判断を経由させる(design spec §35 の「断定禁止」の精神と同じ、
  パイプラインが勝手に強い判断をしない、という既存の設計哲学に合わせる)。
- **「Negative_hitに分類されたが実は高スコアな候補」を常に事前計算して
  見せる**、というより踏み込んだ案(進行中の会話で検討した案)は今回は
  採用しない。理由: `negative_hit` バケツ自体のスコアリングを
  `candidate_sources` の設定に関わらず常時計算する必要があり、実行時間
  (CDD annotationが既に974秒かかっている)への影響と実装コストが、
  「まず気づけるようにする」という今回のゴールに対して過大。将来、
  今回の警告だけでは不十分と分かった場合の拡張案として記録だけしておく。

## テスト計画

- `tests/test_interaction_scoring.py` に新規テスト:
  - クエリの `negative_hit_strength` が `strong`/`medium`/`weak`/`none`
    それぞれの場合で警告の有無を確認(`weak`/`none`では警告なし)。
  - `candidate_sources.negative_hit` が既に `True` の場合は警告が出ない
    ことを確認(「既に対処済みなら黙る」)。
  - 複数クエリのうち1つだけ保存性が高いケースで、警告文にそのクエリの
    IDだけが含まれることを確認(他クエリを巻き込まない)。
- 実データ検証: MA_4115(保存性が低いはずの本命クエリ)では警告が出ず、
  Tier A由来のクエリ(HdrD1=MA_0688など)を試験的に `query_proteins` へ
  加えて実行し、警告が実際に出ることを確認する。

## 見積もり

コード変更は新規関数1つ+既存2箇所への数行の呼び出し追加、新規計算なし
(既存データの参照のみ)なので、実装規模としては
`redundant_negative_hit_sources()`(PR #11)と同程度〜やや大きい程度。
1PRで完結できる。
