# 校正データ拡充パイプライン — 設計書

Status: 設計提案(未実装)。`main`(PR #15まで)および既存の校正関連
成果物(`claude/experimental_interactions_curated.csv`,
`claude/experimental_interactions_calibration_report*.csv/.md`)を
読んで現状の作業実態を確認済み。**2026-09-20: [[run_provenance_design]]
がPR #21として実装済み(未マージ)になったため、下記「実行結果Excel」の
項をサイドカーYAML前提の記述に更新済み。実装順序上の依存はなし
(このツールは新規追加のみで既存コードに触らないため、PR #21のマージを
待たずに実装着手してよい。ただしサイドカーファイル名の規約は
`<output_excelのstem>.run_provenance.yaml`で確定しているので、実装時は
それに合わせること)。**

## 背景・問題

PR #9〜#11のキャリブレーション(design spec §37対応)は、以下の点で
**再利用不可能な一回限りの手作業**として行われた。

- `claude/experimental_interactions_calibration_report.md` 自身が明記して
  いる通り、使用したExcel(`ProteinHunter_results_calibration_check.xlsx`)
  は**リポジトリにコミットされていない**。「configから再生成可能」とは
  書いてあるが、再生成して`_pairs.csv`/`_negatives.csv`を作り直す
  スクリプトは存在せず、当時の突き合わせ作業(curated CSV × Excel出力 ×
  AlphaFold3陰性リスト)は再現不能な手順で行われた。
- 統計的な検証も「平均・中央値を並べて目視で比較する」段階に留まっており
  (同ドキュメントの「Reading the two tables together」節が自ら
  "not conclusions -- 判断は保留" と明記)、n=8という小標本に対する
  有意性検定や信頼区間の計算はしていない。
- `claude/implementation_status_scoring_v2.md` の「進行中・未着手」に
  ある通り、Tier B・要確認19ペアは「論文特定やパラログクラスタの解消が
  進めば追加の校正データとして使える可能性がある」と書かれているが、
  それを実際に追加するための受け皿(スクリプト・スキーマ)が無いため、
  着手のハードルが実質的に高いまま放置されている。
- **2026-09-20に追加で判明した実例**: [[negative_reference_species_mislabeling_finding]]
  で行った前後比較(config.yaml既定設定・キャリブレーション相当設定)も、
  今回の設計が目指す「決定的に再生成できるCLI」ではなく、都度Excelを
  読み直して手計算する形で行われた。まさに本設計が解決しようとしている
  問題そのものが、この直近の調査でも繰り返し発生していたことになる。

**今のままでは、校正データを1件追加するたびに今回と同じ量の手作業
(Excel出力を目で追って突き合わせ、Pythonで都度統計を計算)が発生する。**
これは「拡充」を事実上妨げている。

## 採用する方針

一回限りの分析を、3つの入力ファイルから決定的に再生成できる
**単一のCLIスクリプト**に切り出す。design spec自体を変更するものではなく、
既存の `output/excel.py` が書き出す列(`Interaction_Evidence_Detail`,
`02_Final_Score`)をそのまま読むだけの、**分析専用・スコアリングには
一切影響しない**ツールとして位置づける(既存コードのポリシー
「既存の計算結果を変えない追加ツール」— Wordレポート機能と同じ立ち位置)。

### 入力(3種、いずれもCSV/YAMLで人間がGit管理できる形)

1. **正例キュレーションCSV**: 既存の
   `claude/experimental_interactions_curated.csv` と同じスキーマ
   (`protein_a_old_locus_tag, protein_a_label, protein_b_old_locus_tag,
   protein_b_label, class, confidence, source`)。Tier分類列
   (`tier: A/B/excluded`)を1列追加して、今回別文書だった
   「Tier分け」もこのCSV1本に統合する(現状 curation.md の文章内に
   埋め込まれているTier A/B/除外の判断を、機械可読な列として持たせる)。
2. **陰性リストCSV**: 既存の `alphafold3_calibration_comparison.csv`
   相当のスキーマ(`old_locus_tag, af3_classification, af3_ipTM`)。
   AlphaFold3で確認した陰性以外にも、将来別の方法(実験的陰性等)で
   確認した陰性ペアを同じスキーマで追加できるよう、
   `source`列(例: "alphafold3", "experimental_negative")を持たせる。
3. **実行結果Excel**: 通常の `ProteinHunter_results_*.xlsx`
   (`Interaction_Evidence_Detail` + `02_Final_Score` シート)。
   **[[run_provenance_design]]がPR #21で実装済みのため、`--results`で
   渡されたExcelと同じディレクトリに`<stemと同名>.run_provenance.yaml`
   が存在すれば読み込み、`config_hash`/`git_commit`/`git_dirty`/
   `reference_genome_check_passed`を`calibration_summary.md`に転記する。
   サイドカーが無い場合(古い出力、または`--results`に手作業で用意した
   Excelを渡した場合)はこの節を「provenance情報なし」として省略し、
   処理は続行する(必須依存にしない)。** これにより「この設定で計測」
   という事故(Stage 1検証・PR #20前後比較で実際に起きた設定取り違え)
   の再発をレポート側からも検知しやすくする。

### 出力

- `<curated>_pairs.csv` / `<negatives>_matched.csv`
  (既存の `_pairs.csv`/`_negatives.csv` と同スキーマの後継、
  自動生成に置き換え)。
- `calibration_summary.md`(Markdown、既存の
  `experimental_interactions_calibration_report.md` の「Reading the
  two tables together」相当のセクションを自動生成): 各スコア列
  (`interaction_score`, `interaction_priority_score`,
  `string_neighborhood`, `string_cooccurrence`,
  `coexpression_gse77738`, `coexpression_gse64349`, `final_score`)ごとに、

  - Tier A/B別・陰性別の平均・中央値(現状通り)
  - **追加: Mann-Whitney U検定のp値**(scipyへの依存が増えるが、
    `pandas`/`numpy`は既に依存済みなので `scipy` 追加は妥当な範囲。
    n=8のような小標本でも計算自体は可能で、「有意差なし」という
    結果も含めて正直に報告できる)
  - **追加: ROC-AUC**(`interaction_score` がどれだけ正例/陰性を
    順位で分離できるかを1指標にする。`sklearn` はさすがに重いので、
    自前実装 or `scipy.stats.mannwhitneyu` のU統計量から
    `AUC = U / (n1 * n2)` の関係式で算出すれば追加依存なしで済む
    — この方が既存依存に近く望ましい)
  - 実行に使ったExcelのprovenance情報(サイドカーYAMLがあれば
    config_hash/git commit/reference_genome_check_passedを記載)

### CLIの形

```
python tools/calibration_report.py \
  --curated claude/experimental_interactions_curated.csv \
  --negatives claude/alphafold3_calibration_comparison.csv \
  --results data/output/ProteinHunter_results_calibration_check.xlsx \
  --out claude/calibration/2026-09-19_run \
```

`--out` ディレクトリに `pairs.csv`/`negatives_matched.csv`/
`summary.md` の3点セットを書き出す。日付付きディレクトリにすることで、
過去の校正結果を上書きせず時系列で残せる(「校正データの拡充」は
1回で終わる作業ではなく継続的に積み増すものなので、履歴を残す設計に
しておく)。

### 既存資産の扱い

- `claude/experimental_interactions_curated.csv` は
  `tier` 列を追加した上でそのまま入力として使う(過去の手作業の結果を
  無駄にしない、後方互換)。
- `claude/experimental_interactions_calibration_report.md` は
  「このツール登場以前の、最初の一回限りの分析の記録」として
  そのままリポジトリに残す(過去の意思決定の記録としての価値は
  ドキュメントに書かれた考察・注記(HdrD1-Merの外れ値の考察など)に
  あり、これはツールが自動生成する数値だけでは再現できないので、
  削除や置き換えはしない)。冒頭の訂正ブロック(PR #17で追加、
  [[negative_reference_species_mislabeling_finding]]参照)もそのまま
  残す。

### 校正データの拡充そのものについて(ツールとは別に必要な作業)

ツールが整えば、以下は「入力CSVに行を足すだけ」で着手できるようになる。

- Tier B・要確認19ペア: 2026年ACDS論文の特定、CdhC/CdhDパラログ
  クラスタの解消(どちらのパラログが実際に相互作用するかのGFF上の
  確認)が進み次第、`curated.csv` の該当行の `tier` を `B`→`A` に
  昇格するだけで次回のレポート生成に反映される。
- AlphaFold3陰性セット(現状28件、MA_4115クエリ限定)を、他のTier A
  クエリ(HdrD1, MtpA, NifD/K)についても同様にAlphaFold3で陰性確認
  した候補を追加できると、「陰性側もクエリ非依存に均質化する」ことが
  でき、現状の「陰性はMA_4115限定・正例はMA_4115以外」という
  非対称性(`implementation_status_scoring_v2.md` が明記する制約)を
  緩和できる。これは実験(AlphaFold3実行)そのものはユーザー側の作業だが、
  結果を受け取るCSVスキーマは今回のツールと共通にしておく。

## スコープ外(あえてやらないこと)

- キャップ/重みの自動最適化(ロジスティック回帰等でのフィッティング)は
  今回はやらない。n=8〜数十件規模でパラメータをデータにフィットさせると
  過学習のリスクが高く、「暫定値のまま」という現状の誠実な運用の方が
  この標本サイズでは妥当。統計指標(Mann-Whitney/ROC-AUC)は
  「判断材料を増やす」ためであり、「自動でパラメータを変える」
  ためではない。
- ホールドアウト/交差検証の枠組みも、正例データ数が2桁に達するまでは
  時期尚早と判断し、今回のスコープには含めない。

## テスト計画

- `tests/test_calibration_report.py`(新規): 小さな合成データ
  (数行のcurated CSV・negatives CSV・ダミーのExcel出力)を使い、
  出力される `pairs.csv`/`summary.md` の値を手計算した期待値と突き合わせる。
  特にMann-Whitney U統計量とAUC換算式の一致を確認するテストを含める。
  サイドカーYAMLがある場合/ない場合の両方でエラーなく動作することの
  確認も含める。
- 実データ検証: 今回の設計で、既存の
  `claude/experimental_interactions_calibration_report_pairs.csv` と
  同じ入力データセットを流し、`interaction_score`の平均値(39.55/12.95)
  など既存レポート記載の数値と一致することを確認する(回帰確認)。

## 見積もり

新規スクリプト1本(`tools/calibration_report.py`)+ 新規テスト。
`scipy` の依存追加が必要かどうかは実装時に再検討(AUC計算をU統計量経由の
自前実装にすれば `scipy` すら不要にできるため、依存を増やしたくなければ
その方式を選ぶ)。既存のスコアリング・出力コードには一切手を入れない
独立ツールなので、既存テストへの影響はゼロ。1PRで完結できる。