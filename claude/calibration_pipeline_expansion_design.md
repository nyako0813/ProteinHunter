# 校正データ拡充パイプライン — 設計書

Status: **実装完了・実データ検証済み。PR #22(本体)・PR #23(出力形式の
Excel化)としてpush済み(いずれも未マージ)。** 全634テストパス。設計から
の差分は「実装結果(PR #22)」「実装結果(PR #23)」参照。

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
- 2026-09-20に追加で判明した実例: [[negative_reference_species_mislabeling_finding]]
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
既存の `output/excel.py` が書き出す列をそのまま読むだけの、**分析専用・
スコアリングには一切影響しない**ツールとして位置づける(既存コードの
ポリシー「既存の計算結果を変えない追加ツール」— Wordレポート機能と
同じ立ち位置)。

### 入力(3種、いずれもCSV/YAMLで人間がGit管理できる形)

1. **正例キュレーションCSV**: 既存の
   `claude/experimental_interactions_curated.csv` と同じスキーマ
   (`protein_a_old_locus_tag, protein_a_label, protein_b_old_locus_tag,
   protein_b_label, class, confidence, source`)。Tier分類列
   (`tier: A/B/excluded`)を1列追加して、今回別文書だった
   「Tier分け」もこのCSV1本に統合する。**入力CSVはExcel化の対象外**
   (git管理でdiffが読める形を維持するため、ユーザーの明示的な指示)。
2. **陰性リストCSV**: 既存の `alphafold3_calibration_comparison.csv`
   相当のスキーマ(`old_locus_tag, af3_classification, af3_ipTM`)。
   `source`列(例: "alphafold3", "experimental_negative")を持たせる。
   **これも入力CSVのためExcel化の対象外。**
3. **実行結果Excel**: 通常の `ProteinHunter_results_*.xlsx`。
   [[run_provenance_design]]がPR #21で実装済みのため、`--results`で
   渡されたExcelと同じディレクトリに`<stemと同名>.run_provenance.yaml`
   が存在すれば読み込み、`config_hash`/`git_commit`/`git_dirty`/
   `reference_genome_check_passed`を`calibration_summary.md`に転記する。
   サイドカーが無い場合はこの節を「provenance情報なし」として省略し、
   処理は続行する(必須依存にしない)。

### 出力(PR #23実装後の最終形)

- `<curated>_pairs.xlsx`(シート`Pairs`)。
- `<negatives>_matched.xlsx`(シート`Negatives_Matched`)。
- `calibration_summary.md`(Markdown、変更なし): 各スコア列
  (`interaction_score`, `interaction_priority_score`,
  `string_neighborhood`, `string_cooccurrence`,
  `coexpression_gse77738`, `coexpression_gse64349`, `final_score`)ごとに、

  - Tier A/B別・陰性別の平均・中央値
  - Mann-Whitney U検定のp値
  - ROC-AUC(`U/(n1*n2)`)
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

`--out` ディレクトリに `pairs.xlsx`/`negatives_matched.xlsx`/
`summary.md` を書き出す。日付付きディレクトリにすることで、過去の
校正結果を上書きせず時系列で残せる。

### 既存資産の扱い

- `claude/experimental_interactions_curated.csv` は
  `tier` 列を追加した上でそのまま入力として使う。**Excel化しない。**
- `claude/experimental_interactions_calibration_report.md` は
  「このツール登場以前の、最初の一回限りの分析の記録」として
  そのままリポジトリに残す。冒頭の訂正ブロック(PR #17で追加、
  [[negative_reference_species_mislabeling_finding]]参照)もそのまま
  残す。

### 校正データの拡充そのものについて(ツールとは別に必要な作業)

- Tier B・要確認19ペア: 2026年ACDS論文の特定、CdhC/CdhDパラログ
  クラスタの解消が進み次第、`curated.csv` の該当行の `tier` を
  `B`→`A` に昇格するだけで次回のレポート生成に反映される。
- AlphaFold3陰性セット(現状28件、MA_4115クエリ限定)を、他のTier A
  クエリ(HdrD1, MtpA, NifD/K)についても同様にAlphaFold3で陰性確認
  した候補を追加できると、「陰性はMA_4115限定・正例はMA_4115以外」
  という非対称性を緩和できる。

## スコープ外(あえてやらないこと)

- キャップ/重みの自動最適化は今回はやらない。
- ホールドアウト/交差検証の枠組みも、正例データ数が2桁に達するまでは
  時期尚早と判断し、今回のスコープには含めない。
- 入力CSV(curated.csv・陰性リストCSV)のExcel化(ユーザーの明示的な
  指示により対象外と確定)。

## テスト計画

- `tests/test_calibration_report.py`: 小さな合成データを使い、出力
  される値を手計算した期待値と突き合わせる。Mann-Whitney U統計量と
  AUC換算式の一致、サイドカーYAMLがある場合/ない場合の両方でエラーなく
  動作することの確認、Excel出力の型・書式・シート構成の確認を含む。
- 実データ検証: 既存の`_pairs.csv`と同じ入力データセットを流し、
  `interaction_score`の平均値(39.55/12.95)など既存レポート記載の
  数値と一致することを確認する(回帰確認)。

## 見積もり

新規スクリプト1本(`tools/calibration_report.py`)+ 新規テスト。既存の
スコアリング・出力コードには一切手を入れない独立ツールなので、既存
テストへの影響はゼロ。

## 実装結果(2026-09-20、Claude Code実装・PR #22)

実装・テスト完了(新規22件)、既存分含め全604件パス(PR #21を含まない
mainベース)。PR #21のマージは待たずに実装。CHANGELOGはPR #21と同じ
箇所を編集して衝突するため、この PRでは変更なし(マージ時に統合)。

設計からの差分:

- **scipy不使用。** U統計量・AUC・p値をnumpy/pandasのみで自前実装。
  scipy 1.18と1,500件のランダム標本で全件一致を確認済み(検証用の
  一時インストールのみ、依存には未追加)。
- **シート読み込み対象の修正。** `Interaction_Evidence_Detail`シートは
  現在のExcelに存在しないため、実際の構成(`04_Score_Breakdown`、
  `11_Raw_Audit`、`03_Candidate_Overview`)から読む実装に修正。
- **陰性CSVに`query_old_locus_tag`列(任意)を追加**、既定は
  `--negatives-query`(既定値MA_4115)で補う。
- **`--out`の上書き防止(`--force`なしでは既存出力を上書きしない)**。

データ側の変更(`experimental_interactions_curated.csv`への`tier`列):

- レポートの`tier_final`から機械的に導出(A: 7、B: 20、excluded: 2)。
- **RNAP subunit D・MmcAはexcluded**(問い合わせ可能な遺伝子座番号が
  無いため)。**CdhDクラスターの6ペアはB**
  (`implementation_status_scoring_v2.md`記載の既存の未解決事項を
  そのまま機械可読化したもの)。**両方とも確認・承認済み。**
- AlphaFold3陰性28件を`claude/calibration/af3_negatives_MA_4115.csv`
  として入力用CSV化。

実データ検証:

- 既存レポートとの回帰確認: レポートの`_pairs.csv`/`_negatives.csv`から
  Excelを組み直して流した。Tier Aの`interaction_score`平均39.55・
  中央値43.85(n=8)、陰性12.95・12.27(n=28)、`interaction_priority_score`
  17.97対20.55、`string_neighborhood`0.62対0.04など、レポート記載値と
  一致(レポート記載値の再現であり、パイプラインの独立検証ではないと
  明記)。
- 実際のExcel(5クエリ実行結果)での試走: Tier Aは8ペア一致、
  AF3陰性は28件中4件一致。サイドカー付き実行ではprovenance節が
  サマリに出ることも確認済み。

## 実装結果(2026-09-20、Claude Code実装・PR #23): 出力形式のExcel化

ユーザーからの追加指示(CSV出力を全てExcelに切り替え。対象は本ツールの
出力ファイルのみ、入力CSVは対象外)を受けて実施。全634テストパス。

- `pairs.csv` → `pairs.xlsx`(シート`Pairs`)、`negatives_matched.csv` →
  `negatives_matched.xlsx`(シート`Negatives_Matched`)。入力CSV2種と
  `calibration_summary.md`は変更なし。
- セル型: スコアは数値セル、欠損は空セル(従来は欠損混じりの列が文字列
  扱いになっていた点の副次的な改善)。`found`は真偽値セル。ヘッダ太字・
  ヘッダ行固定・オートフィルタ・列幅上限、といった最低限の書式を適用。
  新規依存なし(`openpyxl`は既存依存)。
- 上書き防止チェックの判定対象を`pairs.xlsx`・`negatives_matched.xlsx`・
  `calibration_summary.md`に変更。古い実行で残った`pairs.csv`等は
  無視され、実行の妨げにならない。
- **2ファイルのまま(1ワークブックに統合しない)と判断。** 理由:
  (1) PairsとNegatives_Matchedは列構成が異なり別々に読まれる、
  (2) ファイル名が設計書のまま維持できる、(3) 上書きチェックが単純に
  保てる。1ワークブック化する案(プロジェクト本体の`ProteinHunter_
  results_*.xlsx`のような単一ワークブック・複数シートの流儀に近づく)
  は将来の選択肢として残す — `write_excel_table`がシート名を引数に
  取る実装になっているため、後から統合する変更コストは低い。
- テスト: CSV読み込みの確認を`pandas.read_excel`に置き換え、セル型・
  書式・1ファイル1シート構成・空表での非破壊・古いCSV残存時の非干渉、
  を新規に確認。
- 回帰確認: 既存レポート数値との回帰テスト(39.55/43.85対12.95/12.27等)
  は変更なく通過。
- 実データ確認: 5クエリの実データExcelで試走し、`calibration_summary.md`
  が変更前の実行結果と(生成時刻の行を除き)完全一致することを確認。
