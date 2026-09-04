# CHANGELOG

ProteinHunter_v5 の変更履歴です。

## 未リリース: 他生物種一致の考慮/非考慮を1スイッチで切り替え(設定.xlsx)

`設定.xlsx`(ユーザー提示の2プリセット比較表)がまとめていた
「他生物種のタンパク質と一致する候補を考慮するか否か」という6項目
(`ranking_metric`、`candidate_sources`の`positive_all_sources`/
`negative_unmatched`/`negative_hit`、`annotation_targets`の`gff`、
`max_candidates_per_query`)を、`config.yaml`冒頭の
`consider_cross_species_matches`という1つのON/OFFスイッチにまとめた。

### Added

- `config.py`: `Config.consider_cross_species_matches: bool = True`を
  追加。`_load_annotation_targets()`/`_load_interaction_scoring()`が
  この値に応じてデフォルト値(`ranking_metric`、
  `candidate_sources.positive_all_sources`/`negative_unmatched`/
  `negative_hit`、`annotation_targets.*.gff`(candidates/
  candidates_relaxedを除く4シート)、`max_candidates_per_query`)を
  切り替える。**個別キーをconfig.yamlで明示した場合は、そちらの値が
  スイッチより優先される**(既存の設定システム全体の慣習と同じ)。
  `candidates`/`candidates_relaxed`/`no_hit`と
  `negative_strong_hit`/`negative_medium_hit`/`negative_weak_hit`は
  このスイッチの対象外(常に既存のデフォルト値のまま)。
  `_validate_consider_cross_species_matches_section()`で型検証を追加。
- `config.yaml`: `project:`直後に`consider_cross_species_matches: true`
  (デフォルト・考慮する)を追加。影響を受ける6項目のうち4項目
  (`annotation_targets`の4シートの`gff`、`candidate_sources`の
  `positive_all_sources`/`negative_unmatched`/`negative_hit`、
  `ranking_metric`)は、スイッチが効くようコメントアウトした参考値に変更
  (`max_candidates_per_query`も同様)。
- `tests/test_config_validation.py`: スイッチがtrue/falseそれぞれの
  プリセットを正しく適用すること、個別キーの明示指定がスイッチより
  優先されること、不正な型を拒否することを検証するテストを追加。

### Changed

- `consider_cross_species_matches`のデフォルト値(true)により、
  `annotation_targets`の`positive_all_sources`/`no_hit`/
  `negative_unmatched`/`negative_hit`の`gff`デフォルトが、
  config.yamlで明示しない場合は`true`から`false`に変わった
  (`candidates`/`candidates_relaxed`は`true`のまま)。これは
  設定.xlsxが定義した「デフォルト(考慮する)」プリセットに合わせる
  意図的な変更。

## 未リリース: Phase 6-8 Stage 2: Wordレポート生成

design spec §29〜45(ユーザー提示の要約に基づく、原文はリポジトリに
存在しない——`claude/phase678_excel_word_redesign_investigation.md`
item 1参照)に基づき、Excelワークブックと並行して単一の`.docx`
Wordレポートを生成する機能を追加。調査・設計・実装記録は
`claude/phase678_stage2_word_report_investigation.md`参照。

### Added

- `requirements.txt`: `python-docx>=1.2`を追加。
- `output/word_narrative.py`(新規): 「Why this candidate ranks highly」
  「Biological Interpretation」を`final_score_tier`・`candidate_source`・
  `negative_hit_strength`・カテゴリ別参照スコア列から条件分岐+穴埋めで
  機械的に生成する純粋関数モジュール。LLM・ネットワーク呼び出しは一切
  行わない(design spec §45再現性: 同一入力なら常に同一文字列)。
  「Biological Interpretation」は必ず「現在利用可能なEvidenceに基づき
  順位付けされた」という枕詞で始まり、「このタンパク質が目的酵素/相互作用
  相手であると確定したものではない」と明記(design spec §35)。両
  `scoring_model`に対応(v2は7カテゴリの参照列、legacy_additiveは
  `scoring_weights`由来の別カテゴリ集合)。
- `output/report_v2.py`: `bookmark_name()`(Word用ブックマーク名を
  `(query_id, candidate_protein_id)`から決定論的に生成、40文字制限は
  ハッシュ付き切り詰めで衝突回避)と`select_top_candidates_per_query()`
  (クエリごとTop-N + Tier1_VeryStrong/Tier2_Strong該当候補は順位に関わらず
  必須採用のセーフティネット)を追加。Excel(`output/excel.py`)とWord
  (`output/word_report.py`)の両方がこの2関数を共有し、互いにデータを
  受け渡すことなく独立して同じ選定・命名結果に到達する。
- `output/word_report.py`(新規): Wordレポート本体の組み立て。
  `python-docx`に高レベルAPIがないブックマーク・外部ハイパーリンク・
  目次(TOC)フィールドを`docx.oxml`直接操作で実装。TOCはWordの
  `w:updateFields`設定(settings.xml)を`true`にすることで、手動の
  「フィールド更新」操作なしに開封時自動更新されるようにした。
  「5. Evidence Architecture」(5.1〜5.7、7エビデンスカテゴリの固定説明文。
  5.5 Evolutionary/5.6 Cellular CompatibilityはPIH bundle未設定時に
  正直に「データ未提供」と表示し、5.7 NegativeはPIH bundleの有無に
  関わらず常に「本バージョンには未実装」と明記——前者はconfigで解消
  可能なデータギャップ、後者は現時点で解消不能な機能ギャップとして
  意図的に書き分けた)、「7. Candidate Ranking」(クエリごとの要約表、
  design spec §29「1パイプライン実行につき1ファイル」を「複数クエリは
  サブセクション分割」で満たす)、「8. Candidate Details」(クエリ×候補
  ごとにブックマーク付き見出し+2つの生成文+Excel参照行)を出力。
- `output/excel.py`: `write_classification_workbook()`に
  `word_report_filename`引数を追加(既定`None`、後方互換)。指定時、
  `select_top_candidates_per_query()`で選ばれた候補の`word_report_link`
  列に`f"{filename}#{bookmark_name(...)}"`形式の実クリック可能な
  Excelハイパーリンクを設定(Stage 1で予約済みだった列、design spec
  Excel⇔Wordリンクのoption A+C)。
- `main.py`: Excel出力の直後にWordレポート出力セクションを追加。
  `paths.output_word`未設定時は`output_excel`と同じディレクトリ・
  同じstemに`.docx`拡張子を付けたパスを既定値として使用(既存の
  config.yamlに新キー追加を強制しない)。
- `config.py` / `config.yaml`: `PathConfig.output_word`(任意)と
  `InteractionScoringConfig.word_report`(`WordReportConfig`:
  `enabled`・`max_candidates_per_query`、既定15)を追加。
  `max_candidates_per_query`は表示件数のみを制御し、スコアリング・
  Excelワークブックには一切影響しない。
- `core/exceptions.py`: `WordReportError`を追加(`ExcelOutputError`と
  同じパターン)。

### Design decisions confirmed by the user (this session)

- クエリの並び順はconfig記載順(query_id昇順ではない)。
- Top-N既定値は**15**、実データでの精密なTier分布再検証は今回は不要
  (表示件数調整用パラメータでありスコアリングに影響しないため、暫定値で
  運用し後日調整)。
- M6の実データ検証はPIH bundleなしで実施(PIH連携自体が別途未着手の
  検証事項であるため)。

### Real-data verification (M6)

Stage 1と同一の5クエリ(MA_0688・MA_4165・MA_3898・MA_3899・MA_4115)・
`v2_evidence_based`・STRING PPI(taxid 188937)+GEO coexpression有効・
PIH bundleなし・`candidate_sources`に`negative_hit`含む・
`ranking_metric: final_score`という設定で、実BLAST/CDD/STRING/GEOデータに
対して実行(`.cache/word_report_verification/`、gitignore対象のため
非コミット、`config.yaml`自体は無変更)。総実行時間974.5秒(うちCDD注釈が
913.5秒、Word出力自体は0.45秒)。

- 4,627件のtarget中151件が`Candidates`シート該当、`interaction_scoring`は
  5クエリ全件で成功し、`02_Final_Score`は2,775行(151候補×クエリ、一部
  除外あり)。Excel・Word両方とも正常終了、クラッシュなし。
- **`final_score_tier`の実分布(2,775行)**: `Tier4_Weak` 1,566件、
  `Tier3_Moderate` 1,201件、`Tier2_Strong` **8件のみ**、
  `Tier1_VeryStrong` **0件**。design spec §29〜35調査時点(実データ分布
  未取得)で示した懸念——「Tier1〜2のみに絞ると空になりかねない」——が
  実データで正確に裏付けられた。今回採用したTop-15(順位ベース)+
  Tier1/2セーフティネットの設計判断は妥当だったことを確認。
- Word側候補数: 5クエリ×最大15 = 75候補(実際も75件、セーフティネットが
  追加した候補はゼロ——今回のTier2該当8件は全てもとから各クエリの
  Top15順位内に収まっていたため)。`.docx`ファイルサイズ49KB(全2,775行を
  ミラーした場合との対比で、「Top-N要約+Excel誘導」設計の妥当性を確認)。
- Excel `word_report_link`列: 選定された75候補全件で
  `f"{filename}#{bookmark}"`形式の実クリック可能ハイパーリンク
  (`cell.style == "Hyperlink"`)を確認、対応するWord側ブックマーク
  (`w:bookmarkStart`)が全件存在することも確認。
- Word設定XML(`settings.xml`)に`w:updateFields=true`が実ファイルにも
  正しく書き込まれていることを確認(手動フィールド更新なしでのTOC自動更新)。
- 生成文のサンプル確認(`candidate_source=No_hit`・`final_score_tier=Tier2_Strong`
  の候補1件): 「Why ranks highly」がカテゴリ別の実スコア・実配点上限を
  正しく列挙、「Biological Interpretation」がdesign spec §35の
  必須枕詞("not a confirmed identification..."）およびEvolutionary/
  Cellular Compatibility/Negativeの正直な限界表示を含むことを目視確認。
- `07_Evolutionary_Evidence`・`10_Negative_Evidence`および対応する
  Word「5.5」「5.6」「5.7」節は、想定通りPIH bundle未設定のため
  「データ未提供」(5.5/5.6)・「未実装」(5.7)の文言で表示。

## 未リリース: Phase 6-8 Stage 1: Excelワークブック12シート再設計

design spec §25.1に準拠し、Excel出力全体を12シートへ統合。従来の
ベース分類シート(Candidates/Candidates_relaxed/...の最大10枚)と
`Interaction_*`バケツシート(最大11枚)を、`candidate_source`列を持つ
統合済み行として再構成した。調査・設計・実装記録は
`claude/phase678_excel_word_redesign_investigation.md`参照。
Wordレポート生成(design spec Phase 6-8のStage 2)は別タスクとして後日
着手予定——本PRには含まれない。

### Added

- `output/report_v2.py`(新規): `candidate_source`統合レイヤー。
  base分類バケツ(Candidates, Candidates_relaxed, ..., Negative_hit)と
  `Interaction_*`バケツの両方について、同じ候補が複数の有効なバケツで
  重複して現れる場合(例: `Candidates`は`Candidates_relaxed`の部分集合、
  PR #11で判明した`negative_hit`と`negative_strong/medium/weak_hit`の
  重複)、固定の優先順位(`CANDIDATE_PRIORITY_BASE`を再利用: Candidates
  > Positive_all_sources > Candidates_relaxed > No_hit >
  Negative_unmatched > Negative_hit)で1行に重複排除する。
  `negative_strong/medium/weak_hit`は表示上`Negative_hit`へ統合され、
  強度は`negative_hit_strength`列(下記)に残る。
- `analysis/interaction_scoring.py`: 全ペア行に`negative_hit_strength`列
  を追加(両モデル、旧`Negative_strong/medium/weak_hit`専用シートの
  代替)。`scoring_model: v2_evidence_based`限定で
  `functional_domain_score`・`evolutionary_score`・
  `cellular_compatibility_score`・`interaction_evidence_score`
  (`external_ppi_evidence`+`coexpression_evidence`+
  `pih_direct_interaction`の合計)をカテゴリ単位の参照列として追加
  (既存の`candidate_priority_score`/`same_gene_neighborhood_score`と
  同じパターン)。クエリなし/interaction_scoring無効時のFinal Score
  フォールバック用に`compute_protein_hunter_only_final_score()`を新設
  (公開関数)。いずれも既存のスコア・順位には影響しない純粋な追加。
- `output/excel.py`: `write_classification_workbook()`を書き換え、
  固定12シート(`01_Index, 02_Final_Score, 03_Candidate_Overview,
  04_Score_Breakdown, 05_Sequence_Evidence,
  06_Functional_Domain_Evidence, 07_Evolutionary_Evidence,
  08_Genomic_Context, 09_Interaction_Evidence, 10_Negative_Evidence,
  11_Raw_Audit, 12_Reserved`)を出力するよう変更。シグネチャを
  `(config, blast_classification, output_path, interaction_result)`
  へ簡素化(旧: 個別バケツ辞書を12個の引数として受け取る形式)。
  `05_Sequence_Evidence`〜`10_Negative_Evidence`は既存の
  `Interaction_Evidence_Detail`(v2ロングフォーマット)を
  カテゴリでフィルタしたビュー(`scoring_model: v2_evidence_based`限定、
  `legacy_additive`ではヘッダーのみ)。`11_Raw_Audit`は旧
  `Interaction_Evidence_Detail`全件に加え、旧`Interaction_query`・
  `Interaction_Neighborhood`の内容を同一シート内に積み重ねたブロックとして
  保持し、情報を失わない。`12_Reserved`はStage 2用の意図的な空シート。
  `Positive_source_summary`シートは廃止(`03_Candidate_Overview`が
  同等以上の内容を含むため)。
- 既存バグ修正: `_add_back_to_index_link`が`"Index"`シート名を
  ハードコードしており、`01_Index`への改名で全シートの「Back to
  Index」リンクが壊れる状態だった。修正済み。

### Real-data verification

Tier A(8ペア)+ AlphaFold3陰性(28件、MA_4115クエリ)を用いた実データ
検証(`.cache/geo_investigation/`、gitignore対象のため非コミット)。

初回の検証はSTRING PPI/GEO coexpressionを無効化した設定で行われ、
Final Score統合フェーズの検証(+17.20)と条件が一致していなかった
(陰性側が外部エビデンスを失い分離幅が見かけ上+32.89まで拡大しただけ
で、正例側やスコアリングロジックの変化ではない)。`git worktree`で
Stage 1適用前のコミット(`466564a`, PR #12マージ後)をチェックアウトし、
Final Score統合フェーズの検証と同一設定(STRING PPI・GEO coexpression
両方有効)で再実行、Stage 1適用後のブランチとも同一設定で実行して
36ペア全件(Tier A 8 + AF3陰性28)を`final_score`/`interaction_score`
まで直接突合した結果、**全件が小数点以下まで完全一致**
(POS平均42.884、NEG平均25.680 → 分離幅 **+17.20**、Final Score統合
フェーズの数値を厳密に再現)。シート統合によるスコア値への影響は
ゼロであることを確認。副産物として、MA_3899→MA_3898ペアが
`Interaction_Candidates_relaxed`と`Interaction_Neg_hit`の両方に
同一スコア(54.033)で重複出現していたこと(PR #11の重複スコアリング
問題そのもの)も実データで直接確認し、統合ロジックが正しく
`Candidates_relaxed`側を採用することを確認済み。重複排除・
"Unclassified"候補ゼロも確認済み。

## 未リリース: Final Score統合(design spec §17-22・§27)

`protein_hunter_score`と`interaction_score`を1つの数値へ統合する
"Final Score"を実装。調査・設計・承認・実データ検証・設計修正の全記録は
`claude/final_score_integration_investigation.md`参照。

### Added

- `analysis/interaction_scoring.py`: `analysis/scoring_engine.py`の
  カテゴリcap方式を再利用し、`protein_hunter_score`(正規化上限18、理論値
  固定)と`interaction_score`(0-100)をそれぞれ独立したトップレベル
  カテゴリ(`protein_hunter`, cap 30 / `interaction`, cap 70、いずれも
  暫定値)として扱う`final_score`を新設。`interaction_score`が存在しない
  場合(クエリ固有証拠なし)は`protein_hunter`カテゴリ単独で自動的に
  再正規化(`score_candidate`既存の「利用可能な証拠だけで再正規化」の
  仕組みをそのまま利用、特別分岐なし)。`final_score`・`final_score_tier`
  (既存`evidence_tier`と同じ閾値ロジックを再利用した別列)を
  `Interaction_*`シートに追加。両スコアリングモデル
  (`legacy_additive`/`v2_evidence_based`)で計算。既存の
  `interaction_priority_score`・`Evidence_Tier`・デフォルトの
  `ranking_metric`は無変更。
- **`negative_hit_strength`ペナルティはFinal Scoreには適用しない
  (実装→実データ検証→方針転換の結果、意図的な設計判断)**。
  `negative_hit_strength`(`analysis/ortholog_filter.py`)は系統特異性/新規性
  のシグナル(このタンパク質が陰性参照ゲノムにも広く存在するありふれた
  ものか)であり、design spec §7.7が定義する「Negative Evidence」
  (functional contradiction、incompatible localization、incompatible
  domain、phylogenetic contradiction——「このペア自体が相互作用として
  矛盾している」という反証)とは概念として別物。当初この2つを混同して
  ペナルティとして組み込んだ結果、古くから保存された中心代謝系
  (Hdr/Mtp/Nif複合体)は`negative_hit_strength = strong`になりやすい
  (＝`Negative_hit`バケツに分類される基準そのもの)ため、真の相互作用
  パートナーを一律に減点する逆効果が実データ検証で判明(PR #11で発見した
  「negative_hitバケツの盲点」と同根の問題がスコア計算の場面で再発)。
  `final_score_negative_penalty`コンポーネント自体(監査列の枠組み)は
  残すが、値は常に`NOT_APPLICABLE`——design spec §7.7が本来意図する真の
  生物学的矛盾シグナルが将来実装された際に使う予約枠であり、
  `negative_hit_strength`はその代用にはならないと判断。
  `interaction_priority_score`側の既存の`negative_hit_strength`適用は
  変更なし(候補全体の妥当性評価として妥当な設計のため)。
- `interaction_scoring.ranking_metric`に`final_score`を追加(既定は
  `interaction_priority_score`のまま変更なし)。既知の制約:
  `protein_hunter_score`は各クエリの候補が既に(別の指標で)ランキング・
  `max_candidates_per_query`件へ切り詰められた後にしか分からないため、
  `final_score`による再ランキングは切り詰め後の候補集合の中でのみ行われる
  (`interaction_priority_score`/`interaction_score`は切り詰め前から利用可能
  なため、この制約を受けない)。
- 後方互換性: Final Score統合以前に作成されたカスタム
  `scoring_engine_config.yaml`(`protein_hunter`/`interaction`の
  cap未定義)を使っていても、`ConfigError`にならず、モジュール既定の
  暫定cap(30/70)にフォールバックする(`_final_score_engine_config`)。
- `Interaction_Evidence_Detail`シート(v2のみ)に、Final Scoreの3成分
  (`protein_hunter_score`, `interaction_score`,
  `final_score_negative_penalty`——常に`NOT_APPLICABLE`)の
  raw/normalized/contribution内訳を追加(design spec §22・§24が求める
  追跡可能性)。
- テスト7件追加(`tests/test_interaction_scoring.py`: Final Scoreの
  基本計算・フォールバック・`negative_hit_strength`が影響しないことの確認
  (`NOT_APPLICABLE`)・後方互換性・legacy_additive対応・
  `ranking_metric: final_score`の再ランキング、`tests/test_config_validation.py`:
  `final_score`が有効な`ranking_metric`値として受理されることの確認)。

### 実データ検証

Tier A正例8ペア・AlphaFold3陰性28件で`interaction_score`単独と
`final_score`を比較しました。

| 指標 | 正例(n=8) mean/median | 陰性(n=28) mean/median | 分離幅 |
|---|---:|---:|---:|
| `interaction_score`単独 | 39.83 / 41.22 | 12.28 / 12.59 | +27.55 |
| `final_score`(ペナルティ適用時、修正前) | 16.05 / 13.86 | 16.00 / 17.94 | +0.05(ほぼ消失) |
| `final_score`(ペナルティ除外後、最終) | 42.88 / 43.86 | 25.68 / 25.23 | **+17.20** |

当初`negative_hit_strength`ペナルティを独立適用した実装では分離幅が
ほぼ消失しましたが(上記「Added」参照、原因分析の通り)、ペナルティ除外後は
`interaction_score`単独の分離幅(+27.55)に対し`protein_hunter_score`混合分
だけやや下がる程度(+17.20)に回復し、想定通りの挙動になったことを確認しました。
cap配分(30/70)自体は変更していません。詳細は
`claude/final_score_integration_investigation.md`の「実装後の実データ検証」
「negative_hit_strengthペナルティの除外」「修正後の再検証」節、生データは
`claude/final_score_verification_positive.csv` /
`final_score_verification_negative.csv`参照。

## 未リリース: 実験的相互作用データによる実データ検証と2件の修正

`Methanosarcina_acetivorans_experimental_protein_interactions.xlsx`(27行、
公開文献ベースの実験的相互作用証拠)をキュレーション・`old_locus_tag`対応付け
した上で(`claude/experimental_interactions_curation.md`)、既存の
AlphaFold3陰性校正データ(28件)と同一条件のパイプライン実行で直接比較する
診断run(`claude/experimental_interactions_calibration_report.md`)を実施。
2つの重要な発見があり、そのうち1件はコード修正、もう1件はドキュメント対応
とした。

### Fixed

- v2の`interaction_score`合算対象(`INTERACTION_SCORE_COMPONENT_NAMES`)から
  `coexpression_gse77738`を除外。実データ検証(陽性8件・陰性28件)で、
  この指標が真の相互作用ペアよりも非相互作用ペアの方が高い値を示す逆転
  (陽性側 平均0.480/中央値0.388 vs 陰性側 平均0.655/中央値0.732。最も低い
  陽性値1件を除いても陽性0.548/0.416 vs 陰性のまま逆転は残る)が確認された
  ため。データ取得・キャッシュ(`analysis/coexpression_bridge.py`)・
  `Interaction_Evidence_Detail`シートへの表示は変更なし——スコアには使わない
  が参考情報としては引き続き見える。`coexpression_gse64349`はこの逆転が
  見られなかった(陽性0.848 vs 陰性0.601)ため現状維持(weight 1/3のまま)。
  実データでの修正確認: Tier A 8ペア全件で`coexpression_gse77738`が
  Evidence_Detailに残ることと、`interaction_score`が期待通り(除外前より
  高い値だったペアは低下、低い値だったペアは上昇)変動することを確認。

### Documented

- 診断runで、HdrD1・Mcr複合体・Nifシステムなど保存性の高い経路/複合体の
  既知相互作用パートナーが軒並み`Candidates`バケツではなく`Negative_hit`
  (またはそのサブバケツ)に分類されることが判明。これらは古くから保存された
  中心代謝酵素のため、陰性参照ゲノム側にも強くBLASTヒットしてしまうことが
  原因。`candidate_sources.negative_hit`は既定で無効なため、**デフォルト
  設定のままでは、このような保存性の高いクエリタンパク質の真の相互作用
  パートナーが`Interaction_*`出力に一切現れない**。既定値は変更せず、
  `config.yaml`の`negative_hit`周辺にこの挙動と対処法(有効化の検討)を
  コメントで追記。`negative_hit`とその3つのサブバケツ(`negative_strong/
  medium/weak_hit`)を同時に有効化すると同一候補が重複してスコア・出力
  される(実例: MtpA-MtpCペアがInteraction_Neg_hitとInteraction_Neg_strong
  の両方に別々の`candidate_rank`で出現)ため、その旨の注意書きも追加。
  `config.redundant_negative_hit_sources()`(`main.py`から呼び出し)で
  この組み合わせを検出した場合はrun時に警告ログを出力する(configエラー
  ではなく警告のみ——意図的な併用を想定したユースケースもあり得るため)。
- 副次的な発見として、`MA_1111`(RNAP subunit D)は現行RefSeq注釈で
  pseudogene(フレームシフト)のため`target.faa`に配列自体が存在せず、
  クエリにも候補にも一切使えないことを確認(前回のCHANGELOGエントリ
  「公開共発現データの統合」に記載のGFF調査を裏付ける実行時の確認)。

## 未リリース: 公開共発現データの統合(Phase 6b)

GEO公開RNA-seqデータ(GSE77738/GSE64349)由来の実測共発現証拠を
`interaction_score`に追加する対応。実装前に両データセットの補助ファイルを
実際にダウンロード・検証し、その結果を`claude/phase6b_coexpression_design.md`
に記録した上で段階実装(M1〜M3、M4は見送り、M5)した。詳細は同ドキュメント
参照。

### Added

- `analysis/coexpression_bridge.py`: GSE77738/GSE64349の加工済み補助ファイル
  (XLS/XLSX)をダウンロード・ローカルキャッシュし、クエリ遺伝子ごとに
  他の全既知遺伝子とのPearson相関 + そのクエリ自身の背景相関分布内での
  パーセンタイル順位を計算、`core/cache.py::JsonCache`
  (名前空間`coexpression_gse77738`/`coexpression_gse64349`)にクエリ単位で
  キャッシュするブリッジを新設。
  - GSE77738は61サンプルのうちアクチノマイシンD処理によるRNA分解タイム
    コースの時系列点を除外し、真に独立な定常状態サンプル13個のみを使用。
  - GSE64349はΔmsrH変異株サブセット(TableS2の一部)を除外し、TableS2の
    「WWM82(親株)」サブセットは追加の野生株レプリカとして採用(TableS1の
    9サンプル+3サンプル=12サンプル)。
  - 遺伝子ID(`MA0001`形式)は`old_locus_tag`(`MA_0001`形式)への単純な
    アンダースコア挿入で変換可能。GSE64349の遺伝子シンボル表記(`cdc6_1`等)
    はGSE77738自身のGene Name列から構築したルックアップテーブルで解決。
- `analysis/interaction_scoring.py`: v2に新カテゴリ`coexpression_evidence`
  (暫定cap 12点)、`coexpression_gse77738`と`coexpression_gse64349`の
  2コンポーネント(STRINGの`string_cooccurrence`/`string_neighborhood`と
  同様、カテゴリcapを共有する別コンポーネントとして追加、プールしない)。
  GSE64349はサンプル数が少なく(12サンプル・4条件)、パーセンタイル正規化
  後もなお統計的信頼性が低いため、weightをGSE77738の1/3に暫定的に低減
  (`V2_COMPONENT_WEIGHTS["coexpression_gse64349"]`)。いずれも
  `interaction_score`に算入(`INTERACTION_SCORE_COMPONENT_NAMES`に追加)。
  legacy_additiveへの適用は今回見送り(Phase 6aのM4と同様、別対応)。
- `config.py`: `interaction_scoring.geo_coexpression_enabled`(既定false=
  無効)。STRINGのtaxidと異なり、データセット自体は固定のためon/offのみ。
- 正規化方式: 固定的な線形マッピングではなく、クエリ遺伝子ごとの背景相関
  分布に対するパーセンタイル順位を採用。実データ調査でGSE64349単独では
  ランダム遺伝子ペアの背景相関平均が0.76という深刻なインフレを起こす
  (GSE77738は0.15〜0.21)ことが判明したため。
- MISSING(遺伝子がデータセット自体に存在しない、または分散ゼロで相関が
  定義できない)とAVAILABLE(存在すれば相関が弱くても評価済みとして扱う)
  の区別は、Phase 6aのSTRING証拠と同じ既存パターンをそのまま踏襲(新しい
  ステータスは追加していない)。
- `output/excel.py`: Indexシートに`coexpression_gse77738`/
  `coexpression_gse64349`の列説明と、GEOの利用条件(NIH公的データベース、
  ライセンスによる強制はないが元論文PMID 27852217/25691524のクレジットを
  慣例として記載)を追加。
- 実データ検証(M_acetivorans、MA_4115、v2_evidence_based、実際の
  Candidatesバケツ148件、taxid 188937のSTRING証拠と実際に相互比較):
  - GSE77738/GSE64349とも実際にダウンロード・パースし、148件中それぞれ
    142件/144件がAVAILABLE(coverage 96%/97%)であることを確認。
  - MA_4115とその遺伝子近傍(MA_4114/4116/4117、STRINGでも上位パートナー)
    の相関を、STRING証拠(Phase 6aで既に取得済み)とは完全に独立な実測
    発現データで計算し直したところ、3件ともGSE77738パーセンタイル
    0.93〜0.98、GSE64349パーセンタイル0.16(MA_4114のみ弱い)〜0.99と、
    STRINGのneighborhood/cooccurrence証拠(MA_4114: 606/531、MA_4116:
    849/0、MA_4117: 696/0)を独立なデータソースから裏付ける結果を得た。
  - 実際の148件Candidatesバケツの中でGSE77738パーセンタイル上位10件は、
    いずれもSTRINGのneighborhood/cooccurrenceが共にゼロ(STRINGが何の
    関連も検出していない候補)であり、共発現証拠がSTRINGと重複しない
    独立シグナルを提供していることを確認。
  - 一方、AlphaFold3校正28件中「確信度の高い陰性」25件(AF3で直接相互作用
    なしと確認済み)のGSE77738パーセンタイルは中央値0.767・平均0.654
    (GSE64349も同水準)と、判定なし(0.5)を上回る水準にとどまり、一部
    (MA_1447: 0.98、MA_4116: 0.98など)はむしろ非常に高いパーセンタイル
    を示した。共発現(転写共制御)は物理的相互作用を保証しないという
    生物学的に予想通りの限界であり、この証拠区分のcap/weightを暫定的に
    抑えめに設定した判断(coexpression_evidence cap 12点、GSE64349は
    さらに1/3)の妥当性を裏付ける結果として記録する。
- テスト追加(`tests/test_coexpression_bridge.py`新設、
  `tests/test_interaction_scoring.py`, `tests/test_config_validation.py`)。

### Notes

- `claude/phase6b_coexpression_design.md`に記載の通り、GSE66445
  (代謝改変株)は今回も除外。GSE64349のΔmsrH変異株サブセットも同じ理由で
  除外。
- GEOの正式な利用条件ページはプログラムからの確認がNCBI側のreCAPTCHAで
  ブロックされたため未検証。STRINGのようなライセンスによる強制表示義務は
  無いという一般的理解のもと、元論文PMIDのクレジット表記のみ対応。

## 未リリース: STRING PPI証拠の統合(Phase 6a)

外部知識ベース(STRING)由来の証拠を`interaction_score`に追加する対応。
実装前にSTRING API/一括ダウンロードファイル/公開GEOデータを実際に調査し、
その結果を`claude/phase6_external_evidence_design.md`に記録した上で
段階実装(M1〜M5)した。詳細は同ドキュメント参照。

### Added

- `analysis/string_ppi_bridge.py`: STRINGの種別一括ダウンロードファイル
  (`protein.links.detailed`/`protein.info`)を取得・ローカルキャッシュし、
  クエリ単位で`core/cache.py::JsonCache`(名前空間`string_ppi`)に結果を
  保存するブリッジを新設。未キャッシュの生物種向けライブAPIフォールバック
  (`interaction_partners`、`caller_identity`はアプリ名、呼び出し間隔1秒)
  も実装。ネットワーク障害時は例外を投げず空の証拠として扱う(PIH bridge
  と同じ「外部証拠は失敗してもローカル実行を止めない」方針)。
- `analysis/interaction_scoring.py`:
  - v2: 新カテゴリ`external_ppi_evidence`(暫定cap 15点)と
    `string_cooccurrence`コンポーネント(STRINGの`cooccurrence`チャンネル)。
    `string_neighborhood`(`neighborhood`チャンネル)は既存の
    `genomic_context`カテゴリのcapを共有する第2コンポーネントとして追加。
    いずれも`interaction_score`に算入(`INTERACTION_SCORE_COMPONENT_NAMES`
    に追加)。
  - legacy: 単一の`string_ppi_score`(cooccurrence+neighborhoodの平均、
    新設`scoring_weights.external_ppi`で正規化)を`interaction_priority_score`
    と`interaction_score`双方に算入。STRING未設定時は既存runの結果を一切
    変えないよう、分母への算入もSTRING有効時のみに限定。
- `config.py`: `interaction_scoring.string_ppi_ncbi_taxon_id`(既定unset=
  無効)。**種レベルのNCBI taxidではなく、STRING独自の株レベルtaxidが
  必要**(この生物種の場合、種レベル2214ではなく株レベル188937)。
- `analysis/scoring_engine_config.py` + `config/scoring_engine.example.yaml`:
  `DEFAULT_CATEGORY_CAPS`に`external_ppi_evidence: 15.0`(暫定値)を追加。
- MISSING(STRINGがそのタンパク質のデータを一切持たない)と
  evaluated-zero(両者ともSTRING既知だがこのペアの行が一括ファイルに
  存在しない→計算済みでスコア0とみなす)を区別。
- `output/excel.py`: IndexシートにSTRINGの列説明とCC BY 4.0クレジット表記
  を追加。
- 実データ検証(M_acetivorans、MA_4115、v2_evidence_based、taxid 188937):
  実際にSTRING一括ファイルをダウンロード・パースし、AlphaFold3校正20件
  (Interaction系シート掲載分)全てがSTRINGに既知(MISSING無し)である
  ことを確認。うち16件(No_hitバケツの4件はEvidence_Detail既定除外のため
  未確認)でSTRING証拠は概ねゼロだったが、4件(MA_1978/MA_1447/MA_0545/
  MA_0164)では非ゼロのSTRING証拠により`interaction_score`が0から
  3.9〜11.75へ上昇。当初の焦点だったMA_0050/MA_0238はSTRING証拠も
  ゼロのままで`interaction_score=0`を維持(悪化なし、想定通り)。
  STRINGを有効化すると、たとえ証拠がゼロでも「評価済みカテゴリ」として
  分母(total_cap)に算入されるため、`interaction_priority_score`は
  STRING既知の全候補でわずかに変動する(例: MA_0050は38.786→32.322)。
  これは設計通りの挙動であり、STRINGを有効化する際の既知の影響として
  記録しておく。
- テスト34件追加(`tests/test_string_ppi_bridge.py`新設、
  `tests/test_interaction_scoring.py`, `tests/test_config_validation.py`,
  `tests/test_exceptions.py`, `tests/test_excel_output.py`,
  `tests/test_scoring_engine_config.py`)。

### Notes

- 当初計画していた「STRINGのfusion evidenceチャンネルでRosetta Stone法を
  代替する」方針は撤回。実データ調査でこの生物種のfusion証拠が0%だった
  ため(`claude/phase6_external_evidence_design.md`参照)。機能がよく
  分かった別のクエリを扱う際に改めて検討する。
- 公開共発現データ(GEO、M. acetivorans C2A株、最大82サンプル)は
  Phase 6bとして別途対応予定。今回は未実装。

## 未リリース: protein_hunter_score / interaction_score の分離(design spec 22章)

`Interaction_Evidence_Detail`シート(直前の変更)を使ったAlphaFold3校正データ
監査の結果、MA_0050/MA_0238(AF3で相互作用なしと確認済み)が高スコアになる
根本原因は、「良いProteinHunter候補である」ことと「クエリと具体的に相互作用
する」ことが`interaction_priority_score`という1本の合成スコアに未分離のまま
混ざっていたことだった。この2本を分離する対応(M1〜M5、段階的に実装)。

### Added

- `analysis/scoring.py`: `build_candidate_score`(`ProteinRecord`を変更せず
  `CandidateScore`を計算する純粋関数)を追加。`score_record`はその薄いラッパー
  に変更(挙動は完全に維持)。
- `analysis/interaction_scoring.py`:
  - `resolve_protein_hunter_scores`(M1)を追加。既存の「Candidate scoring」
    パイプラインステップが`positive_only_records`("Candidates")にしか
    `protein_hunter_score`を計算していなかったスコープの穴(PR #4がCDD注釈で
    修正したのと同型の問題)を、`interaction_scoring`が触れる全候補に拡張。
    `ProteinRecord.score`は変更しない(Candidates_relaxed/No_hit等の既存
    分類シートに影響させないため)。
  - `protein_hunter_score` / `protein_hunter_score_components` /
    `protein_hunter_score_reasons` を全Interaction_*シートの各行に参照列
    として追加(M2)。
  - `INTERACTION_SCORE_COMPONENT_NAMES = {genomic_context, domain_complementarity}`
    を新設し、`interaction_score` / `interaction_evidence_tier`
    (`scoring_model: v2_evidence_based`、M3)を追加。
    `analysis/scoring_engine.py::score_candidate`を無変更のまま、コンポーネント
    リストをこの2つに絞って再度呼ぶだけでcap再正規化済みスコアを取得している。
    `co_occurrence`は意図的に除外: クエリを引数に取るが実質的には各タンパク質
    自身のpositive参照ゲノムへのBLASTヒットパターンを比較しているだけであり、
    かつこのプロジェクトのデフォルト設定ではpositiveソースが2つしかないため
    Jaccard値は理論上0.0/0.5/1.0の3値しか取れない(監査目的で
    `Interaction_Evidence_Detail`には従来通り残す)。
  - legacy_additive向けの`interaction_score`(M4): `(same_gene_neighborhood_score
    + domain_complementarity_score) / (weights.gene_neighborhood +
    weights.domain_complementarity) * 100`。`interaction_evidence_tier`は
    legacyにカテゴリ数の概念がないため常に空欄。
  - `interaction_scoring.ranking_metric`(M5、`config.py`に追加、既定値
    `interaction_priority_score`)。`interaction_score`に設定すると
    `candidate_rank`・行順序だけがクエリ特異的証拠のみで決まるようになり、
    `interaction_priority_score`/`evidence_tier`/`priority_group`等の値は
    どちらの設定でも変化しない。
- `config.py`: `InteractionEvidenceDetailConfig`と対になる形で
  `InteractionScoringConfig.ranking_metric`を追加(validate/load込み)。
- 実データ検証(M_acetivorans、クエリ=MA_4115、v2_evidence_based)で、
  AlphaFold3校正データ28件中20件(残り8件は元々どのInteraction系シートにも
  未掲載)を確認。MA_0050/MA_0238を含む全20件が
  `interaction_score = 0` / `interaction_evidence_tier = Tier4_Weak` に
  統一される一方、`protein_hunter_score`は9〜14(候補としての一般的な質は
  維持)、`interaction_priority_score`/`evidence_tier`は既存値のまま不変
  であることを確認。設計通りの分離が実現された。
- テスト29件追加(`tests/test_scoring.py`, `tests/test_interaction_scoring.py`,
  `tests/test_config_validation.py`)。

## 未リリース: Interaction_Evidence_Detail シート(スコア内訳の監査用出力)

`interaction_priority_score` がどのエビデンスカテゴリ・コンポーネントに
支えられているかを、既存のInteraction系シートを変更せずに追加のシートで
確認できるようにする対応。

### Added

- `analysis/interaction_scoring.py`: 新規シート `Interaction_Evidence_Detail`
  を追加。`scoring_model: v2_evidence_based` では
  `ScoreBreakdown.components`(`analysis/scoring_engine.py` で既に保持されて
  いた各 `EvidenceComponent` の raw_value/normalized_value/weight/status/
  explanation)を1コンポーネント1行のロング形式で展開。`legacy_additive`
  では既存の5つの内訳スコア(`candidate_priority_score` 等)と
  `interaction_score_reasons` を1ペア1行のワイド形式でそのまま射影する
  (新規計算なし)。対象範囲は既存Interaction系シートに掲載される候補
  (上位N件切り詰め後)と同一。
- `config.py`: `InteractionScoringConfig.evidence_detail_sheet`
  (`InteractionEvidenceDetailConfig.include_no_hit`, デフォルト `false`)を
  追加。`no_hit`バケツは件数が最も多く、v2ではスコアがタイになりやすいため
  既定で詳細シートから除外し、必要な場合のみ `include_no_hit: true` で
  含められるようにした。
- `output/excel.py`: `_interaction_dataframes` に
  `Interaction_Evidence_Detail` の書き込みを追加(詳細行が0件の実行では
  シート自体を作らない)。既存シートの列・内容は無変更。
- `tests/test_interaction_scoring.py`, `tests/test_excel_output.py`,
  `tests/test_config_validation.py`: v2のロング形式(1ペア=6コンポーネント)、
  legacyのワイド形式、`no_hit`既定除外、`include_no_hit: true`時の挙動、
  複数バケツ混在時のスコープ一致、既存シートへの非破壊性、config検証の
  テストを追加。

## 未リリース: sequence_evidence (BLAST hit強度のスコアリング反映)

統合設計書3章・9章の指摘「非常に強いBLAST hitと弱いBLAST hitが、同じpositive hit
として同程度に扱われる可能性がある」への対応。`scoring_model: v2_evidence_based`
限定の追加で、デフォルト(`legacy_additive`)挙動には影響しません。

### Added

- `analysis/interaction_scoring.py`: `_build_evidence_components_v2` に
  `sequence_evidence` コンポーネントを追加。`ProteinRecord.positive_hits` の
  代表ヒット(`analysis/candidates.py::get_best_hit` の既存
  `(bitscore, -evalue)` ルールをそのまま再利用、新しい集約ロジックは追加せず)
  のidentity/coverage/evalueを0.0-1.0の強度値へ正規化し、`source_classification`
  カテゴリの30点キャップを `source_classification` コンポーネントと共有します
  (`co_occurrence`/`domain_complementarity` が `functional_annotation` を
  共有している既存パターンと同一)。`positive_hits` が空の候補は `MISSING`
  として扱い、0点として減点しません。bitscoreは今回のスコアリングには
  使用していません(配列長依存でこのプロジェクトに校正済みの絶対閾値が
  ないため)。
- `analysis/scoring_engine_config.py` + `config/scoring_engine.example.yaml`:
  `SequenceEvidenceConfig`(identity/coverage/evalueのfloor・ceiling・
  サブ重み)を追加。identity/coverage floorは `ortholog_filter.weak` の
  既存閾値(25.0% / 50.0%)を、evalue参照上限は `config.blast.evalue` の
  既定カットオフ(1e-5)をそれぞれ転用しています(詳細と根拠は
  `docs/implementation_plan_sequence_evidence.md` を参照)。
- evalue == 0.0(BLASTの浮動小数点丸めによる完全一致相当の値)は
  `-log10` を計算せず、直接最強スコア(1.0)として扱うガードを追加。
- `tests/test_interaction_scoring.py`: 強いhit/弱いhit/hitなし(MISSING)/
  複数hit(代表選択)/evalue=0 の5ケースを追加。
- `tests/test_scoring_engine.py`: zero weight コンポーネント、複数negative
  コンポーネントの合算後cap、の2ユニットテストを追加(41章チェックリスト
  残り2件)。
- `tests/test_scoring_engine_config.py`: `sequence_evidence` の既定値・
  カスタム上書きのテストを追加。

### Notes

- `analysis/ortholog_filter.py` と `legacy_additive` 経路
  (`_score_pair`, `_candidate_priority_score` 等)は一切変更していません。
- identity ceiling(90.0)・evalue参照下限(1e-100)は実データでの校正
  裏付けがない暫定値です(他のv2の重み・cap値と同様)。

## 未リリース: interaction_scoring v2 (evidence-based scoring)

`ProteinHunter_v5 × ProteinInteractionHunter 統合設計書 v1.0` に基づく、相互作用スコアリングの改修。
既存の `legacy_additive` 方式はそのまま残り、デフォルトで有効(挙動は変更なし)。

### Added

- `core/evidence.py`: `EvidenceStatus`(`AVAILABLE`/`MISSING`/`NOT_RUN`/`NOT_APPLICABLE`/`FAILED`/`MALFORMED`/`EXCLUDED`)と `EvidenceComponent` を追加。「評価したが証拠なし」と「評価できなかった」を区別する土台。
- `analysis/scoring_engine.py`: カテゴリ上限(category cap)付きの正規化スコアリングエンジンを追加。欠損証拠を分母から除外し、二重加点を抑制し、Negative Evidenceを上限付きペナルティとして分離し、Evidence Tier(`Tier1_VeryStrong`〜`Tier4_Weak`/`Unclassified`)と決定論的ランキングを算出。
- `analysis/scoring_engine_config.py` + `config/scoring_engine.example.yaml`: スコアリングエンジンの重み・カテゴリ上限・ペナルティ・Tier閾値を設定ファイル化(コードへのハードコード禁止)。
- `analysis/functional_complementarity_rules.py` + `config/functional_complementarity_rules.v1.yaml`: これまで `COMPLEMENTARY_TERM_PAIRS` としてコード内に固定していたキーワード対応表をバージョン管理可能なYAMLへ外部化。
- `analysis/interaction_scoring.py`: `interaction_scoring.scoring_model: v2_evidence_based` を指定すると、上記エンジンを使った新しいスコアリング経路(`_score_pair_v2` 系)が有効になる。`legacy_additive`(デフォルト)は変更なし。
- `analysis/interaction_scoring.py`: `ortholog_filter.py` が既に計算している negative BLAST hit 強度(`record.negative_hit_strength`)を `negative_hit_strength` エビデンスコンポーネント(カテゴリ `source_reliability`、`is_negative=True`)としてv2エンジンに接続。negative hit がない場合は `NOT_APPLICABLE`(ペナルティなし)、ある場合はスコアから上限付きで減算されます。`ortholog_filter.py` 自体は無変更(既存の分類結果を読むだけ)。
- `analysis/pih_evidence_bridge.py`: ProteinInteractionHunter(PIH、別リポジトリ・コード非依存)が出力する `candidate_evidence_bundle.jsonl` をプレーンなJSONとして読み込む、ファイルベースの任意連携ブリッジを追加。PIHのコードは一切importしません。PIHの5カテゴリのうち、v5が既に自前で計算している `genomic_context`/`functional_annotation` は二重加点を避けるため意図的に除外し、v5に相当機能がない `cellular_compatibility`/`evolutionary`/`direct_interaction` の3カテゴリのみを `pih_*` 名で取り込みます。`interaction_scoring.pih_evidence_bundle` で有効化(未設定なら従来どおり無効)。ファイルが存在しない・壊れている場合も実行を止めず、warningとして記録します。
- `docs/scoring_engine_v2.md`: v2スコアリングの設計意図と使い方をまとめたドキュメントを追加。negative_hit_strength とPIHブリッジの節を追記。
- `tests/test_evidence.py`, `tests/test_scoring_engine.py`, `tests/test_scoring_engine_config.py`, `tests/test_functional_complementarity_rules.py`, および `tests/test_interaction_scoring.py` へのv2統合テストを追加(negative evidence、PIHブリッジの正常系・異常系を含む)。

### Notes

- `analysis/ortholog_filter.py`(negative BLAST hit強度分類・homolog除去)は一切変更していません。
- v2モードでも `alphafold_readiness_score` は参考情報として出力されますが、合計スコアには加算されません(構造予測しやすさは相互作用の証拠ではないため)。
- v2の重み・カテゴリ上限・Tier閾値は、旧方式の配点を踏襲した暫定値です。正解データ(既知の相互作用/非相互作用ペア)による校正は未実施です。PIHブリッジのカテゴリ重みも同様に未校正です。
- セルフレビューで発見・修正: `config/scoring_engine.example.yaml` が `DEFAULT_CATEGORY_CAPS` の `pih_*` カテゴリ追加後も更新されておらず、そのままコピーして使うと実データでPIHブリッジのカテゴリが初めて発火した瞬間に `ConfigError` で落ちる状態でした。サンプルファイルを修正し、`tests/test_scoring_engine_config.py::test_example_config_matches_defaults` で今後の再発を防止しています。

## 開発版 v0.1.0

現在の最初の開発版です。FASTA 入力から BLAST 候補探索、注釈、スコアリング、Excel 出力までの最小パイプラインが入りました。

### Added

- 起動前チェックを追加
- 色付きログ表示と実行サマリーを追加
- `config.yaml` の読み込みと設定値の検証を追加
- FASTA 解析ユーティリティを追加
- 入力 FASTA の件数サマリーをログに出す機能を追加
- BLAST データベース作成と `blastp` 実行ユーティリティを追加
- BLAST による positive-only 候補抽出パイプラインを追加
- `ProteinRecord`、`BlastHit`、`DomainHit`、`CandidateScore` モデルを追加
- JSON ファイルを使うキャッシュ機能を追加
- CDD 注釈ユーティリティを追加
- Pfam 注釈ユーティリティを追加
- UniProt 注釈ユーティリティを追加
- AlphaFold URL 注釈ユーティリティを追加
- 候補スコアリング機能を追加
- Excel 出力機能を追加
- Excel の見やすさを整える書式設定を追加
- `--config` によるデモ設定ファイル指定を追加
- `--check-only` による事前確認モードを追加
- デモ用 FASTA ファイルを本番入力とは別に配置
- README、data README、本番実行チェックリストなどのドキュメントを追加

### Changed

- デモ FASTA ファイルを本番用の入力パスから分離
- 本番用の `config.yaml` とデモ用の `config.demo.yaml` を分離
- UniProt 注釈と AlphaFold 注釈の処理を分離
- `annotation.enable_cdd`、`annotation.enable_pfam`、`annotation.enable_uniprot`、`annotation.enable_alphafold` の有効・無効設定を実行時に反映

### Fixed

- CDD/Pfam パーサーが小数のスコアを座標として誤読しないよう修正
- 不正な FASTA ファイルのエラーを `FileValidationError` として分かりやすく扱うよう修正
- 生成される BLAST 一時ファイルが Git で追跡されないよう整理

### Notes

- 本番実行には、次の実データ FASTA が必要です。
  - `data/input/target.faa`
  - `data/databases/positive.faa`
  - `data/databases/negative.faa`
- デモ実行では、次を使います。
  - `config.demo.yaml`
  - `data/demo/*.faa`
- フルテストは `139 passed` で通過済みです。
