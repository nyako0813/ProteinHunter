# CHANGELOG

ProteinHunter_v5 の変更履歴です。

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
