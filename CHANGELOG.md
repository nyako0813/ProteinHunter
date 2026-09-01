# CHANGELOG

ProteinHunter_v5 の変更履歴です。

## 未リリース: interaction_scoring v2 (evidence-based scoring)

`ProteinHunter_v5 × ProteinInteractionHunter 統合設計書 v1.0` に基づく、相互作用スコアリングの改修。
既存の `legacy_additive` 方式はそのまま残り、デフォルトで有効(挙動は変更なし)。

### Added

- `core/evidence.py`: `EvidenceStatus`(`AVAILABLE`/`MISSING`/`NOT_RUN`/`NOT_APPLICABLE`/`FAILED`/`MALFORMED`/`EXCLUDED`)と `EvidenceComponent` を追加。「評価したが証拠なし」と「評価できなかった」を区別する土台。
- `analysis/scoring_engine.py`: カテゴリ上限(category cap)付きの正規化スコアリングエンジンを追加。欠損証拠を分母から除外し、二重加点を抑制し、Negative Evidenceを上限付きペナルティとして分離し、Evidence Tier(`Tier1_VeryStrong`〜`Tier4_Weak`/`Unclassified`)と決定論的ランキングを算出。
- `analysis/scoring_engine_config.py` + `config/scoring_engine.example.yaml`: スコアリングエンジンの重み・カテゴリ上限・ペナルティ・Tier閾値を設定ファイル化(コードへのハードコード禁止)。
- `analysis/functional_complementarity_rules.py` + `config/functional_complementarity_rules.v1.yaml`: これまで `COMPLEMENTARY_TERM_PAIRS` としてコード内に固定していたキーワード対応表をバージョン管理可能なYAMLへ外部化。
- `analysis/interaction_scoring.py`: `interaction_scoring.scoring_model: v2_evidence_based` を指定すると、上記エンジンを使った新しいスコアリング経路(`_score_pair_v2` 系)が有効になる。`legacy_additive`(デフォルト)は変更なし。
- `docs/scoring_engine_v2.md`: v2スコアリングの設計意図と使い方をまとめたドキュメントを追加。
- `tests/test_evidence.py`, `tests/test_scoring_engine.py`, `tests/test_scoring_engine_config.py`, `tests/test_functional_complementarity_rules.py`, および `tests/test_interaction_scoring.py` へのv2統合テストを追加。

### Notes

- `analysis/ortholog_filter.py`(negative BLAST hit強度分類・homolog除去)は一切変更していません。
- v2モードでも `alphafold_readiness_score` は参考情報として出力されますが、合計スコアには加算されません(構造予測しやすさは相互作用の証拠ではないため)。
- v2の重み・カテゴリ上限・Tier閾値は、旧方式の配点を踏襲した暫定値です。正解データ(既知の相互作用/非相互作用ペア)による校正は未実施です。

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
