# CHANGELOG

ProteinHunter_v5 の変更履歴です。

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
- フルテストは `137 passed` で通過済みです。
