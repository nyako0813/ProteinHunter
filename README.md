# ProteinHunter_v5

ProteinHunter_v5 は、タンパク質 FASTA を入力として、BLAST、ドメイン注釈、外部注釈、簡易スコアリングを行い、候補一覧を Excel に出力する解析ツールです。

## 現在の処理の流れ

1. 起動前チェック
2. `config.yaml` の読み込みと検証
3. 入力 FASTA の件数サマリー表示
4. BLAST による positive-only 候補探索
5. CDD 注釈（有効な場合）
6. Pfam 注釈（有効な場合）
7. UniProt 注釈（有効な場合）
8. AlphaFold 注釈（有効な場合）
9. 候補スコアリング
10. Excel 出力

## 必要な入力ファイル

次の 3 つの FASTA ファイルを用意してください。

- `data/input/target.faa`
  - 解析したいタンパク質配列を入れます。
- `data/databases/positive.faa`
  - 目的の機能に近いと考える参照タンパク質配列を入れます。
- `data/databases/negative.faa`
  - 除外したい、または目的の機能とは違う参照タンパク質配列を入れます。

テスト用やデモ用の FASTA が入っている場合は、本番解析の前に実際の解析用ファイルへ置き換えてください。

## 出力ファイル

標準設定では、結果は次の Excel ファイルに出力されます。

- `data/output/ProteinHunter_results.xlsx`

Excel には候補 ID、説明、スコア、BLAST 結果、CDD/Pfam ドメイン情報、UniProt/AlphaFold 情報、メモなどが入ります。

## 主な config.yaml 設定

### `paths`

入力 FASTA、出力 Excel、キャッシュ、ログの場所を指定します。

特に重要な項目です。

- `paths.target_fasta`
- `paths.positive_fasta`
- `paths.negative_fasta`
- `paths.output_excel`
- `paths.cache_dir`
- `paths.log_dir`

### BLAST 設定

- `blast.evalue`
  - BLAST の e-value しきい値です。正の数を指定します。
- `blast.max_target_seqs`
  - BLAST で取得する最大ヒット数です。正の整数を指定します。
- `blast.threads`
  - 使用スレッド数です。`auto` または正の整数を指定します。

### 注釈の有効・無効

次の値を `true` / `false` で切り替えます。

- `annotation.enable_cdd`
- `annotation.enable_pfam`
- `annotation.enable_uniprot`
- `annotation.enable_alphafold`

無効にした注釈ステップはスキップされます。スキップしても、後続のスコアリングと Excel 出力は続行されます。

## 実行方法

WSL / Linux 環境で、プロジェクトのルートディレクトリから実行します。

### 本番解析

通常の本番解析では、何も指定しなければ `config.yaml` が使われます。
本番解析を行う前に、次の 3 つの FASTA ファイルを実データで用意してください。

- `data/input/target.faa`
- `data/databases/positive.faa`
- `data/databases/negative.faa`

```bash
source .venv/bin/activate
.venv/bin/python main.py
```

### デモ実行

動作確認だけをしたい場合は、`config.demo.yaml` を指定します。
この設定では、次のデモ FASTA が使われます。

- `data/demo/target_demo.faa`
- `data/demo/positive_demo.faa`
- `data/demo/negative_demo.faa`

```bash
source .venv/bin/activate
.venv/bin/python main.py --config config.demo.yaml
```

デモ実行の Excel 出力先は次のファイルです。

- `data/output/ProteinHunter_demo_results.xlsx`

デモ FASTA は本番解析には使わないでください。
また、`config.yaml` を `config.demo.yaml` で置き換えず、本番用とデモ用の設定を分けて管理してください。

## 生成されるファイル

次のディレクトリやファイルは実行時に生成されます。

- `data/temp/`
- `data/output/`
- `logs/`
- `.cache/`

これらは解析ごとに変わるため、通常は GitHub にコミットしません。

特に `data/temp/` には BLAST の一時ファイルや作成済みデータベースが入ります。生成された BLAST temp ファイルは GitHub で追跡しないようにしてください。
