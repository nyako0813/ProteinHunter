# data ディレクトリの使い方

このディレクトリには、ProteinHunter が解析に使う入力ファイルと、実行時に作られる出力ファイルを置きます。

## デモ用 FASTA

### `data/demo/`

`data/demo/` には、動作確認用の小さな FASTA ファイルが入っています。

- `data/demo/target_demo.faa`
- `data/demo/positive_demo.faa`
- `data/demo/negative_demo.faa`

これらはテストや操作確認のためのデモデータです。本番解析では使わないでください。

デモ FASTA で動作確認する場合は、`config.demo.yaml` を指定して実行します。

```bash
source .venv/bin/activate
.venv/bin/python main.py --config config.demo.yaml
```

デモ実行の結果は次に出力されます。

- `data/output/ProteinHunter_demo_results.xlsx`

`config.demo.yaml` はデモ専用です。本番用の `config.yaml` と置き換えないでください。

デモ FASTA の件数確認だけをしたい場合は、`--check-only` を使います。

```bash
source .venv/bin/activate
.venv/bin/python main.py --config config.demo.yaml --check-only
```

この確認では、BLAST、注釈、スコアリング、Excel 出力は行われません。

## 必要な入力ファイル

### `data/input/target.faa`

解析したいタンパク質配列を入れる FASTA ファイルです。

ここに入れた各配列が、ProteinHunter の候補探索の対象になります。たとえば、未機能解析タンパク質や、機能を調べたいタンパク質群を入れます。

本番解析では、ユーザーが実際の解析対象 FASTA をこの場所に用意してください。

### `data/databases/positive.faa`

「近い機能を持つ可能性がある」と考える参照タンパク質配列を入れる FASTA ファイルです。

ProteinHunter は、このファイルを positive BLAST の参照として使います。対象タンパク質がここに含まれる配列と似ている場合、候補として前向きな根拠になります。

本番解析では、ユーザーが目的に合った positive 参照 FASTA をこの場所に用意してください。

### `data/databases/negative.faa`

候補から除外したい、または目的の機能とは違うと考える参照タンパク質配列を入れる FASTA ファイルです。

ProteinHunter は、このファイルを negative BLAST の参照として使います。対象タンパク質がここに含まれる配列と似ている場合、目的とは違う候補として扱われる可能性があります。

本番解析では、ユーザーが目的に合った negative 参照 FASTA をこの場所に用意してください。

## 生成されるファイル

### `data/temp/`

BLAST データベースや中間ファイルなど、解析中に自動生成される一時ファイルが入ります。

この中身は再生成できるため、通常は Git にコミットしないでください。

### `data/output/`

Excel レポートなど、解析結果の出力ファイルが入ります。

解析のたびに内容が変わるため、通常は Git にコミットしないでください。

## 本番解析の前に確認すること

通常の本番解析では、`config.yaml` が使われます。
まず、次の 3 つの FASTA を実データで用意してください。

- `data/input/target.faa`
- `data/databases/positive.faa`
- `data/databases/negative.faa`

次に、`--check-only` で入力件数を確認します。

```bash
source .venv/bin/activate
.venv/bin/python main.py --check-only
```

表示された入力 FASTA の件数が正しいことを確認してから、本番解析を実行します。

```bash
source .venv/bin/activate
.venv/bin/python main.py
```

テスト用やデモ用の FASTA ファイルが入っている場合は、本番解析の前に必ず実際の解析用 FASTA ファイルへ置き換えてください。

特に `target.faa`、`positive.faa`、`negative.faa` の内容が目的に合っているかを確認してから実行してください。

`data/demo/` の FASTA を誤って本番解析に使わないように注意してください。

## Directory-based production inputs

For production runs, you may place NCBI-downloaded folders directly under the database directories instead of manually combining FASTA files.

Expected structure:

```text
data/databases/
  target/
    Organism_A/
      ncbi_dataset/
        data/
          GCF_000000001.1/
            protein.faa
  positive/
    Organism_B/
      ncbi_dataset/
        data/
          GCF_000000002.1/
            protein.faa
  negative/
    Organism_C/
      protein.faa
```

Set `input_mode: directory` and configure:

```yaml
paths:
  target_dir: "./data/databases/target"
  positive_dir: "./data/databases/positive"
  negative_dir: "./data/databases/negative"
```

ProteinHunter uses only the immediate child folders as source labels, then searches recursively inside each one for `protein.faa`. Folder names such as `Organism_A` are used as source labels in logs. Missing `protein.faa` files are skipped with a warning when at least one valid source folder exists.

Combined FASTA files are generated under `data/temp/combined/`; this directory is temporary and should not be committed.
