# 本番実行チェックリスト

ProteinHunter_v5 を本番データで実行する前に、このチェックリストを上から順番に確認してください。

## 1. 実行前の確認

- [ ] WSL Ubuntu など、Linux/WSL 環境で作業している
- [ ] プロジェクトの `.venv` を有効化できる
- [ ] `blastp` と `makeblastdb` が使える
- [ ] Python 依存パッケージがインストール済み
- [ ] `config.yaml` の内容が現在の解析に合っている
- [ ] `config.yaml` の値が空欄や不正な値になっていない
- [ ] `config.yaml` を `config.demo.yaml` で置き換えていない

確認例:

```bash
which blastp
which makeblastdb
.venv/bin/python -m pytest tests -q
```

## 2. 入力ファイル

次の 3 つの FASTA ファイルを確認してください。

- [ ] `data/input/target.faa`
  - 解析したいタンパク質配列を入れるファイル
  - 候補探索の対象になる配列を入れます

- [ ] `data/databases/positive.faa`
  - 目的の機能に近いと考える参照タンパク質配列を入れるファイル
  - positive BLAST の参照として使われます

- [ ] `data/databases/negative.faa`
  - 除外したい、または目的の機能とは違う参照タンパク質配列を入れるファイル
  - negative BLAST の参照として使われます

注意:

- [ ] デモ用やテスト用 FASTA が残っていない
- [ ] 本番解析用の FASTA に置き換わっている
- [ ] 各 FASTA に少なくとも 1 件以上の配列が入っている
- [ ] `data/demo/` の FASTA を本番入力として誤って使っていない

デモ FASTA は `data/demo/` に保存されています。動作確認には使えますが、本番解析では必ず `data/input/target.faa`、`data/databases/positive.faa`、`data/databases/negative.faa` を実データで用意してください。

デモ用ファイル:

- `data/demo/target_demo.faa`
- `data/demo/positive_demo.faa`
- `data/demo/negative_demo.faa`

本番解析では、これらのデモ FASTA を使わないでください。

## 3. `config.yaml` のおすすめ確認項目

- [ ] `paths.target_fasta`
- [ ] `paths.positive_fasta`
- [ ] `paths.negative_fasta`
- [ ] `paths.output_excel`
- [ ] `blast.evalue`
  - 正の数になっている
- [ ] `blast.max_target_seqs`
  - 正の整数になっている
- [ ] `blast.threads`
  - `auto` または正の整数になっている
- [ ] `annotation.enable_cdd`
  - CDD 注釈を使うなら `true`
- [ ] `annotation.enable_pfam`
  - Pfam 注釈を使うなら `true`
- [ ] `annotation.enable_uniprot`
  - UniProt 注釈を使うなら `true`
- [ ] `annotation.enable_alphafold`
  - AlphaFold 注釈を使うなら `true`

## 4. 実行方法

### 本番解析

本番解析では、標準で `config.yaml` が使われます。
WSL 側のターミナルで実行します。

```bash
cd /mnt/c/Users/nyako/Documents/GitHub/ProteinHunter_v5
source .venv/bin/activate
.venv/bin/python main.py
```

### デモ実行

動作確認だけをしたい場合は、`config.demo.yaml` を指定します。
この実行では `data/demo/` のデモ FASTA が使われ、結果は `data/output/ProteinHunter_demo_results.xlsx` に出力されます。

```bash
cd /mnt/c/Users/nyako/Documents/GitHub/ProteinHunter_v5
source .venv/bin/activate
.venv/bin/python main.py --config config.demo.yaml
```

本番解析とデモ実行は分けて扱ってください。`config.yaml` を `config.demo.yaml` に置き換えないでください。

## 5. ログで確認すること

実行中または実行後に、画面表示や `logs/latest.log` を確認してください。

- [ ] Startup check passed
- [ ] Input FASTA summary の件数
  - Target proteins
  - Positive references
  - Negative references
- [ ] BLAST positive-only candidate count
- [ ] CDD domain hit count
- [ ] Pfam domain hit count
- [ ] UniProt annotation count
- [ ] AlphaFold annotation count
- [ ] Top candidate
- [ ] Top candidate score
- [ ] Excel output path

## 6. Excel で確認すること

標準設定では、結果は次に出力されます。

```text
data/output/ProteinHunter_results.xlsx
```

Excel で特に見る列:

- [ ] `total_score`
- [ ] `score_components`
- [ ] `score_reasons`
- [ ] `domain_sources`
- [ ] `domain_names`
- [ ] `domain_accessions`
- [ ] `domain_descriptions`
- [ ] `uniprot_accession`
- [ ] `alphafold_url`
- [ ] `notes`

`notes` に注釈失敗やスキップ理由が書かれている場合があります。高スコア候補でも、`notes` は必ず確認してください。

## 7. GitHub にコミットしない生成ファイル

次のディレクトリは実行時に生成されるため、通常は GitHub にコミットしません。

- [ ] `data/temp/`
- [ ] `data/output/`
- [ ] `logs/`
- [ ] `.cache/`

特に `data/temp/` には BLAST の一時ファイルや BLAST データベースが入ります。これらは再生成できるため、GitHub で追跡しないようにしてください。
