# 実行の再現性・トレーサビリティ(run provenance) — 設計書

Status: 実装済み(実装メモは末尾)。`main`(PR #15まで)のコードを読んで
integration pointを確認済み。2026-09-20: スコープ拡張の動機となった
実例([[negative_reference_species_mislabeling_finding]])は、この設計とは
別に`core/reference_genomes.py`として先に実装・PR #17でpush済み(下記
「スコープ拡張の検討」参照、対応方針が確定済み)。

## 背景・問題

このプロジェクトは opt-in の切り替えが非常に多い
(`scoring_model`, `candidate_sources` の9バケツ, `string_ppi_ncbi_taxon_id`,
`geo_coexpression_enabled`, `pih_evidence_bundle`,
`consider_cross_species_matches`, `ranking_metric`, 各種cap/weight...)。
git commit単位でコードは追跡できるが、

- `config.yaml` はgit管理下にあるとは限らず、実験的に手元だけで書き換えて
  実行することが日常的にありうる(このプロジェクトの性質上、複数マシン
  ([[implementation_status_scoring_v2]]記載のWSL/Windows/ノートPC)に
  またがって作業しており、`config.yaml` 自体の差分をコミットせずに
  試行錯誤するケースは十分想定される)。
- 現状、生成される `ProteinHunter_results_*.xlsx` / Word報告書のどこにも、
  「どのconfig設定で」「どのコード版で」生成したかの記録が一切ない。
  `output/word_report.py::write_word_report` は現在
  "Report generated: {datetime}" と "Scoring model: {scoring_model}" の
  2行しかタイトルページに書いていない(該当箇所: line 482-483)。
- 実際に過去起きた事故: Phase 6-8 Stage 1の実データ検証で「分離幅 +32.89」という誤った数値が一度報告され、後で「STRING/GEOを無効化した
  設定で計測していた」ことが判明して自己修正している
  (`claude/implementation_status_scoring_v2.md` の該当節)。これは
  まさに「どの設定で出した数値か」を実行結果自体からは追跡できなかった
  ために起きた誤りであり、この機能があれば再発を防げた具体例。
- もう1件、2026-09-20に実際に起きた事故として追加できる: PR #20の
  GFF既定変更の前後比較で、最初の比較実行がGFF注釈を2バケツにしか
  適用しない設定のまま行われ、STRING/GEOエビデンスが歪んだ状態で比較
  してしまい、破棄してやり直す事態になった
  ([[negative_reference_species_mislabeling_finding]]参照)。これも
  「その実行が実際にどのconfig設定で行われたか」がExcel/Word側から
  追跡できないために起きた誤りで、run provenance機能があれば「あの比較、
  どの設定だっけ」を都度configファイルを読み返さずに確認できていた。

## 採用する方針

2段構えにする。

1. **常設の短い breadcrumb**: Excel `01_Index` シートとWordレポートの
   タイトルページに、git commit SHA(短縮)・dirtyフラグ・configの
   内容ハッシュ(短縮)を必ず1行ずつ追加する。実行のたびに自動で入る、
   人間が読む用の「指紋」。
2. **サイドカーファイル(常時出力)**: 出力Excel/Wordと同じディレクトリに、
   実行に使われた実効設定の全体(`Config` dataclass全体、CLI引数や
   デフォルト適用後の最終形)をYAMLとしてダンプする
   `<output_excelのstem>.run_provenance.yaml` を必ず書き出す。
   commitしていないconfig.yamlの一時的な書き換えでも、出力ファイルさえ
   保存しておけば後から正確に再現できるようにする、実際のペイロード。

短いbreadcrumbだけでは「gitにコミットしていない設定変更」を復元できず、
サイドカーだけでは目視で気づけない(ファイルを開かないと分からない)。
両方揃えて初めて実用的な再現性になる。

## 実装内容

### 1. 新規モジュール `core/provenance.py`

```python
@dataclass(frozen=True)
class RunProvenance:
    app_version: str          # core.constants.APP_VERSION
    git_commit: str | None    # 短縮SHA(7桁)。gitリポジトリでない/取得失敗時 None
    git_dirty: bool | None    # 追跡ファイルに未コミット差分があるか。取得不能なら None
    config_hash: str          # 下記「configハッシュの対象」参照。短縮8桁
    reference_genome_check_passed: bool | None  # 下記「スコープ拡張の検討」参照
    generated_at: datetime

def collect_run_provenance(config: Config, repo_root: Path) -> RunProvenance: ...
def write_provenance_sidecar(provenance: RunProvenance, config: Config, output_path: Path) -> Path: ...
```

- `git_commit`/`git_dirty` の取得は
  `subprocess.run(["git", "rev-parse", "--short=7", "HEAD"], cwd=repo_root, ...)`
  と `git status --porcelain` を使う(新規の外部依存を増やさない —
  GitPython 等は導入しない)。`git` コマンド自体が無い/`.git` が無い
  (配布zip等)場合は例外を握りつぶして `None` にする、best-effort。
  ログには warning ではなく info レベルで「git provenance unavailable」
  とだけ出す(git無しでも実行は完全に正常なので、警告は過剰)。
- configハッシュの対象: `config.yaml` 本体に加え、スコアリング挙動を
  左右する2つの外部YAML — `scoring_engine_config`(`config/scoring_engine.example.yaml`
  相当)と `functional_complementarity_ruleset`
  (`config/functional_complementarity_rules.v1.yaml`) — も
  実際に使われたパスのファイル内容を連結してSHA-256を取る。
  `pih_evidence_bundle`/STRINGキャッシュ等の「データ」ファイルは対象外
  (これらはconfigではなくデータであり、パス+更新時刻程度の軽い記録に
  留める。フルハッシュ化するとGB級のキャッシュを毎回読むことになり
  実行時間への影響が無視できない)。陰性/陽性参照ゲノムについては
  下記「スコープ拡張の検討」の通り、`core/reference_genomes.py`の
  起動時チェック結果を再利用する形にし、ここでの独自ハッシュ化は
  行わない。
- サイドカーの中身: `dataclasses.asdict(config)` をそのままYAML
  ダンプする。`Path` 型フィールドは `str` に変換する adapter が必要
  (`yaml.safe_dump` はデフォルトで`pathlib.Path`を扱えないため、
  `default=str` 相当の変換をシリアライズ前に噛ませる)。

### 2. `main.py` への配線

- Excel/Word両方のwriterを呼ぶ直前、1回だけ
  `provenance = collect_run_provenance(config, PROJECT_ROOT)` を計算し、
  `write_classification_workbook(...)` と `write_word_report(...)` の
  両方に `provenance=provenance` として渡す(両関数のシグネチャに
  optional引数を1つ追加。デフォルト `None` で後方互換 —
  例えばこの2関数を直接呼ぶ既存テストが引数無しでも壊れない)。
- 実行の最後に
  `write_provenance_sidecar(provenance, config, config.paths.output_excel)`
  を呼ぶ。Word/Excelどちらかしか有効でない
  実行でも、`output_excel` は必須設定なので基準パスとして使える。
- `core/reference_genomes.py`の起動時チェック(既にPR #17で`main.py`に
  配線済み)の結果(pass/fail、警告件数)を`collect_run_provenance`に
  渡し、`RunProvenance.reference_genome_check_passed`に格納する
  (チェック自体の再実行はしない、結果を受け取るだけ)。

### 3. Excel側 (`output/excel.py`)

`01_Index` シートの `INDEX_ROWS_V2` テーブルの上、既存の
`_format_index_worksheet` が書き込むヘッダより前の数行に、
provenanceの3行(app_version/git/config_hash + generated_at)を
追加する小さな `_write_provenance_header(worksheet, provenance)` を
新設。既存のテーブル書き込み開始行をそのぶん下にずらす(既存の
座標決め打ちがあれば、offset変数化が必要)。
`reference_genome_check_passed`が`False`の場合は、この見出しに
"Reference genome check: WARNING (see log)"のような1行を追加する
(詳細はログ/`core/reference_genomes.py`の警告に譲り、ここでは
気づけることだけを目的にする)。

### 4. Word側 (`output/word_report.py`)

`write_word_report` のタイトルページ部分(line 481-488)に3〜4行追加する:

```python
document.add_paragraph(f"Report generated: {datetime.now():%Y-%m-%d %H:%M}")
document.add_paragraph(f"Scoring model: {scoring_model}")
if provenance is not None:
    document.add_paragraph(
        f"Code version: {provenance.app_version}"
        f" (git {provenance.git_commit or 'unknown'}"
        f"{'+dirty' if provenance.git_dirty else ''})"
    )
    document.add_paragraph(f"Config fingerprint: {provenance.config_hash}")
    document.add_paragraph(
        f"Full effective configuration saved alongside this report as "
        f"{sidecar_filename}"
    )
    if provenance.reference_genome_check_passed is False:
        document.add_paragraph(
            "Reference genome check: WARNING — see run log for details "
            "(core/reference_genomes.py)"
        )
```

### 5. サイドカーファイルへの参照を相互に埋め込む

Excel/Wordどちらのprovenance行にも、サイドカーファイル名を明記する
(「このExcelを見ている人が、対応する完全な設定をすぐ見つけられる」ため)。
Excel↔Word相互リンクの設計(PR #14)と同じ思想 — ファイル名ベースの
クロスリファレンスに倣う。

## スコープ外(あえてやらないこと)

- STRING/GEOキャッシュ・PIHバンドル・BLASTデータベース自体のバージョン
  管理(NCBI側のデータ更新でスコアが変わりうる、という別種の再現性問題)
  は今回は扱わない。将来必要になれば「データソース鮮度」の記録として
  別途設計する。
- 参照ゲノムディレクトリの構成チェックそのものは実装しない
  (`core/reference_genomes.py`がPR #17で既に実装済み。本設計は
  その結果を1フィールドとして取り込むだけで、ロジックを重複させない)。
- 過去に生成済みの出力ファイルへの遡及的なprovenance付与はしない
  (対象は「これから生成する出力」のみ)。
- config復元用のCLIコマンド(サイドカーYAMLを読んで実際に
  `config.yaml` を上書きする機能)は今回は作らない。サイドカーは
  「人間が見て手で反映する」ドキュメントとして扱う。自動復元は
  誤って本番configを上書きする事故のリスクの方が大きいと判断。

## テスト計画

- `tests/test_provenance.py`(新規): git無し環境(一時ディレクトリを
  `cwd` にしてgitリポジトリでない状態を作る)で
  `collect_run_provenance` が例外を出さず `git_commit=None` を返すこと。configハッシュが
  ファイル内容の変更で変わり、同一内容では再現すること(決定的)。
  `reference_genome_check_passed`を明示的に渡した場合にそのまま
  `RunProvenance`に反映されることの確認も含める。
- `tests/test_excel_output.py` / `tests/test_word_report.py` に
  provenance行がIndex/タイトルページに出ることの確認を追加。
  `provenance=None` を渡した場合(後方互換パス)に既存の出力が
  変わらないことも回帰テストとして残す。`reference_genome_check_passed
  =False`のときに警告行が出ることの確認も追加。
- 実データ検証: 実際に `config.yaml` を1箇所書き換えて2回実行し、
  config_hashが変わること、gitにコミットしていない変更でも
  サイドカーYAMLの中身が変更後の値になっていることを目視確認する。

## 見積もり

新規モジュール1つ(小)、Excel/Word writerへの数行ずつの追加、
`main.py` の配線数行。既存の `output/excel.py`/`output/word_report.py`
のシグネチャ変更は後方互換なdefault引数なので、既存呼び出し元
(テスト含む)への影響は最小。1PRで完結できる規模。

## スコープ拡張の検討(2026-09-20追記、対応方針確定)

[[negative_reference_species_mislabeling_finding]]で、`data/databases/
negative/Sulfolobus_solfataricus/`に実際にはMethanococcus maripaludis
のデータが入っているという取り違えが発覚した。この件でもう1つ判明した
のが、`data/databases/negative/thermoplasma_acidophilum/`が
`.gitignore`で除外されており、除外理由の記録が残っていないという点。
別マシンでこのリポジトリをcloneすると、陰性参照ゲノムが(意図した3種
ではなく)2種になり、エラーも警告も出ないまま異なる分類結果になる、
という再現性の抜け穴だった。

**対応済み(2026-09-20、PR #17): この抜け穴自体は、run provenance機能を
待たずに`core/reference_genomes.py` + `config/reference_genomes.v1.yaml`
+ `main.py`起動時チェックとして先に実装・push済み。** 期待する陽性・
陰性ゲノムのフォルダ名/accession/生物種のマニフェストと実データを
突き合わせ、フォルダ欠落・リスト外フォルダ・accession不一致・FASTA
生物種不一致・別フォルダ間の完全一致・陽性陰性間の同一種混入の6種を
検出しWARNINGログを出す。`tests/test_reference_genomes.py`でテスト済み。

したがって、本設計(run_provenance)を実装する際は、この機能を
再実装・重複させず、その結果(チェックが通ったか)を`RunProvenance`の
1フィールド(`reference_genome_check_passed`)として取り込むだけに
留める。
上記「実装内容」1〜4節に、この統合方法を反映済み。

## 実装順序に関する注意(2026-09-20追記)

現時点でPR #17(ゲノム修正)・#18(クローズ済み)・#19(ドキュメント)・
#20(GFF既定変更)がマージ待ち。本設計は`output/excel.py`・
`output/word_report.py`・`main.py`を触るが、[[conserved_query_visibility_design]]
(PR #16)も`output/word_report.py`(5.7節)と`main.py`の警告出力経路に
手を入れている。ファイル内の変更箇所は別セクションなので直接の
コンフリクトは想定しにくいが、実装着手はPR #16〜#20が全てmainに
マージされた後にすること。 特に`reference_genome_check_passed`の
統合は、`core/reference_genomes.py`(PR #17)がmainに入っていることが
前提。

## 実装メモ(2026-09-20、実装時に追記)

PR #16〜#20のマージ後に実装。設計からの差分と判断:

- **`collect_run_provenance` のシグネチャ**: `config_path`・
  `reference_genome_check_passed`・`generated_at` をキーワード引数で追加した。
  `Config` は読み込んだ設定ファイルのパスを保持しておらず、configハッシュには
  実際に使った `config.yaml` の内容が要るため。
- **configハッシュの対象に `domain_family_map_path` も追加**。設計の原則
  (「スコアリング挙動を左右する外部YAML」)に照らすと、
  `domain_complementarity` を決める `config/domain_family_map.v1.yaml` も
  該当するため。パスが未設定(組み込み既定)・ファイル欠落・読み取り不能は
  それぞれ別のプレースホルダとしてハッシュに入り、状態の違いが区別される。
- **CRLF→LF 正規化してからハッシュ**。複数マシン(WSL/Windows)で同一内容の
  ファイルが改行違いで別ハッシュになるのを防ぐ。
- **git のトップレベル一致を確認**: `git rev-parse --show-toplevel` が
  `repo_root` と一致しない場合(別リポジトリの中に展開した配布zip等)は
  そのリポジトリのSHAを拾わず `None` にする。`git status` は
  `--untracked-files=no`(追跡ファイルの差分のみ、設計どおり)。
- **`reference_genome_check_passed` の意味**: `main._log_reference_genome_check` が
  `True`(マニフェストと一致)/`False`(指摘あり)/`None`(検証できず: マニフェスト
  無し+指摘無し、チェック自体が失敗、ファイル入力モード)を返す。警告行は
  `False` のときだけ出る。
- **サイドカー**: 設計の `dataclasses.asdict` の代わりに、Path・tuple・Enum・
  datetimeを扱う再帰変換(`_to_plain`)を使った。ペイロードは
  `run_provenance`(指紋)・`config_file`(使ったconfigのパス)・
  `effective_config`(実効設定全体)。サイドカー書き出しの失敗は
  出力が既に書かれた後なので警告ログのみで、実行は失敗させない。
- **Excel の Index**: 指紋ブロック(4〜5行)+空行の後にテーブルを書く
  (`to_excel(startrow=...)`)。`_format_index_worksheet` は `header_row`
  引数でリンクとフォーマットの開始行をずらす。`provenance=None` のときは
  従来どおり1行目がテーブルのヘッダ。
