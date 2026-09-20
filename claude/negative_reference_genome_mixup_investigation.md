# 陰性参照ゲノム取り違え調査(Sulfolobus_solfataricus フォルダ)

Status: 修正済み(作業ツリー上、未コミット)。BLAST分類レベルの再計算まで実施、
スコアリング全体(interaction_score / rank)の再実行は未実施。

## 何が起きていたか

`data/databases/negative/Sulfolobus_solfataricus/` の中身は *Sulfolobus
solfataricus* ではなく **Methanococcus maripaludis DSM 2067
(GCF_002945325.1)** だった。しかも `data/databases/positive/Methanococcus_maripaludis/`
の `protein.faa` と**バイト単位で同一**(md5 `b1285b09...`)。
つまり陽性参照ゲノムが陰性参照にも二重に登録されていた。NCBI付属の
`md5sum.txt` とも一致するため、ファイル破損ではなく、ダウンロード時に
誤ったアセンブリを取得して置いたものと判断できる。

- 追加されたコミット: `4f02818` (2026-07-02, "Fix GFF parent old locus tag
  mapping") で Aeropyrum_pernix と同時に追加。以降内容は変更されていない
  (`git log --follow` で確認)。コミットメッセージと無関係の大きなデータ追加で、
  レビューで気づきにくい形だった。
- 同じフォルダ名を使うテスト(`tests/test_fasta_sources.py`)はtmp上の合成
  フィクスチャで実データを読まないため、検知できなかった。

## 他の参照ゲノムの確認

| フォルダ | data_summary.tsv | FASTAヘッダ | 判定 |
|---|---|---|---|
| negative/Aeropyrum_pernix | Aeropyrum pernix K1, GCF_000011125.1 | Aeropyrum pernix ×1689 | 一致 |
| negative/thermoplasma_acidophilum (未追跡) | Thermoplasma acidophilum DSM 1728, GCF_000195915.1 | Thermoplasma acidophilum ×1494 | 一致 |
| positive/Methanocaldococcus_jannaschii | M. jannaschii DSM 2661, GCF_000091665.1 | (同) | 一致 |
| positive/Methanococcus_maripaludis | M. maripaludis DSM 2067, GCA/GCF_002945325.1 | M. maripaludis | 一致 |
| negative/Sulfolobus_solfataricus (修正前) | **Methanococcus maripaludis** | **Methanococcus maripaludis ×1813** | **不一致** |

## 修正

`GCF_000007005.1`(*Saccharolobus solfataricus* P2 = 旧名 *Sulfolobus
solfataricus* P2, RefSeq reference genome, 2,859 protein-coding)を NCBI
Datasets API から取得し、`negative/Sulfolobus_solfataricus/` を差し替えた。
`md5sum.txt` の検証は全ファイルOK、FASTAヘッダは Saccharolobus solfataricus
×2541 ほか。フォルダ名(ラベル)は `[source=Sulfolobus_solfataricus]` として
他所で使われるため変更していない。API版のzipには `data_summary.tsv` が
含まれないため、`assembly_data_report.jsonl` から同じ列形式で手作成した
(整合性テストがこのファイルを読む)。旧フォルダのZone.Identifierファイルは
削除した(作業ツリー上の削除、未コミット)。BLAST DBは実行のたびに
`makeblastdb` で再構築される(`blast/runner.py`)ため、古いキャッシュが
残る心配はない。

## 影響: BLAST分類の再計算(同一config・同一target/positive、negativeのみ差し替え)

`main.py` と同じ手順(`prepare_directory_fasta` → `run_blast_classification_pipeline`)を
差し替え前後で実行。差し替え前の5クエリのstrengthは検証で実測された値
(MA_4115 strong / HdrD1 medium / NifD medium / NifK strong / MtpA weak)を再現。

バケツ件数(全 4,627 タンパク質):

| バケツ | 修正前 | 修正後 |
|---|---|---|
| Candidates(positive-only) | 127 | **709** |
| Candidates_relaxed | 1261 | 1668 |
| Negative_hit | 2253 | 1737 |
| Negative_strong_hit | 834 | 435 |
| Negative_medium_hit | 593 | 558 |
| Negative_weak_hit | 493 | 454 |
| Negative_unmatched | 2374 | 2890 |

strength 内訳: none 2707→3180 / weak 493→454 / medium 593→558 / strong 834→435。
陽性ゲノムの相同体を持つタンパク質は、陰性側の M. maripaludis コピーにも強くヒットして
`strong_only` 除外に掛かっていた。**Candidates が約5.6倍になる**のはこのため。

### 5クエリの negative_hit_strength

| クエリ | 修正前 | 修正後 |
|---|---|---|
| MA_4115 | strong | none |
| HdrD1 (MA_0688) | medium | weak |
| MtpA (MA_4165) | weak | none |
| NifD (MA_3898) | medium | none |
| NifK (MA_3899) | strong | none |

### Tier A/B パートナー(calibration report のペア表の候補側)

調べた21タンパク質(5クエリ + Tier A/B のパートナー)のうち19が分類変化。
修正前に strong だった14件のうち11件(MA_4115, NifK, CdhC×2, MtpC, NifI1, NifI2,
MA_4548, MA_3992, MA_4546, MA_4550)は修正後 strength=none で **`Candidates` に入る
(Negative_unmatched)**。残る3件は Mer(MA_3733: none だが Candidates_relaxed +
Negative_hit)、MA_3998(strong→weak)、MA_1478(strong のまま)。変化なしは
MA_1478 と MA_4574(medium のまま)のみ。

### AlphaFold3 陰性セット(28件)

10件で strength/バケツが変化。strong 9→5、strict `Candidates` 収容は 9→11。
陰性セット全体の傾向は大きくは変わらないが、個別には入れ替わりがある
(例: MA_0165 strong→none、MA_0050 none→medium)。

## 影響: フルパイプライン再実行(スコア・ランキング)

修正前の Excel/Word は `data/output/before_negative_genome_fix/`(MA_4115_2.xlsx / .docx /
実行ログ。`data/output/` は git 管理外)に退避済み。以下は2通りで比較した。

### (1) リポジトリの config.yaml そのまま(5クエリ、候補ソースは Candidates / relaxed / No_hit)

修正前は退避Excel、修正後は 2026-09-20 01:33 の再実行(`data/output/MA_4115_2.xlsx`)。

- `02_Final_Score` の行数(クエリあたり): MA_4115 494→458、MA_0688 492→455、MA_4165 509→432、
  MA_3898 493→437、MA_3899 503→445(候補プールの組成が変わったため)。
- Tier2_Strong 件数: MA_4165 0→1、MA_3898 0→3、MA_3899 1→3、MA_0688 1→1、MA_4115 0→0。
- config.yaml のクエリ(MA_0688/4165/3898/3899)に対応する Tier A 既知ペア13行のうち、
  **修正前に最終スコア表へ載っていたのは2行(MtsF, NifD)のみ、修正後は11行**
  (載らない残りは DnaK / Hsp20。Tier B の12行はそのクエリ自体が config.yaml に無いため対象外)。
  修正前は negative_hit を有効にしていない限り Negative_hit 側に落ちて非表示になっており、
  保存性クエリ警告の設計書が想定した「真陽性が候補プールに入らない」問題そのものだった。
  取り違えが主因だったことになる。
- 副作用: `Candidates` が 127→709 件になり、`max_candidates_per_query: 200` の切り捨てが
  効くようになった。AF3 陰性28件のうちクエリ MA_4115 の表に載るのは 12→4 件。
- MA_4115 の Top10 は 1〜3位が同じ(MA_0363, MA_4110, MA_0361)、4位以降が入れ替わる
  (MA_4112/MA_1471 など relaxed+medium 由来が消え、Candidates 由来が上がる)。
- 保存性クエリ警告: この実行のログには**1件も出ない**(5クエリが全て weak/none になったため)。

### (2) キャリブレーション相当の設定で同一config・negativeフォルダのみ差し替え

10クエリ(MA_1111 は unresolved で除外)、全 candidate_sources 有効、
`max_candidates_per_query: 5000`、全バケツで GFF 注釈 ON。修正前は旧 Sulfolobus フォルダ
(git HEAD の内容)を別パスに再現して実行。Excel 書き出しは重すぎてハングしたため
(4,626行×10クエリ×全ソース)、最終スコア行のみをダンプして比較した。

| 指標 | 修正前 | 修正後 |
|---|---|---|
| 行数(クエリ×候補) | 46,260 | 46,260 |
| candidate_source = Candidates 行 | 1,270 | 7,083 |
| Tier A(n=13) interaction_score 平均 | 22.48 | 22.48 |
| Tier B(n=12) interaction_score 平均 | 29.67 | 29.67 |
| AF3 陰性(n=28) interaction_score 平均 | 10.66 | 10.66 |
| Tier A final_score 平均 | 30.09 | **36.51** |
| Tier B final_score 平均 | 35.77 | **44.10** |
| AF3 陰性 final_score 平均 | 25.14 | 25.14 |
| AUC(Tier A vs AF3陰性)interaction_score | 0.607 | 0.607 |
| AUC(Tier A vs AF3陰性)final_score | 0.563 | **0.684** |

- **`interaction_score` は全ペアで完全に不変**(negative genome に依存しないため。
  サニティチェックとして期待どおり)。したがって、キャリブレーション報告書の
  `interaction_score` に関する結論(Tier A と AF3 陰性の分離)は、この取り違えの影響を
  受けていない。
- **`final_score` は Negative_hit(strong等)だったペアで約 +8 点上がる**(例: MtpC 45.9→54.3、
  Tier3→Tier2)。final_score が negative_hit_strength を独立に読むため。既知の真陽性は
  上がり、AF3 陰性(元々 Negative_hit ではない)は動かないので、final_score の
  Tier A / AF3 陰性の分離は改善する(AUC 0.563→0.684。n が小さいので参考値)。
- クエリ別 Tier1-2 件数: MA_3997 1→6、MA_3998 1→6、MA_3992 4→6、MA_3898/MA_3899 1→3、
  MA_4165 0→1、MA_4548 4→4、MA_4115 0→0、MA_0688 1→1、MA_0658 0→0。
- 注意: 報告書に載っている絶対値(Tier A 平均 39.55 など)は今回の再実行では再現できない
  (報告書作成後にコードが変わっている: GEO の Mutual Rank 正規化、old_locus_tag
  フォールバックなど)。比較は同一コード・同一 config の前後のみで有効。

### 副次的な発見(未修正)

`annotation_targets.<bucket>.gff` が既定では candidates / candidates_relaxed だけ true のため、
Negative_hit バケツに入った候補は GFF 由来の `old_locus_tag` が付かず、STRING/GEO の
証拠が「not found」(MISSING)になる。実際に、GFF 注釈を全バケツで ON にしなかった
最初のキャリブレーション再実行では、修正前の Negative_hit 側の MtpC 等が STRING/GEO を
欠いていた(その結果を破棄して全バケツ ON で再実行し、上表はそちら)。
つまり Negative_hit バケツの真陽性は、BLAST 分類だけでなく証拠の取得面でも
不利になっている。
(2026-09-20 追記: 対応は別ブランチ `link-negative-hit-gff-annotation` で実装。
`_ANNOTATION_TARGETS_GFF_ALWAYS_ON` の2バケツ限定は性能上の制限ではなく、`beec3c8` で
`設定.xlsx` の既定プリセット(「考慮する」= Candidates/Candidates_relaxed のみ true)を
そのまま反映したものだった。GFF 注釈は 4,627 レコードでも 0.21 秒(既定の 1,261 件で 0.20 秒)。)キャリブレーション報告書の「STRING MISSING」の一部はこれが原因の可能性がある。

## 既存の分析結論への含意

- `claude/experimental_interactions_calibration_report.md` の**「Headline finding:
  既知の真陽性は Candidates に入らず Negative_hit に落ちる」は、この取り違えに
  大きく依存していた**可能性が高い。修正後は Tier A パートナーの大半が
  Candidates に入る。同レポートの rank 分母(2088/780/151)、
  「Negative_hit の解釈」に関する記述、および `config.yaml` の
  candidate_sources コメント(「保存性の高い経路のクエリでは negative_hit
  の有効化を検討」)の根拠の一部も同様に見直しが必要。
- `interaction_score` の Tier A vs AF3 陰性の比較(平均 39.55 vs 12.95)など
  スコア値そのものは、本調査では**再計算していない**。バケツの母集団と rank
  分母が変わるため、値の再確認が必要(再実行が必要。下記)。
- 一方で、保存性クエリ警告(`conserved_query_visibility`)機能は
  「保存性の高いクエリを持つ場合に警告する」という仕組み自体は有効だが、
  **修正後は上記5クエリのいずれも警告対象(strong/medium)にならない**。
  この機能の動機づけとなった実例はこの取り違えの産物だった。真に保存性の高い
  クエリ(例: 中核代謝酵素)で意味のある機能かどうかは、修正後データで
  改めて実例を探して確認するのが望ましい。

## 未実施・要判断

- フルパイプライン再実行(Excel/Word再生成、interaction_score/rank の再測定)は
  未実施。`output_excel` を上書きするため実行していない。
- calibration report 本体は書き換えていない(冒頭に注意書きのみ追記)。

## thermoplasma_acidophilum の .gitignore 除外の解消と起動時チェック(2026-09-20)

- 除外理由は記録になかったため、除外を解除して git 管理下に戻した。
  ローカルにあった `protein.faa` を NCBI から取得し直した `GCF_000195915.1` と比較し、
  **md5 が完全一致**(`ddebb59d...`)、`md5sum.txt` の検証もOK、ラベルとFASTA内容も
  一致していたため、ローカルのファイルをそのまま `git add`(ステージのみ、コミットは未実施)。
- `.gitignore` から該当行を削除。
- 起動時チェックを追加(`core/reference_genomes.py`、`main.py` の Configuration 節):
  - 期待される構成は `config/reference_genomes.v1.yaml`(フォルダ名 → RefSeq
    accession + FASTAが名乗るべき「Genus species」)に記載。
  - 検出するもの: 期待フォルダの欠落(マシン間差の検知)、リスト外フォルダ、
    accession不一致、FASTAの生物種不一致、**別フォルダのprotein.faaのbyte重複**、
    **同一種が陽性と陰性の両方に存在**(後者2つはmanifest無しでも動く)。
  - すべて WARNING でログに出すのみで実行は止めない(意図的な変更を妨げないため)。
    manifest が無い場合もその旨を WARNING に出す。チェック自体が例外を出しても実行は継続。
  - 各ゲノムの `category/label accession organism proteins=N md5=...` を INFO で
    毎回ログに残す。将来の再現性トラッキング(run provenance)で「どのゲノム構成で
    走ったか」の記録にそのまま使える。
  - テスト: `tests/test_reference_genomes.py`(各検出種別の単体テスト、tmp上で今回の
    取り違えを再現する結合テスト、リポジトリ実データ+実manifestの一致確認、
    manifest記載ゲノムが .gitignore されていないことの確認)。
- `run_provenance_design.md` はこのチェックアウトからは見つからなかった(別の場所/未同期?)
  ため、追記の確認・統合はしていない。上記ログ形式がそこの想定と合うか要確認。
