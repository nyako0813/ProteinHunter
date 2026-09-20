# Claude Code向け指示書: オペロン構造の直接予測(Rockhopper) — 実データ調査(Phase 6e候補、実装前の投資判断)

## 位置づけ・目的

これも`claude/phase6c_kegg_pathway_investigation.md`・
`claude/phase6d_interolog_investigation.md`と同じ**調査指示書**です。

## 背景・方針転換の経緯

当初は「オペロン予測データベース(DOOR等)」を調査対象として想定していた。
Cowork側の予備調査で以下が判明したため、方針を転換する。

1. **DOOR**(`csbl.bmb.uga.edu/DOOR/`)は接続タイムアウトで応答なし。
   **MicrobesOnline**もSSL証明書エラーで正常にアクセスできない。いずれも
   長期間メンテナンスされていない学術サイトの典型的な症状で、実質的に
   廃止されている可能性が高い(Cowork側からの単発確認のため、念のため
   ローカルからも一度だけ確認してほしい。下記「調査1」)。
2. 後継の**OperomeDB**も見つかったが、公式説明が"bacterial genomes"と
   明記されておりアーキア対応が不明瞭。
3. 一方、このプロジェクト自身の`claude/genomic_distance_weight_finding.md`
   (2026-09-04)で、「オペロン間隔の閾値が緩すぎる」という核心的な問題は
   **既に実際のGFF座標を使って直接解決済み**(150bp閾値、Mcr/Nif/Mtp複合体
   の実測値で検証済み)であることを再確認した。外部の汎用オペロン予測DBの
   出番は、少なくとも「隣接遺伝子ペアのオペロン判定」に関してはほぼ無い。
4. ただし同文書の「残課題」に明記されている通り、**同じオペロンの非隣接
   メンバー同士**(例: NifI1とNifD、間にNifI2・NifKを挟む)は単純な距離
   では捉えられず、現状はSTRING neighborhoodスコアに一部頼っている状態。
   ここは埋まっていない穴として残っている。

この穴を埋める候補として、外部データベースではなく**Rockhopper**という
ツールを見つけた。ゲノム配列+RNA-seqデータから直接オペロン構造(転写単位
全体、非隣接メンバーも含む)を予測する、無料で現在も公開されている
ソフトウェア(https://cs.wellesley.edu/~btjaden/Rockhopper/)。主にバクテリア
10種で検証された論文([PMC6776731](https://pmc.ncbi.nlm.nih.gov/articles/PMC6776731/))
だが、他グループによりアーキアにも適用された実績があるとの記述がある。

このプロジェクトは既にPhase 6b(coexpression)で*M. acetivorans*の
RNA-seqデータ(GSE77738, GSE64349)を扱っている。Rockhopperの入力として
再利用できる可能性が高く、**外部DBのカバレッジに頼らず、この生物種専用の
オペロンマップを直接作れる**という点で、STRING/KEGG/interolog推定より
筋が良い可能性がある。

## 調査してほしいこと

### 1. レガシーDBの生死確認(軽く、念のため)

```bash
curl -sI --max-time 15 http://csbl.bmb.uga.edu/DOOR/ ; echo "---"
curl -sI --max-time 15 https://www.microbesonline.org/ ; echo "---"
```

タイムアウト・接続エラー・証明書エラーになれば、Cowork側の確認(廃止と
推定)が裏付けられる。もし意外にもアクセスできた場合は、その中身
(*M. acetivorans*が収録されているか)を軽く確認してほしい。深追いは
不要(5分程度で切り上げてよい)。

### 2. Rockhopperの動作確認

```bash
java -version   # Rockhopperの実行にJavaが必要
```

Java未導入の場合はその旨を報告(導入コストとして記録)。

Rockhopperの配布物(JARファイル)をダウンロードし、簡単な動作確認を行う。
正式な使い方は公式サイト・マニュアルを参照
(https://cs.wellesley.edu/~btjaden/Rockhopper/)。

### 3. 入力データの準備状況の確認

Phase 6b(coexpression)でGSE77738/GSE64349のどのファイルを既にダウンロード
済みか確認してほしい。Rockhopperはリード配列(FASTQ、可能ならペアエンド)を
必要とする可能性が高く、Phase 6bで使ったのがGEOの処理済みサマリファイル
(XLSX)だけだった場合、SRA(GSE77738の生リードはSRP069835に登録済み)から
改めてFASTQを取得する必要がある(`prefetch` + `fasterq-dump`、SRA toolkit)。

```bash
# 既存のダウンロード済みファイルを確認
find data -iname "*GSE77738*" -o -iname "*GSE64349*" -o -iname "*SRP069835*"
```

生リードが無ければ、GSE77738のサンプルのうち**1〜2件だけ**を試験的に
取得する(全61サンプルを取得する前に、まず動くかどうかの確認が目的)。

### 4. 試験実行と既知オペロンでの検証

ゲノム配列(`data/input/genome.gff`に対応するFASTA。既存の候補生成
パイプラインが使っているものを流用できるはず)と、取得したRNA-seqリードを
Rockhopperに投入し、オペロン予測を実行する。

実行できたら、`claude/genomic_distance_weight_finding.md`で既に実測・
検証済みの以下の既知オペロンが、Rockhopperの予測でも正しく1つの転写単位
として出力されるか確認する(これが今回の投資判断の核心):

| 既知オペロン | 遺伝子 | 備考 |
|---|---|---|
| Mcr活性化複合体クラスタ | MA_4546, MA_4547, MA_4548, MA_4549, MA_4550 | 隣接gap 2〜25bp |
| Nif複合体 | MA_3896(NifI1), MA_3897(NifI2), MA_3898(NifK), MA_3899(NifD) | 隣接gap 3〜28bp。**NifI1-NifDは非隣接**、これがまさに埋めたい穴 |
| Mtp複合体 | MA_4165(MtpA), MA_4164(MtpC) | gap 70bp |

特に**NifI1(MA_3896)とNifD(MA_3899)が同じオペロンとして正しくグループ化
されるか**を重点的に確認すること。ここがRockhopper導入の主な価値になる。

### 5. コスト評価

- 実行時間(1サンプルあたり、全ゲノムスケールでどの程度かかるか)
- 計算資源(メモリ、ディスク)
- 全61サンプル(GSE77738)を処理する場合の概算時間
- STRING(Phase 6a)と同様、「一度だけ計算してキャッシュする」運用が
  現実的かどうか

## 成果物として欲しいもの

これまでと同じ体裁のレポート。以下を含むこと:

1. DOOR/MicrobesOnliveの生死確認結果
2. Rockhopperの動作確認結果(Java環境、ダウンロード、実行可否)
3. 入力データ準備状況(Phase 6bの既存データで足りるか、SRAから追加取得が
   必要か)
4. 既知オペロン(特にNifI1-NifD)での検証結果
5. コスト評価
6. 「実装する価値があるか」についての率直な所見。特に、既に
   `_gene_neighborhood_v2`で解決済みの「隣接ペアのオペロン判定」とは別の
   価値(非隣接メンバーの検出)がどれだけ実際に得られそうか。

## 非対応・スコープ外

- 実際のコード実装(`analysis/operon_bridge.py`相当、`genomic_context`
  カテゴリへの新規コンポーネント追加等)は、この調査が完了しCowork側で
  結果を確認してから別途指示する。
- 全61サンプルの本番実行は今回の調査範囲外(1〜2サンプルでの動作確認・
  精度検証までで十分)。
