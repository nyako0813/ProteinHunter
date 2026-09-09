# GEO共発現スコアのMutual Rank正規化(coexpression_gse77738/gse64349)

2026-09-05、MA_4115クエリでの実データ確認を発端に実装した記録。

## 背景・きっかけ

`consider_cross_species_matches`スイッチ対応の一環で、クエリ(MA_4115)の
`resolved_old_locus_tag`が空になっていた不具合(target FASTAにold_locus_tag
情報が存在しないデータセットで、protein_id解決時にold_locus_tagが引き継がれ
ない問題)を修正し、`query_proteins[].old_locus_tag`をフォールバックとして
使うようにした(別コミット)。この修正によりSTRING/GEO共発現の両方が初めて
MA_4115クエリに対して実際に機能するようになったが、その結果順位が大きく
入れ替わり、ユーザーから「STRINGは疎(正しく機能)なのに対し、GEO共発現の
値が広範囲に高すぎるのではないか」という指摘があった。

## 実データでの確認

MA_4115クエリ、351候補(interaction_scoring評価対象)で集計:

| 指標 | 中央値 | 平均 | 0.7超の割合 |
|---|---:|---:|---:|
| `coexpression_gse77738`(既存の一方向percentile) | 0.663 | 0.621 | 47.7% |
| `coexpression_gse64349`(既存の一方向percentile) | 0.648 | 0.605 | 42.4% |
| `string_cooccurrence`(参考、STRING pscore) | 0.0 | 0.073 | 0.0% |

無シグナルなら中央値0.5付近に分布するはずが、GSE77738/GSE64349はいずれも
半数近くが0.7超という明らかな底上げが見られた。一方STRING側は健全(疎)に
機能していた。さらに、GSE77738とGSE64349は独立したデータセットにも
関わらず、両方に値を持つ291候補間でPearson r=0.549(両方0.7超が33.7%、
両方0.3未満が5.8%)という強い相関があり、「クエリ固有の共発現シグナル」
というよりも「多くの候補に共通する何らかの性質(coexpressionハブ、または
測定精度由来のバイアス)」が両データセットに共通して乗っていることを示唆
していた。

これはPhase 6b(`claude/phase6b_coexpression_design.md`)で既に指摘されていた
「GSE64349の小標本による背景相関インフレ」および「GSE77738のAlphaFold3
陰性セットでの逆転」と同根だが、既存の一方向(クエリ側のみ)percentile-rank
正規化だけでは解消しきれていなかったことが、MA_4115の実データで改めて
定量的に確認された形。

## ユーザーの要望

「除外とまではいかないまでも、スコアを正規化する必要はあるように思う」との
方針。既存のweight調整(`coexpression_gse64349`の1/3重み)や`interaction_score`
からの`coexpression_gse77738`除外(PR #11)とは別に、スコアそのものの計算式を
見直す方向で対応した。

## 原因分析

既存の実装(`analysis/coexpression_bridge.py::_build_bundle_from_matrix`)は、
相関値をクエリ遺伝子自身の背景分布(クエリ対他の全遺伝子)に対する
percentile rankに変換していた。これは「クエリから見てこの候補は他の遺伝子
より際立って相関が高いか」だけを見ており、「候補から見てクエリは他の遺伝子
より際立って相関が高いか」を見ていない。そのため、単に多くの遺伝子と広く
相関する遺伝子(coexpressionの「ハブ」。実際には成長期に伴う発現量変化を
強く受ける遺伝子や、測定精度が高くノイズが少ない高発現遺伝子である場合が
多い)が、クエリとの関係とは無関係に高スコアを得てしまう。

## 対応: Mutual Rank正規化

`normalized_value`を、クエリ側percentileと候補側percentile(候補遺伝子自身の
背景分布に対する、同じ相関値のpercentile)の幾何平均に変更した。
CoNekT/ATTED-IIなど共発現ネットワーク解析ツールが標準的に採用している
"Mutual Rank"の考え方を、既にpercentile化された値に適用したもの。

- 真にクエリ特異的なパートナー(候補自身は他の遺伝子とあまり相関しない)は、
  候補側から見てもクエリとの相関が際立つため、両側percentileが高く保たれる。
- ハブ遺伝子(候補が多くの遺伝子と広く相関する)は、候補自身の背景分布も
  底上げされているため、候補側percentileは相対的に下がり、結果として
  mutual rankも下がる。

除外ではなく正規化なので、GSE77738/GSE64349とも計算対象からは外していない
(Phase 6bで`interaction_score`から除外されたGSE77738の扱いは無変更)。

### 実装

- `analysis/coexpression_bridge.py`:
  - 内部専用の`_OneSidedCoexpression`(旧`CoexpressionPairValue`相当、
    ディスクキャッシュの保存形式は無変更)と、公開用の新しい
    `CoexpressionPairValue`(`query_percentile`/`candidate_percentile`/
    `percentile`(mutual)の3値を保持)に分離。
  - `CoexpressionBundle`に相関行列(`_z_scores`/`_gene_index`/`_gene_order`)を
    保持するフィールドを追加し、`lookup()`が候補側の背景分布をオンデマンドで
    計算できるようにした。候補ごとの背景分布はプロセス内メモリのみで
    メモ化し(`_background_cache`)、ディスクキャッシュ(JsonCache)には
    書き込まない。理由: 候補は1クエリあたり数百件になり得るため、
    クエリ用ディスクキャッシュと同じ方式で全候補分の背景分布(遺伝子ごと
    ~数千エントリ)を永続化すると、キャッシュファイルが数百倍に膨れ上がる。
    行列演算自体はマイクロ秒オーダーで安価なため、プロセスごとの再計算で
    十分と判断した。
  - 既存のディスクキャッシュ形式(`{correlation, percentile}`)は変更して
    いないため、この修正適用後もユーザーの既存キャッシュファイルは
    そのまま再利用される(再ダウンロード不要)。
- `analysis/interaction_scoring.py::_coexpression_status_and_value`:
  説明文(`explanation`)にmutual/query_side/candidate_sideの3値を全て表示する
  ように変更(監査性のため)。

### テスト

`tests/test_coexpression_bridge.py`に、クエリ特異的パートナーとハブ遺伝子を
明示的に作り分けた合成データセットで検証するテストを追加
(`test_gse77738_mutual_rank_deflates_hub_genes`)。ハブは一方向percentileの
時点で既にパートナーより低いが(0.6 vs 1.0)、mutual rank適用後はさらに差が
開くこと(候補側percentile 0.4 vs 1.0、mutual 0.49 vs 1.0)を確認済み。
既存テストは全てpassすることを確認(フルスイート503 passed、既存の2件の
失敗はGSE77738キャッシュ再利用に関する別の既知の不具合で本変更とは無関係、
下記参照)。

## ついでに発見した別の不具合(追記: 2026-09-05に対応済み)

`load_gse77738_coexpression_bundle`/`load_gse64349_coexpression_bundle`は、
ディスクキャッシュにクエリの結果が既に存在する場合でも、ソースファイル
(`GSE77738_ReadCounts.xls`等)の存在チェック・読み込みを無条件に先に実行して
しまう。ソースファイルが削除/未取得の状態だとオフライン環境で再ダウンロード
を試みて失敗し、キャッシュが存在するにも関わらず空バンドルにフォールバック
してしまう。`tests/test_coexpression_bridge.py::test_result_is_cached_and_reused_without_reparsing`
と`tests/test_interaction_scoring.py::test_coexpression_gse77738_reuses_cache_across_runs`
が既にこれを検出していた(本コミット適用前のmainでも同じく失敗することを
確認済み)。当初はMutual Rank対応とは別の独立した不具合として修正を見送って
いたが、ユーザーからの追加依頼(「キャッシュ関連の調査もお願いします」)を
受けて調査・修正した。

### 原因(補足)

キャッシュ判定自体は`_build_bundle_from_matrix`のクエリごとのループ内で
行われていたが、そのループに到達する前に、両ローダー関数がソースファイルの
`.exists()`チェック→(無ければ)ダウンロード→`pd.read_excel(...)`を無条件に
実行していたため、キャッシュが完全に揃っているケースでもソースファイルへの
依存が残っていた。

さらに、今回のMutual Rank対応で候補側percentileの計算に相関行列
(`_z_scores`/`_gene_index`)そのものが必要になった点も考慮が必要だった。
この行列はクエリ単位ではなくデータセット単位で1つだけの、比較的小さい
オブジェクト(GSE77738で遺伝子数×サンプル数 ≒ 4540×13 ≒ 数万要素、JSON化
してもおよそ1-2MB程度)であり、候補ごとの背景分布(`_background_cache`、
候補は1クエリあたり数百件になり得るため意図的にプロセス内メモリのみに
留めている)とは性質が異なる。そのため、この行列だけは既存の
クエリ単位キャッシュとは別の専用ネームスペース
(`{cache_namespace}_matrix`)に1エントリとして永続化することにした。

### 対応

- `analysis/coexpression_bridge.py`に以下を追加:
  - `_cache_matrix` / `_load_cached_matrix`: データセット単位の相関行列
    (`gene_order`/`z_scores`/`n_samples`)を専用ネームスペースに保存・復元。
  - `_try_fully_cached_bundle`: 行列キャッシュと、設定された全クエリタグ分の
    ディスクキャッシュが両方揃っている場合に限り、ソースファイルに一切
    触れずに完全な`CoexpressionBundle`(Mutual Rankに必要な相関行列込み)を
    再構築する。クエリのうち1つでも未キャッシュならNoneを返し、
    呼び出し側は従来通りファイル読み込みにフォールバックする(中途半端な
    バンドルを返さないための安全策)。
  - 両ローダー関数の冒頭(`enabled`チェック直後)で`_try_fully_cached_bundle`
    を呼び出し、成功すればソースファイルの存在チェック/ダウンロード/
    読み込みを完全にスキップする。
  - `_build_bundle_from_matrix`内で相関行列を計算した直後に`_cache_matrix`
    を呼び出し、ファイル読み込みが発生した回は必ず行列キャッシュを更新する
    (既存のクエリ単位キャッシュの保存形式・タイミングは無変更)。
- 既存のディスクキャッシュ形式(クエリ単位の`{correlation, percentile}`)は
  変更していないため、この修正適用後もユーザーの既存キャッシュファイルは
  そのまま再利用される。

### テスト

`test_result_is_cached_and_reused_without_reparsing`と
`test_coexpression_gse77738_reuses_cache_across_runs`が、追加コードなしで
そのままpassするようになったことを確認(いずれもソースファイルを
`unlink()`した後にキャッシュのみで再実行するテスト)。フルスイート
505 passed(既知の失敗は解消、新規の失敗なし)。
