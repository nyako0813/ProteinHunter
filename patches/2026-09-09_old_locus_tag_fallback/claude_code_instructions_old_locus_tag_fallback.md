# Claude Code 実装指示書: resolved_old_locus_tag フォールバック未適用問題の修正

## 背景・目的

`domain_complementarity` v3 (UniProt Pfam/InterProベースのドメインファミリー照合、コミット `81b4c44`) を実データ (MA_4115 クエリ) で検証したところ、`06_Functional_Domain_Evidence` / `08_Genomic_Context` 両シートで `query_old_locus_tag` が `None` のままで、MA_0361 / MA_0363 の `domain_complementarity_score` が期待通り 0→正の値 に変化しませんでした (`explanation` も "generic-only description overlap ignored" のまま)。

調査の結果、原因は今回実装した機能のバグではなく、**もっと手前の段階 — クエリ解決 (`_resolve_query`) が `resolved_old_locus_tag` を空文字列のまま返している** ことでした。`uniprot_domain_map.get(query["resolved_old_locus_tag"])` は空文字列では絶対にヒットしないため、新しいドメインファミリー照合層は静かにスキップされ、既存の自由文一致ロジックにフォールスルーしていました。

この `resolved_old_locus_tag` が空になる問題自体は **以前のセッションで一度修正・テスト済みでした** が、実装したブランチが GitHub の `origin/main` に一度も push されていなかったため、現在のリポジトリには存在していません (`git log --all --grep` で見つからなかったのはこれが理由です)。今回、その未pushの修正一式を発掘し、現在の `origin/main` (`81b4c44`, domain_complementarity v3 適用後) に対してクリーンに適用できることを確認しました。添付のパッチをこのまま適用してください。

## 添付ファイル

- `0001-old_locus_tag_fallback.patch` — **必須**。本題の修正。
- `0002-geo_mutual_rank_normalization.patch` — 任意 (推奨)。詳細は下記「付随する2つの改善」参照。
- `0003-geo_coexpression_cache.patch` — 任意 (推奨)。同上。

3つとも `origin/main` (`81b4c44`) にクリーンに適用できることを確認済みです。適用後、3パッチ込みで `pytest` 全件 (523件) がパスすることも確認済みです。

## 実装タスク

### 1. パッチの適用

```bash
cd ~/projects/ProteinHunter
git status   # 作業ツリーがクリーンであることを確認してから実行してください

# 必須
git apply --3way 0001-old_locus_tag_fallback.patch

# 推奨 (任意。適用する場合は0001の後に)
git apply --3way 0002-geo_mutual_rank_normalization.patch
git apply --3way 0003-geo_coexpression_cache.patch
```

`git am` ではなく `git apply --3way` を使ってください — こちらの検証環境では `git am` が index blob 不一致エラーで失敗しましたが (autocrlf 関連の環境差異と思われます)、`git apply --3way` は3パッチとも "cleanly" で適用できました。万一 `git apply --3way` でも競合する場合は、パッチ内容を参考に手動で反映してください (変更点は小さく、`analysis/interaction_scoring.py` の `_resolve_query` 関数が中心です)。

### 2. 0001 の変更内容 (必須修正の要点)

`analysis/interaction_scoring.py` の `_resolve_query` 関数:

- **Before**: ターゲットレコードが見つかった場合、`resolved_old_locus_tag = matched_record.old_locus_tag or ""` — レコード自体に old_locus_tag が無ければ常に空文字列。
- **After**: レコードに old_locus_tag が無い場合、`query["old_locus_tag"]` (= config の `query_proteins[].old_locus_tag`) をフォールバックとして使用。レコード自体の old_locus_tag が存在する場合は常にそちらが優先されます (config値がレコードの実データを上書きすることはありません)。両方とも空なら従来通り空文字列のまま。

`config.yaml` にはこの挙動を説明するコメントが追加されます (`query_proteins.old_locus_tag` は任意入力だが、target FASTA が old_locus_tag を持たない生物種では STRING 連携のために必須、という内容)。

`tests/test_interaction_scoring.py` に単体テストが3件追加されます (フォールバックが効くケース/レコード優先のケース/両方空のケース)。

### 3. 【重要・パッチとは別に必要な作業】config.yaml の `old_locus_tag` 設定

このパッチは「config に `old_locus_tag` が設定されていればフォールバックとして使う」というロジックを実装するだけです。**MA_4115 クエリの実行に使っている実際の config.yaml (パッチが当たる、リポジトリ管理下のテンプレートではなく、実行時に読み込んでいる方) の `query_proteins` エントリに `old_locus_tag: "MA_4115"` を設定しないと、修正の効果が出ません。**

```yaml
interaction_scoring:
  query_proteins:
    - protein_id: "WP_011024006.1"   # 既存の値のまま
      old_locus_tag: "MA_4115"        # ← これを追加/設定してください
      sequence: ""
```

パッチ適用後、この設定が漏れていないか必ず確認してください。

### 4. 動作確認

```bash
python -m pytest -q   # 523 (パッチ全部適用時) または該当件数がパスすることを確認
```

その後、MA_4115 クエリで実際にパイプラインを再実行し、以下を確認してください:

1. 出力 Excel の `06_Functional_Domain_Evidence` シートおよび `08_Genomic_Context` シートで、`query_old_locus_tag` (or 相当のクエリ側 old_locus_tag 列) が `None` ではなく `MA_4115` になっていること。
2. `06_Functional_Domain_Evidence` シートで MA_0361 / MA_0363 の `domain_complementarity_score` が `0` から正の値に変化し、`explanation` に "domain family match (UniProt Pfam/InterPro): ..." が含まれること (domain_complementarity v3, コミット `81b4c44` の効果がここで初めて実データに反映されるはずです)。
3. (0002/0003 も適用した場合) STRING/GEO 系のカラム (`string_cooccurrence` / `string_neighborhood` / `coexpression_gse77738` / `coexpression_gse64349` 等) が MA_4115 に対して `MISSING` ではなく実際の値を返すようになっていること。

## 付随する2つの改善 (0002 / 0003, 任意)

調査中に、`old_locus_tag` フォールバック修正と同じ未push作業ブランチ上に、それを土台にした2つの追加改善が既に実装・テスト済みの状態で見つかりました。今回の domain_complementarity 検証には必須ではありませんが、内容だけ共有します。適用するかどうかはユーザー判断でお願いします。

### 0002: GEO coexpression スコアの Mutual Rank 正規化

`old_locus_tag` フォールバック修正後、実際に MA_4115 で GEO coexpression (`coexpression_gse77738` / `coexpression_gse64349`) の評価が届くようになったところ、候補の約半数が 0.7 パーセンタイル超という偏った分布になっており (正しく判別できているなら 0.5 付近が中心のはず)、かつ独立した2つのデータセット間でパーセンタイルの相関が r=0.549 と高すぎる (クエリ固有のシグナルというより何らかの系統的バイアス) ことが分かりました。STRING側 (`string_cooccurrence`/`string_neighborhood`) は同じ候補セットで中央値0という健全な疎さを保っていたため、対比としてこれが問題として浮上した形です。

原因は `normalized_value` がクエリ遺伝子側のパーセンタイルランクのみで計算されており、「候補側から見てもクエリが外れ値的に高相関か」という対称な問いを見ていなかったことでした。Mutual Rank (両方向のランクの幾何平均) に切り替えることで解消しています。

### 0003: GEO coexpression 相関行列のキャッシュ

`load_gse77738_coexpression_bundle` / `load_gse64349_coexpression_bundle` が、ディスクキャッシュを確認する前に無条件で元データソースファイルの存在チェック/ダウンロード/パースを行っていたため、キャッシュが揃っているクエリでもソースファイルが (例: 初回使用後に削除されているなどで) 見つからないと失敗、または空のバンドルに黙って劣化していました。データセット全体の相関行列を1回だけキャッシュし、要求された全クエリタグが既にキャッシュ済みならソースファイルを一切読まずに済ませる `_try_fully_cached_bundle` を追加し、それ以外は既存の読み込みパスにフォールバックする形です。

## 非破壊要件

- `legacy_additive` スコアリングモデルの挙動には一切触れません (今回の変更は `v2_evidence_based` 経路のクエリ解決部分のみ)。
- `resolved_old_locus_tag` は「ターゲットレコード自体の値」が常に「config 側の値」より優先されます。config 側の値がレコードの実データを上書きすることはありません。
- config で `old_locus_tag` を設定していないクエリの挙動は完全に従来通りです (フォールバック元が空なら結果も空)。

## 対象外

- `sulfur_carrier` / `radical_sam` の Pfam ID 拡充 (`config/domain_family_map.v1.yaml`) — 別途対応予定、今回は対象外。
- STRING/GEO 以外のエビデンスソースへの `resolved_old_locus_tag` の影響範囲の網羅的な再検証。

## 参考

- `claude/domain_complementarity_v3_design.md` (設計ドキュメント本体)
- `claude/coexpression_mutual_rank_normalization.md` (0002 の設計メモ。0002適用でこのファイルがリポジトリに追加されます)
