# Claude Code向け指示書: KEGGパスウェイ共起の実データ調査(Phase 6c候補、実装前の投資判断)

## 位置づけ・目的

これは**実装指示書ではありません**。`claude/phase6_external_evidence_design.md`
(Phase 6a: STRING PPI evidence)と同じやり方で、コードを書く前に実データで
価値を検証するための**調査指示書**です。

Cowork側でKEGG REST APIをWebFetch経由(要約モデル通過)で予備調査したところ、
以下の概算値が得られました(いずれも精度の低い概算です):

- ゲノム全体(約4,500遺伝子)のうちKO(KEGG Orthology)が付与されているのは
  約59%
- PATHWAYが付与されているのは約15%
- MA_4115自身にはKO(`K07585`, tRNA methyltransferase)は付与されているが、
  具体的なPATHWAYへの割り当ては無し("Unclassified: genetic information
  processing"止まり)
- MA_0361(electron_carrier候補、フラボドキシン)はKO・PATHWAYともに無し

Phase 6aのSTRING調査と似た構図(生物種全体では一定のカバー率があるが、
個々のクエリでは信号が乏しい可能性がある)に見えるため、正確な数字を
ローカル環境から直接KEGG REST APIを叩いて確認し、実装する価値があるか
判断したい。

## 生物種・識別子(確認済み)

- KEGG生物種コード: `mac`(Methanosarcina acetivorans C2A、T00080)
- 遺伝子ID形式: `mac:MA_XXXX`(アンダースコア付き旧ロカスタグ。STRINGの
  `188937.MA_XXXX`と同じ`old_locus_tag`をそのままキーに使える。新しい
  識別子変換ロジックは不要)

## 調査で確認したいこと

### 1. カバー率の正確な集計

STRING調査(Phase 6a)と同様、**個別の遺伝子ごとにAPIを叩くのではなく、
一括ダウンロード型のエンドポイントを使う**こと(KEGGの利用規約上も推奨されて
いる方式。過度な個別リクエストは避ける)。

```python
import requests
import time

BASE = "https://rest.kegg.jp"

def fetch(path: str) -> str:
    resp = requests.get(f"{BASE}/{path}", timeout=30)
    resp.raise_for_status()
    time.sleep(1)  # KEGGの利用ガイドラインに配慮し、リクエスト間隔を空ける
    return resp.text

# 生物種全体の遺伝子一覧(総遺伝子数の確認。既存のUniProtバルクデータの
# 4492件と比較して整合性も確認する)
genes = fetch("list/mac").strip().splitlines()
print(f"total genes (KEGG list/mac): {len(genes)}")

# 遺伝子→KOのリンク(タブ区切り: mac:MA_XXXX<TAB>ko:KXXXXX)
ko_links = fetch("link/ko/mac").strip().splitlines()
genes_with_ko = {line.split("\t")[0] for line in ko_links if line}
print(f"genes with KO: {len(genes_with_ko)} ({len(genes_with_ko)/len(genes)*100:.1f}%)")

# 遺伝子→PATHWAYのリンク(タブ区切り: mac:MA_XXXX<TAB>path:macXXXXX)
pathway_links = fetch("link/pathway/mac").strip().splitlines()
genes_with_pathway = {line.split("\t")[0] for line in pathway_links if line}
print(f"genes with pathway: {len(genes_with_pathway)} ({len(genes_with_pathway)/len(genes)*100:.1f}%)")

# パスウェイごとの遺伝子数分布(大きすぎるパスウェイは特異性が低く、
# 「同じパスウェイ」というだけで無関係な遺伝子ペアを大量に拾ってしまう
# リスクがある。domain_family_mapのfold-level過剰一致と同種の懸念)
from collections import Counter
pathway_sizes = Counter(line.split("\t")[1] for line in pathway_links if line)
print("largest pathways by gene count:")
for pw, count in pathway_sizes.most_common(10):
    print(f"  {pw}: {count} genes")
```

### 2. 現在の主要クエリでの実際のシグナル

以下4クエリについて、それぞれのPATHWAY割り当てと、既存のCandidates/
Candidates_relaxedバケツの候補との共起(同じPATHWAYに属する候補が何件
あるか)を確認する。

- MA_4115(`mac:MA_4115`)
- MA_0795(`mac:MA_0795`)
- MA_0826(`mac:MA_0826`)
- MA_0074(`mac:MA_0074`)
- MA_0266(`mac:MA_0266`)

```python
def query_pathways(old_locus_tag: str) -> set[str]:
    gene_key = f"mac:{old_locus_tag}"
    return {line.split("\t")[1] for line in pathway_links if line.startswith(gene_key + "\t")}

for tag in ["MA_4115", "MA_0795", "MA_0826", "MA_0074", "MA_0266"]:
    pws = query_pathways(tag)
    print(f"{tag}: pathways = {pws or 'none'}")
```

各クエリについて、`pws`が空でなければ、そのクエリの現在の
Candidates/Candidates_relaxedバケツ(既存のExcel出力や候補生成ロジックから
取得)の中に同じPATHWAYを持つ候補が何件あるかも確認すること。0件でも
構わない(STRING同様、「今回のクエリでは信号が乏しいが将来のクエリの
ための汎用機能」という結論もあり得る)。

### 3. 利用規約・レート制限の確認

STRING(Phase 6a)と同様、KEGGの利用規約ページ(`https://www.kegg.jp/kegg/rest/`
または`https://www.kegg.jp/kegg/legal.html`)を確認し、以下を明記すること:

- REST APIの学術利用における制限(無料か、商用利用の条件、アカデミック
  ライセンスとの関係)
- 推奨されるリクエスト間隔・one-time bulk downloadの可否(STRINGのような
  まとめてダウンロード可能なファイル配布があるか、REST APIを都度叩く
  しかないか)
- 引用・帰属表示の要否

### 4. キャッシュ方針の検討

STRING(Phase 6a decisions 6)と同様の設計にできるか検討する:
`link/ko/mac`・`link/pathway/mac`の結果をまるごとローカルファイルとして
一度だけ取得・保存し(生物種全体で数千行程度、STRINGの76MBに比べれば
遥かに小さいはず)、以後はそのローカルファイルを参照する形にできないか。

## 成果物として欲しいもの

`claude/phase6_external_evidence_design.md`の「Investigation: STRING」
セクションと同じ体裁で構いません。以下を含むレポート(Markdown、コード
不要、調査結果のみ):

1. カバー率の正確な数字(KO・PATHWAY、それぞれ件数とパーセンテージ)
2. 上位パスウェイの遺伝子数分布(特異性が低すぎるパスウェイが無いか)
3. 5クエリ(MA_4115, MA_0795, MA_0826, MA_0074, MA_0266)それぞれの
   PATHWAY割り当てと、候補バケツ内での共起件数
4. 利用規約・レート制限のまとめ
5. 「実装する価値があるか」についての率直な所見(実装コストに見合う
   シグナルの強さがあるか、STRINGのケースとの比較)

## 非対応・スコープ外

- この調査結果を元にした実際のコード実装(`analysis/kegg_pathway_bridge.py`
  相当のもの)は、この調査が完了し、Cowork側で調査結果を確認してから
  別途指示する。
- パスウェイの意味的な近さ(例えば同じ上位カテゴリの姉妹パスウェイ)を
  考慮した高度な類似度計算は今回の調査範囲外。単純な「同じPATHWAY IDを
  共有しているか」の有無のみを見る。
