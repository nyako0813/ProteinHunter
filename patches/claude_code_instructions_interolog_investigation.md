# Claude Code向け指示書: BioGRID/IntActのinterolog推定 — 実データ調査(Phase 6d候補、実装前の投資判断)

## 位置づけ・目的

これも`claude/phase6c_kegg_pathway_investigation.md`と同じ**調査指示書**です。
実装前に、実データで価値とコストを検証する。

Cowork側で公開ドキュメント・APIをWebFetch経由で予備調査した結果は以下の通り。
一部はレート制限(HTTP 429)でリアルタイム確認できなかった箇所があり、
ローカル環境からの直接確認が必要。

## interolog推定とは(前提の確認)

*M. acetivorans*自体には直接の実験的PPIデータがほとんど無い
(`claude/phase6_external_evidence_design.md`のSTRING調査で確認済み: 全
チャンネルで信号が乏しい)。interolog推定は、よく研究された他生物種
(大腸菌・酵母・ヒト、あるいはよく研究されたアーキア)で**既に実験的に
確認されている**物理的相互作用ペアを、オルソログ関係を介して
*M. acetivorans*の遺伝子ペアに"転移"させる手法。

このプロジェクトの`claude/domain_family_map_v2_sulfur_relay_expansion.md`
に既に整理されている「他生物種で構造的に確認されている相互作用ドメイン
ペア」(ThiS×ThiF、IscS×IscU、TusA×TusBCD×TusE、SufS×SufE等、いずれも
主に大腸菌での構造既知の相互作用)は、まさにinterolog推定の元ネタになる
候補群。

**重要な注意**: interolog推定は「オルソログが相互作用する」ことからの
**推論**であり、*M. acetivorans*自身での直接的な証拠ではない。domain_
complementarityの「ドメインの組み合わせとしては妥当」という限界
(`claude/alphafold3_candidates_MA_0826.md`のMA_0826×MA_0361/MA_0363の件)
と同種の弱さを持つ点を、スコアリング設計上も明示すること(過信させない)。

## Cowork側で確認できたこと(公開情報、要検証)

### IntAct(EBI)

- 旧IntActウェブサイトは2025年12月15日に閉鎖され、新しい「IntAct Portal」
  に移行済み([参考](https://www.ebi.ac.uk/legacy-intact/))。
- REST APIは4系統([ドキュメント](https://intact-portal.gitbook.io/intact-portal-doc/)):
  - Interactor検索: `https://www.ebi.ac.uk/intact/ws/interactor/` (Swagger UI: `.../swagger-ui.html`)
  - Interaction検索: `https://www.ebi.ac.uk/intact/ws/interaction/`
  - Network Viewer: `https://www.ebi.ac.uk/intact/ws/network/`
  - Data/graph API: `https://www.ebi.ac.uk/intact/ws/graph/`
- APIキー不要という記述は見当たらず(BioGRIDと違い登録不要の可能性が高いが、未確認)。
- バルクダウンロードあり。MITAB/PSI-MI形式、"Species-based Datasets"という
  生物種別のダウンロードカテゴリが存在することは確認できたが、正確な
  URL・ファイル一覧はCowork側では取得できなかった(該当ページが404、
  もしくはJavaScriptレンダリングでWebFetchから内容を取得できず)。
- `findInteractions/iscS?taxonId=511145`(大腸菌K12 MG1655、IscSで検索)を
  実際に試みたが、EBI側のレート制限(HTTP 429)で複数回失敗し、実データでの
  動作確認ができていない。**ローカル環境で改めてこのクエリを試すこと**。

### BioGRID

- REST APIは**無料のAPIキー登録が必須**([登録フォーム](https://webservice.thebiogrid.org))。
  Cowork側では登録できないため、**ユーザー本人による登録が必要**。
- クエリ例: `/interactions/?geneList=MDM2&taxId=9606&searchNames=true`
  のような形式(`geneList`はパイプ区切りで複数指定可、`taxId`はNCBI
  分類ID)。1リクエストあたり`max`パラメータの上限は10,000件、それ以上は
  ページネーションが必要。月次更新(毎月4日)。
- データライセンス: **MITライセンス**(帰属表示以外の制約がほぼ無く、
  商用・自動化パイプラインへの組み込みも許容範囲内)。STRING(CC BY 4.0)・
  KEGG(学術利用無料だがレート制限厳格)より扱いやすい可能性がある。
- バルクダウンロード: `https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-5.0.261/`
  に`BIOGRID-ORGANISM-5.0.261.mitab.zip`という、複数生物種のMITABファイルを
  まとめたzipがある。**中身(生物種ごとのファイル一覧、アーキアが含まれる
  かどうか)はzipを展開しないと確認できず、Cowork側からは未確認**。

## 調査してほしいこと

### 1. 既存のBLASTオルソログ基盤の確認(最優先)

このパイプラインは既に`positive_dir`配下の参照ゲノムに対してBLASTを
行っている(`config.yaml`の`paths.positive_dir`)。

```python
# data/databases/positive/ 配下にどの生物種のゲノムが置かれているか確認
import os
for entry in sorted(os.listdir("data/databases/positive")):
    print(entry)
```

もし大腸菌(*E. coli*)や酵母(*S. cerevisiae*)が既にこの中に含まれて
いれば、interolog推定に必要な「M. acetivorans遺伝子 → 大腸菌/酵母の
オルソログ」というマッピングは、**新たなBLAST実行なしで既存の候補生成
ロジックの副産物としてほぼ無料で手に入る**可能性が高い。含まれていなければ、
新たにBLAST参照ゲノムを追加する作業が必要になり、実装コストが一段上がる。

### 2. IntAct API疎通確認(ローカルから)

Cowork側でレート制限にかかった以下のクエリを、ローカル環境から再実行して
結果を確認してほしい。

```python
import requests
import time

# 大腸菌(taxonId=511145, K12 MG1655)でIscSを検索
resp = requests.get(
    "https://www.ebi.ac.uk/intact/ws/interaction/findInteractions/iscS",
    params={"taxonId": 511145},
    timeout=30,
)
print(resp.status_code)
print(resp.text[:2000])
time.sleep(1)

# ThiS-ThiF(claude/domain_family_map_v2_sulfur_relay_expansion.mdでPDB 1ZUDとして既出)も試す
resp2 = requests.get(
    "https://www.ebi.ac.uk/intact/ws/interaction/findInteractions/thiS",
    params={"taxonId": 511145},
    timeout=30,
)
print(resp2.status_code)
print(resp2.text[:2000])
```

IscS-IscU、ThiS-ThiFのような既知ペアが実際にヒットするか、レスポンス
形式(JSON構造、どのフィールドに相互作用相手の識別子が入っているか)を
確認すること。APIキーが必要だった場合はその旨も報告してほしい。

### 3. BioGRIDのAPIキー取得とバルクファイルの中身確認

- `https://webservice.thebiogrid.org` で無料のAPIキーを取得できるか
  (ユーザー本人の対応が必要な場合は、その旨を報告して一旦保留でよい)。
- `BIOGRID-ORGANISM-5.0.261.mitab.zip`をダウンロードし、中に含まれる
  生物種ファイルの一覧を確認する。特に:
  - *E. coli*、*S. cerevisiae*(酵母)が含まれているか
  - アーキア(*Haloferax volcanii*, *Sulfolobus*, *Pyrococcus*,
    *Methanocaldococcus jannaschii*等)が1種類でも含まれているか
    (含まれていれば、より近縁な生物種からのinterolog推定ができ、
    信頼性が上がる)
  - zip全体のサイズ、および該当する生物種ファイル単体のサイズ

```bash
curl -sL -o /tmp/biogrid_organism.zip \
  "https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-5.0.261/BIOGRID-ORGANISM-5.0.261.mitab.zip"
unzip -l /tmp/biogrid_organism.zip | head -100
```

### 4. 実データでの検証(オルソログが見つかった場合)

上記1で大腸菌等が既存のBLAST参照ゲノムに含まれていた場合、以下の代表的な
クエリについて、実際にinterolog推定でどれだけ候補が拾えるか試算する:

- MA_4115(現行の主要クエリ)
- MA_0795 / MA_0826(radical SAM/Elp3型)
- MA_0074(RadA/RadB)
- MA_0266(RtcB)

それぞれについて、「大腸菌/酵母のオルソログが特定でき、かつそのオルソログが
IntAct/BioGRID上で既知の相互作用パートナーを持っている」候補が何件出るかを
数える。0件でも構わない(STRING/KEGGと同様、正直に報告すること)。

## 成果物として欲しいもの

`claude/phase6c_kegg_pathway_investigation.md`と同じ体裁のレポート
(Markdown、英語で構わない)。以下を含むこと:

1. 既存BLAST参照ゲノムの構成(どの生物種が既に使えるか)
2. IntAct API疎通確認の結果(実際のレスポンス例、認証の要否)
3. BioGRIDのAPIキー取得可否、バルクファイルの生物種カバレッジ
4. 上記4クエリでのinterolog候補の試算件数
5. 「実装する価値があるか」についての率直な所見。特に、STRING(既に
   実装済み)・KEGG(不採用)との比較で、実装コスト(オルソログマッピング
   ロジックの新規実装、場合によってはBLAST参照ゲノムの追加、BioGRIDの
   APIキー管理)に見合うシグナルの強さがあるかどうか。

## 非対応・スコープ外

- 実際のコード実装(`analysis/interolog_bridge.py`相当)は、この調査が
  完了しCowork側で結果を確認してから別途指示する。
- オルソログ判定の厳密な閾値設計(BLAST identity/coverage/evalueの
  基準をどう置くか)は今回の調査範囲外。既存の`ortholog_filter`設定
  (`config.yaml`のstrong/medium/weak閾値)がそのまま使える可能性が
  高いという前提でよい。
