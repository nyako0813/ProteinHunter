# Word出力 日本語化 — 設計書

Status: **設計中(未実装)。実装前にClaude Codeへこの設計書全文を渡すこと。**

## 背景・要望

ユーザーから: 「出力のwordファイルは日本語にできますか？」。あわせて
[[notion_export_design]] でNotion出力を第三の出力として追加する計画がある
ため、日本語化はWord単体ではなく、両出力が共有するナラティブ生成層に
実装し、将来Notionにもそのまま反映されるようにする。

## 方針

- **英語のまま維持するもの**: 遺伝子座番号(locus tag、例: MA_4165)、
  遺伝子名・タンパク質名(gene symbol、例: MtpA、cmtA)、スコア列名
  (`interaction_score`, `final_score`など。Excelとの対応を保つため)、
  `candidate_source`/`Evidence_Tier`/`negative_hit_strength`などの
  内部分類ラベルの**データ値そのもの**(Excelの列値と完全一致させる)。
- **日本語化するもの**: 見出し・ラベル文言(例: "Final Score:" →
  "最終スコア:")、自動生成される説明文・接続詞・評価コメントの文章表現、
  節タイトル(例: "Candidate Details" → "候補の詳細")。
- 分類ラベルは「英語表記(日本語)」の併記とする(例:
  "Evidence Tier: Strong(強)")。これによりExcelとの用語一致を壊さずに
  可読性を上げる。対応表は本書末尾の「訳語対応表」に固定する。

## 実装方法

- `config.yaml` に `report_language: "en"`(既定値。後方互換のため既定は
  英語のまま)を追加し、`"ja"`を指定すると日本語出力に切り替わる仕組みに
  する。
- 新規 `output/report_i18n.py` に、見出し・固定文言の `en`/`ja` 辞書と、
  分類ラベルの訳語対応表(下記)を持たせる。
- `output/word_narrative.py` の候補ごとの説明文組み立て処理を、
  [[notion_export_design]] で述べた中間表現(`NarrativeSection`のリスト:
  見出し/段落/箇条書き、テキストは`report_i18n`経由で言語切り替え)を
  返す形にリファクタリングする。
- `output/word_report.py` は、この中間表現をdocxにレンダリングする側に
  徹し、文言のハードコードを`report_i18n`参照に置き換える。
- 出力ファイルは言語ごとに分けない。`report_language`設定1つで、生成
  される1つのWordファイルの言語が決まる(Excelは対象外・英語のまま
  変更しない — 列名の一貫性を壊さないため)。

## 訳語対応表(案。実装時に固定)

| 英語(データ値・変更しない) | 日本語表記(説明文中でのみ併記) |
|---|---|
| Strong | 強 |
| Medium | 中 |
| Weak | 弱 |
| Candidates | 候補 |
| Candidates_relaxed | 候補(緩和条件) |
| No_hit | ヒットなし |
| Negative_hit | 陰性ヒット |

(候補ソース・エビデンス階層の全ラベルを実装時に洗い出して追記する)

## スコープ外

- Excel出力の日本語化(列名・シート名は変更しない。データ突き合わせ・
  [[calibration_pipeline_expansion_design]]との整合性を壊すため)。
- 論文由来の遺伝子名・タンパク質名自体の翻訳(そのまま英語表記を維持)。

## テスト計画

- 既存のWord出力テスト(英語、`report_language`未指定/`"en"`)は変更なく
  パスすることを確認(回帰)。
- `report_language: "ja"` 時に、見出し・固定文言が辞書通りに置き換わり、
  遺伝子座番号・スコア数値・分類ラベルの英語表記(データ値部分)が変わらない
  ことを検証する新規テストを追加する。

## 見積もり・実装順序

推奨実装順序: **本設計(ナラティブの中間表現化・i18n辞書導入・日本語Word
出力)を先に実装し、その上に [[notion_export_design]] のNotion出力を
実装する。** 中間表現とi18n辞書はNotion出力からもそのまま再利用できる
ため、二重実装を避けられる。

新規: `output/report_i18n.py`。変更: `output/word_narrative.py`
(中間表現へのリファクタリング)、`output/word_report.py`(文言参照の
置き換え)、`config.py`(`report_language`設定の追加)。既存の
スコアリング・Excel出力コードには手を入れない。