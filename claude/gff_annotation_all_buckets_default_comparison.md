# annotation_targets.gff を全バケツ既定 true に戻す — 前後比較

Status: 変更提案(PR)。既定の出力が変わるため、前後比較を記録する。

## 経緯

`annotation_targets.<bucket>.gff` の既定は、`beec3c8`(2026-09-04、`consider_cross_species_matches`
導入)より前は全バケツ true だった。同コミットが `設定.xlsx` の既定プリセットを反映して
candidates / candidates_relaxed 以外を false にしたが、性能上の制限ではなかった
(GFF 注釈は辞書引きで、全4,627レコードでも約0.21秒、既定の1,668件で0.20秒)。
GFF 注釈は `old_locus_tag` を付与する処理で、STRING/GEO の証拠は `old_locus_tag` で
引くため、gff=false のバケツの候補は STRING/GEO 証拠が静かに MISSING になる。
`no_hit` は既定で候補ソースが ON の主要バケツなので、既定実行に影響する。

変更: `ANNOTATION_TARGET_DEFAULTS`(全バケツ gff=true)をそのまま既定にし、
`_ANNOTATION_TARGETS_GFF_ALWAYS_ON` と preset 連動を削除。config.yaml の明示値は従来どおり優先。

## 比較条件

同一コード基盤・同一 config(リポジトリの config.yaml + クエリ5件 MA_4115 / MA_0688 /
MA_4165 / MA_3898 / MA_3899、候補ソースは既定の Candidates / Candidates_relaxed / No_hit、
`max_candidates_per_query: 200`)。陰性参照は修正済みの Saccharolobus solfataricus を使う
(取り違え修正 #17 とは独立に比較するため、両方の実行で同じ陰性ゲノムを使用)。
前: 変更前の挙動(GFF 注釈対象 1,668 レコード)。後: 変更後(4,627 レコード)。
Excel は書かず最終スコア行(`02_Final_Score` と同じ行)を比較した。

## 結果

### 証拠が付く候補

| バケツ(行数は前→後) | old_locus_tag あり | STRING/GEO 証拠あり | 平均証拠カテゴリ数 |
|---|---|---|---|
| Candidates(1,000→1,000) | 956→956 | 950→950 | 4.88→4.88 |
| Candidates_relaxed(227→227) | 227→227 | 227→227 | 5.00→5.00 |
| **No_hit(1,000→1,000)** | **0→719** | **0→712** | **3.00→4.39** |

Candidates / Candidates_relaxed は行も値も完全に同一。変化は No_hit のみ。

### No_hit の候補プール

`No_hit` は「クエリごとに priority score 上位200件」をスコアする。証拠が付くとその選別が
変わり、No_hit 1,000行のうち **851行が入れ替わる**(両方に残るのは149行)。
以前は No_hit 全員が同じ3カテゴリしか持たず、スコアが強く同点になっていたため、
上位200件の選別は実質的に恣意的だった。

### スコアの変化(前後の両方に残った No_hit 149行)

- 113行で interaction_score / final_score が変化。interaction_score の変化量は
  平均 +8.0(範囲 -15.1〜+16.6)、final_score は平均 +5.6(範囲 -10.5〜+11.6)。
- 証拠が付いた行の多くは上がり、証拠が付いた結果かえって下がる行もある。
- No_hit 全体(プール入れ替えを含む): interaction_score 平均 3.92→9.98、
  final_score 平均 14.97→19.21、final_score 最大 60.1→50.2。
- **証拠が欠けているとスコアが有利に出ていた**: 3カテゴリしか無い No_hit 候補が、
  証拠が付くと下がる例。MA_0725 / MA_2514(MA_3898/MA_3899 クエリで final_score
  46.1、3カテゴリ)→ 37.0 / 35.6(5カテゴリ)、順位は top10 圏外(13〜17位)へ。
  MA_0687(HdrD1 = MA_0688 の隣接遺伝子)は 60.1→50.2 だが同クエリの1位のまま。

### ランキングへの影響(クエリ別)

| クエリ | Tier1-2 件数 | top10 の重なり | 備考 |
|---|---|---|---|
| MA_4115 | 0→0 | 10/10 | 変化なし |
| MA_0688 | 1→1 | 10/10 | No_hit の MA_0687 が 1 位のまま(60.1→50.2) |
| MA_4165 | 1→1 | 10/10 | 変化なし |
| MA_3898 | 3→3 | 8/10 | No_hit 2件(MA_0725, MA_2514)が top10 から抜け、Candidates 系が繰り上がる |
| MA_3899 | 3→3 | 8/10 | 同上 |

全体の Tier: Tier2_Strong 8→8、Tier3_Moderate 771→923、Tier4_Weak 1,448→1,296
(No_hit の入れ替えとスコア変化による)。

## 解釈と注意

- Candidates / Candidates_relaxed の結果は変わらないので、既存の主要な結論には影響しない。
- 変わるのは No_hit の扱い。証拠が欠けたまま高得点に見えていた No_hit 候補が
  是正される一方、No_hit のスコア対象プールが大きく入れ替わる(既定で ON のバケツ
  なので、既定実行の出力が変わる)。
- n(クエリ5件)は小さい。No_hit 内での真陽性の増減は本比較では評価していない。
- `negative_hit` 系だけを条件連動で対応する案(PR #18)は、この変更で冗長になるため
  クローズした。
