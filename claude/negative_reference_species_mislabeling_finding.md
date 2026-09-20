# 陰性参照ゲノムのラベル誤り — 調査記録

Status: **修正・影響評価・再発防止・派生した全事項が対応済み。マージ待ちの
みのPRが4本。** 全564テストパス。PR #17(ゲノム修正・thermoplasma復帰・
起動時チェック)・PR #19(patches/の7ファイル追加、MtpA/MtpC訂正含む)・
PR #20(GFF既定を全バケツtrueに戻す)としてpush済み。PR #18は#20採用に
伴いクローズ済み(コメント付き)。詳細な調査ログは
`claude/negative_reference_genome_mixup_investigation.md`、GFF既定変更の
前後比較は`claude/gff_annotation_all_buckets_default_comparison.md`
(いずれもリポジトリ内)。

残るオープンな作業(コード変更ではない): PR #19の`rockhopper_multisample_
validation.md`訂正注記に付けた「確度は高いが確定ではない」を、原著論文
(J Bacteriol 2019, doi:10.1128/jb.00130-19)の遺伝子座対応表そのものを
読める機会があれば「確認済み」に格上げできる(下記末尾参照)。急ぐ理由は
ないため、対応するPRのマージ判断とは切り離してよい。

## 確定した事実(再掲)

`data/databases/negative/Sulfolobus_solfataricus/`の中身は実際には
`GCF_002945325.1` = *Methanococcus maripaludis* strain DSM 2067で、
`positive/Methanococcus_maripaludis/`のprotein.faaとmd5が完全一致していた
(陽性参照と同一ファイルが陰性参照にも混入)。混入はコミット`4f02818`
(2026-07-02)。*Saccharolobus solfataricus* P2(`GCF_000007005.1`、
旧名Sulfolobus solfataricus)の正しいデータに差し替え済み、md5照合済み。

## 影響評価(フルパイプライン再実行・完了)

### config.yaml既定設定(5クエリ、candidate_sources.negative_hit=False)での前後比較

- 対応するTier A既知ペア13行のうち、Final Score表に載っていたのは
  修正前2行→修正後**11行**。誤ったデータではTier A真陽性の多くが
  Negative_hitバケツに落ちて(=既定でこのバケツはOFFなので)非表示に
  なっていたが、正しいデータではCandidates/Candidates_relaxedに入り、
  **negative_hitを有効化しなくても見えるようになった。**
- Tier2_Strong該当数: MA_4165 0→1、MA_3898 0→3、MA_3899 1→3。
- 副作用: `Candidates`バケツが127→709件に急増し、
  `max_candidates_per_query: 200`の切り捨てが効くようになった
  (MA_4115のFinal Score表に載るAF3陰性が28件中12→4件に減少 — バグでは
  なく、真の候補が大量に増えたことで上位200件の構成が変わったため)。
- [[conserved_query_visibility_design]]の警告は、この5クエリでは
  修正後は1件も出ない(全てweak/noneに転じたため、想定通り)。

### キャリブレーション相当の設定(全バケツ有効、上限5000、10クエリ)での比較

同一コード・同一config、陰性参照フォルダの中身だけ差し替えて比較
(最初の比較はGFF注釈が2バケツにしか効いておらずSTRING/GEOが歪んでいた
ため破棄し、全バケツGFF onで取り直し済み)。

| 指標 | 修正前 | 修正後 |
|---|---:|---:|
| Tier A の interaction_score 平均 | 22.48 | 22.48(不変) |
| AF3陰性 の interaction_score 平均 | 10.66 | 10.66(不変) |
| Tier A の final_score 平均 | 30.09 | **36.51** |
| Tier B の final_score 平均 | 35.77 | **44.10** |
| AF3陰性 の final_score 平均 | 25.14 | 25.14(不変) |
| Tier A vs AF3陰性 分離(AUC, final_score) | 0.563 | **0.684** |

- `interaction_score`(genomic_context+domain_complementarity、負の
  エビデンスを含まない)はnegative genomeに依存しないため完全に不変
  — 期待通りで、**既存の`experimental_interactions_calibration_report.md`
  のinteraction_scoreに関する結論(旧来のブレンド済みスコアより
  遥かに優れた分離を示す、という発見1)は今回の取り違えの影響を
  受けていない。**
- `final_score`はNegative_hit(strong等)だったペアで約8点上昇
  (例: MtpC 45.9→54.3、Tier3→Tier2)。AF3陰性はNegative_hit扱いでは
  ないため変化なし。**分離(AUC)も0.563→0.684に改善しており、
  取り違えを直したことで、むしろ元のPR #9-11の仮説(negative_hit系
  ペナルティを外した設計判断)がより支持される結果になった。**
- n が小さいため AUC は参考値。また今回の再実行値は元の
  `experimental_interactions_calibration_report.md`の絶対値
  (Tier A平均39.55等)とは直接比較できない — その後にコード自体が
  変わっている(Final Score統合等)ため、**有効な比較は「同一コード・
  同一configでの前後」のみ。**

## 再発防止(実装完了)

- `thermoplasma_acidophilum`: ローカルの`protein.faa`をNCBIの
  `GCF_000195915.1`と照合しmd5一致を確認、`.gitignore`から除外を解除、
  git管理に復帰(PR #17)。
- `config/reference_genomes.v1.yaml`: 期待する陽性・陰性ゲノムの
  フォルダ名/accession/生物種を記載したマニフェストを新設(PR #17)。
- `core/reference_genomes.py` + `main.py`起動時チェック: マニフェストと
  実データを突き合わせ、(1)期待フォルダの欠落 (2)リスト外フォルダ
  (3)accession不一致 (4)FASTA生物種の不一致 (5)別フォルダ間でのFASTA
  完全一致 (6)陽性/陰性への同一種混入、の6種を検出しWARNINGログを出す
  (実行は止めない)。各ゲノムのaccession/organism/protein数/md5を
  毎回INFOログに残す。今回の取り違えを4件の警告として再現できること、
  修正後は"Reference genome check passed"になることを確認済み。
  `tests/test_reference_genomes.py`で各検出種別・実データとの整合を
  テスト済み(PR #17)。
- この起動時チェックは[[run_provenance_design]]の「スコープ拡張の検討」
  で提案したデータ構成フィンガープリントの実質的な実装にあたる。
  今後run provenance機能を実装する際は、このマニフェスト/チェック結果を
  再利用し、ロジックを重複させないこと。

## 未決定事項1(対応済み): calibration reportとconfig.yamlコメントの書き直し

PR #17で対応完了。

- `claude/experimental_interactions_calibration_report.md`: 冒頭に訂正
  ブロックを追記(本文は当時の記録としてそのまま保持、
  `implementation_status_scoring_v2.md`のMA_RS16635/MA_RS20505の訂正と
  同じ流儀)。17パートナー中、Candidatesが0→14、Negative_hitが17→3
  (DnaK・Mer・Hsp20)に変わったことを明記。interaction_scoreへの影響は
  ないことも明記。
- `config.yaml`の`candidate_sources`コメント: 一般的な助言(保存性の
  高いクエリではnegative_hitの有効化を検討)は残しつつ、HdrD1-Mer等の
  具体例への言及は「要検証」+取り違えの説明に置き換え。
- 付随して、[[conserved_query_visibility_design]]の警告文言と当該設計書
  自身がこの具体例を実証済みの根拠として引用していた箇所も、PR #16で
  同様に訂正済み(依頼になかった追加対応、仕組みと影響だけを述べる
  表現に変更)。

## 未決定事項2(対応済み): GFF注釈バケツ差

### 性能面の確認結果

`_ANNOTATION_TARGETS_GFF_ALWAYS_ON`が`candidates`/`candidates_relaxed`
の2バケツだけなのは、性能上の制限ではなかった。コミット`beec3c8`で
「設定.xlsx」の既定プリセットをそのまま反映したものであり、それ以前は
全シートがtrueだった。実測でもGFF注釈は4,627件で0.21秒
(既定の1,261件では0.20秒)とほぼ無視できるコストであり、2バケツ限定に
積極的な理由はなかったことが確認された。

### 決定: PR #20(全バケツGFF既定trueに戻す)を採用、PR #18はクローズ済み

`no_hit`(既定ONの主要バケツ)にもnegative_hit系と同種の問題
(BLASTヒットが無い新規性の高い候補が既定でSTRING/GEOエビデンスを
持てない)があることが判明したため、バケツごとの個別連動(PR #18)を
積み増すより、`_ANNOTATION_TARGETS_GFF_ALWAYS_ON`をbeec3c8以前の
「全バケツtrue」に戻す方針(PR #20)を採用した。PR #18は#20に置き換える
コメントを付けてクローズ済み。

前後比較(同一コード・同一config、クエリ5件、陰性ゲノムは修正済み):

| 項目 | 前 | 後 |
|---|---|---|
| Candidates / Candidates_relaxed | 基準 | 完全に同一(行・値とも不変) |
| No_hitで`old_locus_tag`が付いた行(1,000行中) | 0 | 719 |
| No_hitでSTRING/GEO証拠が付いた行 | 0 | 712 |
| No_hitの平均証拠カテゴリ数 | 3.00 | 4.39 |

- 主要な候補選定ロジック(Candidates/Candidates_relaxed)には一切影響なし。
- No_hitプールは1,000行中851行が入れ替わる(以前は全員同じ3カテゴリで
  同点が多く、選別が実質的に恣意的だった)。
- MA_0725・MA_2514は、証拠3カテゴリ時のfinal_score 46.1が、5カテゴリでは
  37.0/35.6に下がりtop10から外れる — **エビデンス不足がスコアを
  過大評価していた具体例**(negative_hit系で見つかったのと同種の問題が
  `no_hit`にも実在していたことの直接証拠)。
- top10: MA_4115/MA_0688/MA_4165は同一、MA_3898/MA_3899は8/10が一致。
  Tier1-2の件数は5クエリとも不変。
- 未評価: クエリ数が5件と少なく、`no_hit`の真陽性検出力自体が
  向上したかどうかまでは評価していない。

## 副次的な整理(対応済み): `patches/`配下の未追跡ファイル7件

PR #19で対応完了。`claude/phase6c/6d/6e_*.md`と`claude_code_instructions_*`
4件を、削除せず追記訂正の方針で追加。

- **6c**: 候補バケツ件数(Candidates 127/1,261など)が取り違えたゲノムでの
  出力だったため、修正後の値(709/1,668)を注記として追記。結論(全クエリが
  PATHWAYを持たず共起計算不能)はバケツサイズに依存しないため有効なまま。
- **6d**: 参照ゲノム一覧に、Sulfolobusフォルダの中身が実はM. maripaludis
  だった旨の注記を追記。E. coli/酵母が必要という結論には影響なし。
- **6e・指示書4件**: 矛盾なし、原文未変更。

### 解決: MtpA/MtpCの遺伝子座割り当てが1ファイルだけ逆だった問題

`claude_code_instructions_rockhopper_multisample_validation.md`だけ
MA_4164=MtpA、MA_4165=MtpCとなっており、`config.yaml`・キャリブレーション
報告書・operon指示書(MA_4165=MtpA)と逆だった。GFFに遺伝子名が無く判断
できなかったため、以下の根拠で多数派に合わせて訂正済み(7行目・20行目に
日付入り訂正注記、本文は当時の記録として不変)。

- UniProt: MA_4165 → 登録遺伝子名**cmtA**("Methylcobalamin:CoM
  methyltransferase isozyme A", Q8TII4)。MA_4164 → 遺伝子名の登録なし、
  "Corrinoid protein"とのみ登録(Q8TII5)。
- Europe PMC経由で読めた一次文献(J Bacteriol 2019、
  doi:10.1128/jb.00130-19)の要旨: MtpAは精製されたメチルトランス
  フェラーゼ、MtpCはコリノイドタンパク質、と明記。UniProtの機能情報
  (MA_4165=メチルトランスフェラーゼ、MA_4164=コリノイドタンパク質)と
  一致し、MA_4165=MtpAの裏付けが強まった。
- 出版社サイト(403)・PMC(本文XML取得不可)経由での遺伝子座番号
  対応表そのものの直接確認は取れていない。注記・PRコメントには
  「確度は高いが確定ではない」と明記済み。論文本文を読める環境が
  あれば、対応表で最終確認し注記を「確認済み」に格上げできる
  (優先度低、急ぐ必要なし)。
- 本文中でMA_4164/MA_4165という遺伝子座番号そのものを直接使っている
  既存の手順・コードは、この訂正の影響を受けない
  (影響を受けるのはMtpA/MtpCという「ラベル」の解釈のみ)。
