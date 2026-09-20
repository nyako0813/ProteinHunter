# スコアリング再較正(2026-09-21)

`tools/scoring_recalibration.py`による分析と、その結果の採否の記録です。最適化・グリッドサーチはしていません
(Tier A n=11、AF3陰性 n=28 は小標本のため)。候補は事前に決めた少数だけを比較しました。

## 結論

| 項目 | 結果 |
|---|---|
| Tier3閾値(`tiers.tier3_min_score`) | **25 → 35 に変更** |
| Tier1(70)・Tier2(50)の閾値 | 変更なし |
| `external_ppi_evidence` cap 15 / `coexpression_evidence` cap 12 / `coexpression_gse64349`重み 1/3 / Final Scoreのcap 30/70 | **変更なし**(baselineを明確に上回る候補が無かった) |

## Phase 1: Tier閾値

- 定義(コード): `analysis/scoring_engine.py::_classify_tier`。Tier1 ≥70かつ3カテゴリ、Tier2 ≥50かつ2カテゴリ、Tier3 ≥25かつ1カテゴリ、
  それ未満は`Tier4_Weak`、スコアなしは`Unclassified`。`final_score_tier`・`evidence_tier`・`interaction_evidence_tier`の**3列が同じ`TierThresholds`を共有**する。
- `final_score`の分布: Tier A(n=11)は最小19.7・中央値44.3・最大62.5、AF3陰性(n=28)は最小9.7・中央値23.8・最大49.0。AUC 0.896。
- 基準: Youden's Jが最大になるのは`final_score ≥ 39.7`(J=0.802)。「既知の陽性の下限のすぐ下」基準では19.7だが、
  その水準(t=20)では陰性の57%が通ってしまう(J=0.34)。Tier3は「弱い証拠がある」という下限の意味合いなので、
  Jが最大の付近で**Tier A 11件中10件を残し、かつ陰性の通過率が14%まで下がる35**を採用した(t=35: 感度0.91・偽陽性率0.14・J=0.77)。
  Youden最適の40は、Tier Aの通過が10件から9件に減るため採らなかった。
- 詳細表: `recalibration_summary.md`。

### 副作用の実測(実パイプラインを新旧の既定で再実行して比較。27,756行)

スコア・順位は完全に同一。変わるのはTier3→Tier4のラベルだけ。

| 列 | Tier3の行数 旧→新 | Tier A(n=11)のTier3以上 | 陰性(n=28)のTier3以上 |
|---|---|---|---|
| `final_score_tier` | 4,305 → 194 | 10 → 10 | 12 → 4 |
| `evidence_tier` | 11,730 → 1,532 | 11 → 11 | 13 → 4 |
| `interaction_evidence_tier` | 99 → 27 | 9 → 4 | 3 → 1 |

- `select_top_candidates_per_query()`(Tier1/Tier2は順位に関係なく常に表示): 対象はTier1/Tier2だけで、その閾値は変えていないため、
  表示される候補は**旧と完全に同一**(上位15件/クエリで計90行、順位圏外の追加表示は全クエリで0)。
- **注意(要判断)**: `interaction_evidence_tier`は`interaction_score`(相互作用カテゴリだけの副スコア)由来で、Tier3のTier Aが
  7件中5件Tier4に落ちる。順位・表示対象には使われず参考列だが、感度は下がる。この列だけ旧値を保ちたい場合は、
  `scoring_engine_config`の`tiers.tier3_min_score: 25`で全体を旧値に戻せる(3列は一括で動く)。列ごとに分ける設計変更はしていない。
- Word/Notion/Excelのラベル表示: Tier文言は`tier_opening.*`(Tier3/Tier4で文面が違う)から出るだけで、Tier3/Tier4を条件にした
  絞り込みは無い。既存のラベル表示テストはすべて通る(790件)。

## Phase 2: cap/重みの候補(leave-one-out)

Tier Aは10ペア(一方向ずつではなく無向ペアを1 foldとして除外)、陰性28件。評価指標は`final_score`のAUC。
事前に決めた候補(YAMLは`variants/`):

| 候補 | 内容 |
|---|---|
| A | `coexpression_gse64349`の重み 1/3 → 1/2 |
| B | Final Scoreのcap 30/70 → 20/80 |
| C | `external_ppi_evidence` 15→20、`coexpression_evidence` 12→16 |
| D | `genomic_context` 25→20 |

結果(`final_score`、LOO平均AUC ± jackknife SE、baseline 0.896 ± 0.052): A 0.896、B 0.893、C 0.896 ± 0.071、D 0.889。
どの候補も「ΔAUC ≥ 0.02 かつ z ≥ 2.4 かつ差のあるfoldの75%以上で改善」を満たさない。入れ子選択(各foldで他のペアだけから最良候補を選ぶ)も
baselineの0.886に対して0.861で、選択しても得にならない。**変更なし**。

- 候補Aは構造上`final_score`/`interaction_score`を動かせない: `coexpression_gse64349`は`interaction_score`内の
  `coexpression_evidence`カテゴリのただ1つの構成要素で、カテゴリ内の重みは利用可能な重みの合計で割られて打ち消される。
  効くのは`interaction_priority_score`だけだが、その指標は現在のTier Aと陰性でAUC 0.987に飽和しており(全候補で同一)、判別に使えない。
- `interaction_priority_score`は順位付けに使われる既定の指標なので、飽和は「この小標本では区別できない」という意味であり、
  「最適」という意味ではない。

## 再現

```
PYTHONPATH=. python tools/scoring_recalibration.py --curated claude/experimental_interactions_curated.csv \
  --negatives claude/calibration/af3_negatives_MA_4115.csv --baseline baseline=rows_baseline.csv \
  --candidate A=rows_A_gse64349_half.csv --candidate B=rows_B_final_20_80.csv \
  --candidate C=rows_C_ppi_coexp_up.csv --candidate D=rows_D_genomic_down.csv --out out --tier3 35
```

各`rows_*.csv`は、`variants/engine_*.yaml`(候補Aだけは`harness.py`が`V2_COMPONENT_WEIGHTS["coexpression_gse64349"]`を差し替える)で
パイプラインを実行し、`harness.py`で最終スコアの行を書き出したもの(各約2.5MBなのでリポジトリには含めていない)。
