# PIHブリッジ `interaction_priority_score` cap再較正(2026-09-24)

PIH連携ブリッジ(`analysis/pih_evidence_bridge.py`)Stage B1〜B3(実bundle動作確認・Tier Aフルカバレッジ化・
`interaction_score`スコープ拡大検証)の続き。`interaction_score`/`final_score`は`INTERACTION_SCORE_COMPONENT_NAMES`に
`pih_*`カテゴリを含まないため構造的に無関係(Stage B3で確定・現状維持)。一方、既定の`ranking_metric`である
`interaction_priority_score`(v2 evidence-basedの全カテゴリ加重合計)は`pih_cellular_compatibility`/`pih_evolutionary`の
影響を受け、旧cap(5/10)ではTier A/Bの分離をPIH無効時より悪化させていた。この記録はそのcap値の再較正。

## 結論

| 項目 | 旧値 | 新値 |
|---|---:|---:|
| `pih_cellular_compatibility`(cap) | 5 | **0** |
| `pih_evolutionary`(cap) | 10 | **3** |
| `pih_direct_interaction`(cap) | 20 | 20(**変更なし**、未curationのため今回のスコープ外) |

## データセット

- Tier A 10ペア全件(`claude/experimental_interactions_curated.csv`)、Tier B 5ペア(HdrD1/MtpAが関わる分)、
  AlphaFold3陰性28件(`claude/calibration/af3_negatives_MA_4115.csv`)。
- PIH側: 7クエリ(MA_4115 + HdrD1・MtpA・MtpC・NifD・NifK・AtwA)の`candidate_evidence_bundle.jsonl`
  (PSORTb/OrthoFinder再実行なし、既存の中間TSVを再利用して生成)。
- ProteinHunter側: `negative_hit`有効の標準Tier A比較設定(`candidates`/`candidates_relaxed`/`no_hit`/`negative_hit`、
  `max_candidates_per_query: 2600`)。

## グリッドサーチ(`interaction_priority_score`のAUC、Tier A優先・Tier B参考)

`pih_cellular_compatibility` ∈ {0, 2, 5, 8} × `pih_evolutionary` ∈ {0, 1, 2, 3, 5, 10} の24通り。
`analysis/scoring_engine.py::score_candidate`の実コードを直接呼び出し(独自再実装ではない)、
cap=(5, 10)で実際のパイプライン出力(`interaction_priority_score`列)と全15ペアで完全一致(誤差<0.0005)することを
検証した上で実施。

| cap_cc | cap_evo | AUC(Tier A) | AUC(Tier B) | Tier A vs. 無効時 | Tier B vs. 旧既定 |
|---:|---:|---:|---:|---:|---:|
| **0** | **3** | **0.9923** | **0.6346** | **+0.0077** | **+0.0654** |
| 0 | 5 | 0.9923 | 0.6269 | +0.0077 | +0.0577 |
| 0 | 10 | 0.9923 | 0.6038 | +0.0077 | +0.0346 |
| 2 | 10 | 0.9923 | 0.5962 | +0.0077 | +0.0269 |
| 0 / 0(PIH無効相当) | 0 | 0.9846 | 0.6346 | ±0.0000 | +0.0654 |
| 5 | 10(旧既定) | 0.9808 | 0.5692 | -0.0038 | ±0.0000 |
| 8 | 5 | 0.9692(グリッド最悪) | 0.5692 | -0.0154 | ±0.0000 |

全24通りの完全な表は本ディレクトリの`grid_results.txt`を参照。

**選定根拠**: 「PIH無効時の水準(Tier A 0.985・Tier B 0.635)に近い、または上回る」組み合わせのうち、
`cap_cc=0, cap_evo ∈ {1,2,3}` のみが両指標で無効時と同等以上(Tier B完全一致、Tier Aは上回る)。
その中でTier Aが最良の`cap_evo=3`を採用した。`cap_cc=0, cap_evo=5または10`も同じTier A AUCだが、
Tier Bがより低い(0.6269 / 0.6038)ため`cap_evo=3`を優先。

## 既知の限界(参考程度、断定しない)

- n=10(Tier A)・n=5(Tier B)・n=28(陰性)の小標本。近傍グリッドセル間の差は1〜2ペアの入れ替わりに相当し、
  「最適値」の断定はしない。
- `pih_cellular_compatibility=0`は、カテゴリを`INTERACTION_SCORE_COMPONENT_NAMES`やスコープから除外するのとは異なる:
  スコアへの寄与はゼロになるが、該当ペアで証拠がAVAILABLEなら`evidence_category_count`/`available_weight_total`には
  引き続きカウントされる(`analysis/scoring_engine.py::_score_categories`の"available_weight>0なら active"という
  判定はcapを見ない)。これはTier閾値・formal_score_available判定に極めて限定的な影響を与えうるが、今回の
  AUC検証はスコア値自体の分離のみを見ており、この副次効果は未検証(将来の課題)。
- `pih_direct_interaction`(cap 20)はTier A/B/陰性のいずれにも実データが流れておらず(`known_interactions`/
  `fusion`が未curation)、このグリッドサーチでは評価不能。今回は変更していない。

## 実データでの前後比較(パイプライン実行、`interaction_priority_score`)

同一config・同一bundleで、cap変更前(5/10)と変更後(0/3)のパイプラインを実行し比較(reimplementationではなく実出力):

| 指標 | 変更前 | 変更後 | 差分 |
|---|---:|---:|---:|
| `interaction_score`(全15ペア) | - | - | **完全一致(diff=0.0)** |
| `final_score`(全15ペア) | - | - | **完全一致(diff=0.0)** |
| `interaction_priority_score` AUC(Tier A, n=10, 陰性n=28) | 0.9821 | 0.9929 | +0.0107 |
| `interaction_priority_score` AUC(Tier B, n=5, 陰性n=28) | 0.5714 | 0.6321 | +0.0607 |

(グリッドサーチ表の値と実行のたびに完全一致しないのは、陰性データ26件 vs 28件など抽出時の母集団の違いによる。
方向性・改善幅は一致。)

## 再現

```bash
# PIH側(既存の中間TSVを再利用、7クエリ)
.venv/bin/python -m protein_interaction_hunter generate-candidates --config <7-query config>

# ProteinHunter側(negative_hit有効、pih_evidence_bundleに上記出力を指定)
python main.py --config <tier-a-full config>
```

`grid_results.txt`にグリッドサーチの生ログ(バリデーション+24通り全件)を保存。
