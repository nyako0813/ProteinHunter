"""Report wording in English and Japanese (Word report today, Notion export later).

Two things live here, both plain data with no ``python-docx`` or scoring
knowledge:

* ``STRINGS``: every heading, label and narrative sentence of the report, keyed
  by a stable name, once per language. ``en`` is byte-for-byte the wording the
  report always had (it is the default and the regression baseline); ``ja`` is
  the Japanese rendering. Both use the same ``str.format`` placeholders.
* ``LABEL_JA``: Japanese glosses for classification labels (candidate sources,
  tiers, negative-hit strengths, evidence categories). In Japanese prose a label
  is written ``English(日本語)`` via :func:`bilingual`, so the English data value
  that also appears in the Excel workbook stays visible. Data values themselves
  (locus tags, scores, ``candidate_source`` in tables, tier names) are never
  translated.

Conventions for the Japanese text (claude/word_report_japanese_localization_design.md):
plain "である" style, hedged wherever the English is hedged ("consistent with"
-> "整合的であるが ... 立証するものではない"), no sentence asserting that a
candidate *is* the interaction partner.
"""

from __future__ import annotations

from typing import Any

SUPPORTED_LANGUAGES: tuple[str, ...] = ("en", "ja")
DEFAULT_LANGUAGE = "en"

_EN: dict[str, str] = {
    # -- title page --------------------------------------------------------
    "report.title": "ProteinHunter Candidate Report",
    "report.generated": "Report generated: {timestamp}",
    "report.scoring_model": "Scoring model: {scoring_model}",
    "report.code_version": "Code version: {version}",
    "report.config_fingerprint": "Config fingerprint: {fingerprint}",
    "report.sidecar": "Full effective configuration saved alongside this report as {location}",
    "report.sidecar_unknown_location": "<Excel workbook stem>.run_provenance.yaml",
    "report.reference_genome_check": "Reference genome check: {text}",
    "report.reference_genome_warning": "WARNING - see run log for details (core/reference_genomes.py)",
    "report.queries_evaluated": (
        "Queries evaluated: {n_queries}. Candidates shown per query: up to "
        "{max_per_query}, plus any additional Tier1_VeryStrong/Tier2_Strong "
        "candidate regardless of rank."
    ),
    "report.toc_placeholder": "Right-click and choose “Update Field” to generate the table of contents.",
    "sentence_joiner": " ",
    # -- section 5: evidence architecture ---------------------------------
    "arch.title": "5. Evidence Architecture",
    "arch.intro": (
        "This run used the {scoring_model} scoring model. The seven "
        "evidence categories below are the pipeline's full evidence "
        "vocabulary; which ones actually contributed evidence for any "
        "given candidate depends on what data and configuration were "
        "available for this specific run (see each candidate's own "
        "“why this candidate ranks highly” text in section 8)."
    ),
    "arch.seq.heading": "5.1 Sequence Evidence (cap {cap:.0f})",
    "arch.seq.body": (
        "BLAST-based positive/negative classification and best-hit "
        "identity/coverage/E-value strength. Populated for essentially "
        "every candidate that has any BLAST hit at all -- this is the "
        "pipeline's most consistently available evidence category."
    ),
    "arch.func.heading": "5.2 Functional/Domain Evidence (cap {cap:.0f})",
    "arch.func.body": (
        "Shared or complementary functional annotation (CDD/Pfam domains, "
        "description terms) between query and candidate. Populated "
        "whenever domain annotation is enabled and available for both "
        "proteins."
    ),
    "arch.genomic.heading": "5.3 Genomic Context (cap {cap:.0f})",
    "arch.genomic.body": (
        "Genomic proximity between query and candidate genes, from GFF "
        "coordinates -- used as positive evidence only (a distant "
        "candidate is never penalized for being far away, only not "
        "credited for being close). Populated whenever GFF neighborhood "
        "data is available for both genes."
    ),
    "arch.interaction.heading": "5.4 Interaction Evidence (cap {cap:.0f})",
    "arch.interaction.body": (
        "External protein-protein interaction evidence from up to three "
        "independent sources: STRING PPI, GEO transcript coexpression, "
        "and the optional ProteinInteractionHunter (PIH) direct-interaction "
        "bridge. Populated whenever the corresponding optional data source "
        "(STRING taxon ID, GEO coexpression, or a PIH evidence bundle) is "
        "configured for the run; each source is independent and any subset "
        "may be available."
    ),
    "arch.evolutionary.heading": "5.5 Evolutionary Evidence (cap {cap:.0f})",
    "arch.evolutionary.body": (
        "Phylogenetic/evolutionary profile consistency between candidate "
        "and query, sourced entirely from the optional PIH bridge. {closer}"
    ),
    "arch.cellular.heading": "5.6 Cellular Compatibility (cap {cap:.0f})",
    "arch.cellular.body": (
        "Subcellular localization / compatibility evidence, also sourced "
        "from the optional PIH bridge. {closer}"
    ),
    "arch.negative.heading": "5.7 Negative Evidence (reserved)",
    "arch.negative.body": (
        "Reserved for evidence that directly contradicts a candidate/query "
        "pairing -- e.g. incompatible cellular localization, phylogenetic "
        "inconsistency, or functionally contradictory annotation. No such "
        "signal is implemented in this pipeline version. This category is "
        "reported as not evaluated for every candidate, in every run, with "
        "no exceptions. This is a deliberate design decision, not an "
        "oversight: an earlier implementation attempt used "
        "negative_hit_strength (shown elsewhere in this report, under "
        "candidate_source) as a stand-in for this category and found, on "
        "real-data verification, that it measures a different thing -- how "
        "broadly a candidate protein is conserved across negative "
        "reference genomes -- and using it as a contradiction penalty "
        "incorrectly punished well-conserved true interaction partners. "
        "The two signals are kept visibly separate in this report for "
        "that reason, not merged back together for convenience."
    ),
    # -- section 7: ranking ------------------------------------------------
    "ranking.title": "7. Candidate Ranking",
    "ranking.empty": "No query-specific candidates were produced by this run.",
    "ranking.query_heading": "7.{index} Query: {query_id}",
    "ranking.col.rank": "Rank",
    "ranking.col.candidate": "Candidate",
    "ranking.col.final_score": "Final Score",
    "ranking.col.tier": "Tier",
    "ranking.col.source": "Candidate Source",
    # -- section 8: details ------------------------------------------------
    "details.title": "8. Candidate Details",
    "details.empty": "No query-specific candidates were produced by this run.",
    "details.query_heading": "8.{index} Query: {query_id}",
    "details.why_label": "Why this candidate ranks highly: ",
    "details.interpretation_label": "Biological Interpretation: ",
    "details.excel_ref": (
        "Full data for this candidate: see {excel_filename}, sheet "
        "02_Final_Score, query_id={query_id}, "
        "candidate_protein_id={candidate_id}."
    ),
    # -- narrative: why this candidate ranks highly ------------------------
    "tier_opening.Tier1_VeryStrong": (
        "This candidate reached the highest confidence tier for this query "
        "(Tier 1 — Very Strong), with a Final Score of {final_score:.1f}/100 "
        "supported by evidence from {evidence_category_count} independent "
        "evidence categories."
    ),
    "tier_opening.Tier2_Strong": (
        "This candidate reached the second-highest confidence tier for this "
        "query (Tier 2 — Strong), with a Final Score of {final_score:.1f}/100 "
        "supported by evidence from {evidence_category_count} independent "
        "evidence categories."
    ),
    "tier_opening.Tier3_Moderate": (
        "This candidate reached a moderate confidence tier for this query "
        "(Tier 3 — Moderate), with a Final Score of {final_score:.1f}/100. "
        "The supporting evidence is present but comes from fewer independent "
        "categories ({evidence_category_count}) than higher-tier candidates."
    ),
    "tier_opening.Tier4_Weak": (
        "This candidate's Final Score of {final_score:.1f}/100 places it in "
        "the weakest confidence tier evaluated (Tier 4 — Weak). Its position "
        "in this list reflects a relative rank among the candidates scored "
        "for this query, not strong independent support."
    ),
    "tier_opening.unclassified": (
        "This candidate did not receive a formal Final Score for this query "
        "(insufficient evidence to meet the pipeline's minimum evidence "
        "requirement). It appears here only because it was among the "
        "higher-ranked candidates by whatever partial evidence was available; "
        "treat its inclusion as informational, not as a ranked, scored result."
    ),
    "why.contributing": "Contributing evidence categories: {listing}.",
    "why.listing_separator": ", ",
    "why.none": (
        "No individual evidence category scored above zero for this "
        "candidate under the current run's configuration."
    ),
    "source.Candidates": (
        "This protein was classified as a strict positive candidate — a "
        "BLAST hit to a positive reference sequence with no negative-reference "
        "hit at all, the most stringent classification this pipeline produces."
    ),
    "source.Positive_all_sources": (
        "This protein hit every configured positive reference source, with "
        "no negative-reference hit."
    ),
    "source.Candidates_relaxed": (
        "This protein was classified as a positive candidate under relaxed "
        "criteria — it has a positive-reference hit, and tolerates a "
        "non-strong (medium or weak) negative-reference hit."
    ),
    "source.No_hit": (
        "This protein had no BLAST hit to either the positive or negative "
        "reference sets. It is a lineage-specific or novel candidate not "
        "detected by sequence homology alone — its ranking here depends "
        "entirely on non-sequence evidence (genomic context, domain "
        "annotation, interaction evidence), not on a BLAST match."
    ),
    "source.Negative_unmatched": (
        "This protein has no negative-reference hit, though it did not meet "
        "the stricter positive-match criteria above."
    ),
    "source.Negative_hit": (
        "This protein has at least one negative-reference BLAST hit "
        "(strength: {negative_hit_strength}). It is included in this "
        "ranking despite that hit because of the other evidence enumerated "
        "above; interpret its position with additional caution."
    ),
    "negative_hit_caveat": (
        "Note: this candidate also has a {negative_hit_strength} BLAST hit to a "
        "negative-reference sequence, which does not change its candidate_source "
        "classification but is a relevant caveat."
    ),
    # -- narrative: biological interpretation ------------------------------
    "interpretation.opening": (
        "Based on the evidence currently available to this pipeline, "
        "{candidate_id} ranked {rank} of {n_candidates} evaluated "
        "candidates for query {query_id}. This reflects the evidence "
        "categories described below as currently available to the "
        "pipeline — it is not a confirmed identification of "
        "{candidate_id} as the target enzyme or interaction partner, and "
        "should be read as 'best-supported by currently available "
        "evidence,' not as a settled conclusion."
    ),
    "color.genomic_context": (
        "Its genomic proximity to the query gene is consistent with — but "
        "does not by itself establish — a shared operon or functional "
        "module."
    ),
    "color.functional_domain": (
        "Shared or complementary functional domain annotations were found "
        "between this candidate and the query, consistent with a related or "
        "complementary biochemical role; domain similarity alone does not "
        "establish direct interaction."
    ),
    "color.interaction": (
        "External interaction evidence (protein-protein interaction "
        "database and/or coexpression data) supports a functional "
        "association between this candidate and the query, independent of "
        "sequence similarity or genomic position."
    ),
    "color.evolutionary": (
        "Evolutionary/phylogenetic profile evidence from the optional PIH "
        "bridge supports consistency between this candidate and the query "
        "across the genomes compared."
    ),
    "color.cellular_compatibility": (
        "Cellular-compatibility evidence from the optional PIH bridge "
        "supports this candidate and the query being compatible with the "
        "same subcellular context."
    ),
    "negative_evidence_closer": (
        "Negative (biological contradiction) evidence is a reserved category "
        "with no implemented signal in this pipeline version (see Evidence "
        "Architecture, section 5.7) — its absence here does not mean this "
        "candidate was checked for contradicting evidence and cleared."
    ),
    "interpretation.closer": (
        "Categories not listed above were not evaluated for this candidate "
        "in this run — a blank category means no evidence was available "
        "to assess, not that the evidence was checked and found absent. "
        "{evolutionary_closer} {negative_evidence_closer}"
    ),
    "evolutionary_closer.legacy": (
        "This run used the legacy_additive scoring model, which does "
        "not evaluate Evolutionary or Cellular Compatibility evidence "
        "at all — those categories exist only under the "
        "v2_evidence_based scoring model."
    ),
    "evolutionary_closer.pih": (
        "Evolutionary and Cellular Compatibility evidence were "
        "evaluated for this run using a supplied ProteinInteractionHunter "
        "(PIH) data bundle; a blank value for either means no matching "
        "evidence was found for this specific candidate, not that the "
        "category was skipped."
    ),
    "evolutionary_closer.no_pih": (
        "Evolutionary and Cellular Compatibility evidence require an "
        "optional external data bundle (ProteinInteractionHunter/PIH) that "
        "was not supplied for this run, so neither category was evaluated "
        "for any candidate in this report."
    ),
    # Japanese-only rendering target for the conserved-query note; English
    # keeps the analysis module's own message verbatim (see conserved_query_note).
    "conserved_query.note": (
        "Query {query_id} itself has a {strength} negative-reference BLAST hit "
        "(negative_hit_strength={strength_raw})."
    ),
}

_JA: dict[str, str] = {
    # -- title page --------------------------------------------------------
    "report.title": "ProteinHunter 候補レポート",
    "report.generated": "レポート生成日時: {timestamp}",
    "report.scoring_model": "スコアリングモデル: {scoring_model}",
    "report.code_version": "コードのバージョン: {version}",
    "report.config_fingerprint": "設定フィンガープリント: {fingerprint}",
    "report.sidecar": "実効設定の全体は、このレポートと同じ場所に {location} として保存されている",
    "report.sidecar_unknown_location": "<Excelブックの名前>.run_provenance.yaml",
    "report.reference_genome_check": "参照ゲノムのチェック: {text}",
    "report.reference_genome_warning": "警告 — 詳細は実行ログを参照(core/reference_genomes.py)",
    "report.queries_evaluated": (
        "評価したクエリ数: {n_queries}。クエリごとの表示候補数: 最大 {max_per_query} 件。"
        "加えて、順位にかかわらず Tier1_VeryStrong / Tier2_Strong の候補はすべて表示する。"
    ),
    "report.toc_placeholder": "右クリックして「フィールド更新」を選択すると、目次が生成される。",
    "sentence_joiner": "",
    # -- section 5: evidence architecture ---------------------------------
    "arch.title": "5. 証拠の構成",
    "arch.intro": (
        "この実行では {scoring_model} スコアリングモデルを使用した。以下の7つの証拠カテゴリは、"
        "パイプラインが扱う証拠の全体である。個々の候補について実際にどれが証拠を提供したかは、"
        "この実行で利用できたデータと設定に依存する(セクション 8 の各候補の"
        "「この候補が上位にある理由」の記述を参照)。"
    ),
    "arch.seq.heading": "5.1 配列の証拠 (Sequence Evidence, 上限 {cap:.0f})",
    "arch.seq.body": (
        "BLAST に基づく陽性・陰性の分類と、ベストヒットの同一性・カバレッジ・E-value の強さ。"
        "BLAST ヒットが少しでもある候補のほぼすべてで値が得られ、パイプラインの中で最も安定して"
        "利用できる証拠カテゴリである。"
    ),
    "arch.func.heading": "5.2 機能・ドメインの証拠 (Functional/Domain Evidence, 上限 {cap:.0f})",
    "arch.func.body": (
        "クエリと候補の間で共通または相補的な機能注釈(CDD/Pfam ドメイン、説明語)。"
        "ドメイン注釈が有効で、両方のタンパク質について利用できる場合に値が得られる。"
    ),
    "arch.genomic.heading": "5.3 ゲノム文脈 (Genomic Context, 上限 {cap:.0f})",
    "arch.genomic.body": (
        "GFF 座標に基づく、クエリ遺伝子と候補遺伝子のゲノム上の近接。正の証拠としてのみ用いる"
        "(遠い候補は、遠いことを理由に減点されることはなく、近いことを評価されないだけである)。"
        "両方の遺伝子について GFF の近傍データが利用できる場合に値が得られる。"
    ),
    "arch.interaction.heading": "5.4 相互作用の証拠 (Interaction Evidence, 上限 {cap:.0f})",
    "arch.interaction.body": (
        "最大3つの独立した情報源による外部のタンパク質間相互作用の証拠。STRING PPI、"
        "GEO の転写共発現、および任意の ProteinInteractionHunter(PIH)直接相互作用ブリッジである。"
        "対応する任意のデータソース(STRING の taxon ID、GEO 共発現、または PIH の証拠バンドル)が"
        "その実行で設定されている場合に値が得られる。各ソースは独立しており、どの組み合わせが"
        "利用できてもよい。"
    ),
    "arch.evolutionary.heading": "5.5 進化的証拠 (Evolutionary Evidence, 上限 {cap:.0f})",
    "arch.evolutionary.body": (
        "候補とクエリの間の系統・進化プロファイルの整合性。すべて任意の PIH ブリッジに由来する。{closer}"
    ),
    "arch.cellular.heading": "5.6 細胞内適合性 (Cellular Compatibility, 上限 {cap:.0f})",
    "arch.cellular.body": (
        "細胞内局在・適合性の証拠。これも任意の PIH ブリッジに由来する。{closer}"
    ),
    "arch.negative.heading": "5.7 負の証拠 (Negative Evidence, 予約)",
    "arch.negative.body": (
        "候補とクエリの組み合わせと直接矛盾する証拠のために予約されたカテゴリである。"
        "例として、互換性のない細胞内局在、系統的な不整合、機能的に矛盾する注釈が挙げられる。"
        "このバージョンのパイプラインには、そのような指標は実装されていない。"
        "このカテゴリは、すべての実行で、例外なく、すべての候補について未評価として報告される。"
        "これは見落としではなく、意図的な設計判断である。以前の実装では、negative_hit_strength"
        "(このレポートの別の箇所で candidate_source の下に表示される)をこのカテゴリの代用として"
        "用いたが、実データでの検証により、それが別のもの、すなわち候補タンパク質が陰性参照ゲノム"
        "にどれだけ広く保存されているかを測っていることが分かった。矛盾のペナルティとして用いると、"
        "よく保存された真の相互作用パートナーを誤って減点してしまう。このため、この2つの指標は、"
        "便宜のために再び混ぜ合わせることはせず、このレポートでは明確に分けて示している。"
    ),
    # -- section 7: ranking ------------------------------------------------
    "ranking.title": "7. 候補の順位",
    "ranking.empty": "この実行では、クエリ別の候補は得られなかった。",
    "ranking.query_heading": "7.{index} クエリ: {query_id}",
    "ranking.col.rank": "順位",
    "ranking.col.candidate": "候補",
    "ranking.col.final_score": "最終スコア",
    "ranking.col.tier": "Tier(信頼度階層)",
    "ranking.col.source": "候補ソース",
    # -- section 8: details ------------------------------------------------
    "details.title": "8. 候補の詳細",
    "details.empty": "この実行では、クエリ別の候補は得られなかった。",
    "details.query_heading": "8.{index} クエリ: {query_id}",
    "details.why_label": "この候補が上位にある理由: ",
    "details.interpretation_label": "生物学的な解釈: ",
    "details.excel_ref": (
        "この候補の全データ: {excel_filename} のシート 02_Final_Score、"
        "query_id={query_id}、candidate_protein_id={candidate_id} を参照。"
    ),
    # -- narrative: why this candidate ranks highly ------------------------
    "tier_opening.Tier1_VeryStrong": (
        "この候補は、このクエリに対して最高の信頼度階層(Tier 1 — Very Strong(非常に強い))に達しており、"
        "最終スコアは {final_score:.1f}/100 である。{evidence_category_count} 個の独立した証拠カテゴリの証拠に"
        "支えられている。"
    ),
    "tier_opening.Tier2_Strong": (
        "この候補は、このクエリに対して2番目に高い信頼度階層(Tier 2 — Strong(強))に達しており、"
        "最終スコアは {final_score:.1f}/100 である。{evidence_category_count} 個の独立した証拠カテゴリの証拠に"
        "支えられている。"
    ),
    "tier_opening.Tier3_Moderate": (
        "この候補は、このクエリに対して中程度の信頼度階層(Tier 3 — Moderate(中程度))に達しており、"
        "最終スコアは {final_score:.1f}/100 である。裏付けとなる証拠はあるが、上位階層の候補よりも少ない"
        "独立カテゴリ({evidence_category_count} 個)に由来している。"
    ),
    "tier_opening.Tier4_Weak": (
        "この候補の最終スコア {final_score:.1f}/100 は、評価された中で最も低い信頼度階層"
        "(Tier 4 — Weak(弱))に位置づけられる。この一覧での位置は、このクエリについてスコア付けされた"
        "候補の間の相対的な順位を反映したものであり、独立した強い裏付けを意味するものではない。"
    ),
    "tier_opening.unclassified": (
        "この候補は、このクエリに対して正式な最終スコアを得ていない(パイプラインが求める最低限の証拠を"
        "満たせなかった)。利用できた部分的な証拠によって上位に入ったためここに掲載されているにすぎず、"
        "順位付けされたスコア結果としてではなく、参考情報として扱うこと。"
    ),
    "why.contributing": "寄与した証拠カテゴリ: {listing}。",
    "why.listing_separator": "、",
    "why.none": "現在の実行設定では、この候補についてゼロを超えるスコアを得た個別の証拠カテゴリはなかった。",
    "source.Candidates": (
        "このタンパク質は、厳密な陽性候補(Candidates(候補))として分類された。陽性参照配列への BLAST ヒットが"
        "あり、陰性参照へのヒットが一切ない、このパイプラインが出力する最も厳しい分類である。"
    ),
    "source.Positive_all_sources": (
        "このタンパク質は、設定されたすべての陽性参照ソースにヒットし、陰性参照へのヒットはない"
        "(Positive_all_sources(全陽性ソース一致))。"
    ),
    "source.Candidates_relaxed": (
        "このタンパク質は、緩和条件の陽性候補(Candidates_relaxed(候補(緩和条件)))として分類された。"
        "陽性参照へのヒットがあり、強くない(中程度または弱い)陰性参照へのヒットは許容される。"
    ),
    "source.No_hit": (
        "このタンパク質は、陽性・陰性のいずれの参照セットにも BLAST ヒットがなかった(No_hit(ヒットなし))。"
        "配列相同性だけでは検出されない、系統特異的または新規の候補であり、ここでの順位は BLAST の一致では"
        "なく、配列以外の証拠(ゲノム文脈、ドメイン注釈、相互作用の証拠)のみに依存している。"
    ),
    "source.Negative_unmatched": (
        "このタンパク質は陰性参照へのヒットがないが、上記の、より厳しい陽性一致の基準は満たさなかった"
        "(Negative_unmatched(陰性未ヒット))。"
    ),
    "source.Negative_hit": (
        "このタンパク質は、少なくとも1つの陰性参照 BLAST ヒットを持つ(Negative_hit(陰性ヒット)、"
        "強度: {negative_hit_strength})。そのヒットがあるにもかかわらず、上に列挙した他の証拠により、"
        "この順位付けに含まれている。順位の解釈には追加の注意が必要である。"
    ),
    "negative_hit_caveat": (
        "注: この候補は、陰性参照配列に対する {negative_hit_strength} の BLAST ヒットも持つ。これは"
        "candidate_source の分類を変えるものではないが、留意すべき点である。"
    ),
    # -- narrative: biological interpretation ------------------------------
    "interpretation.opening": (
        "このパイプラインが現時点で利用できる証拠に基づくと、{candidate_id} は、クエリ {query_id} について"
        "評価された {n_candidates} 個の候補のうち {rank} 位であった。これは、以下に述べる、パイプラインが"
        "現時点で利用できる証拠カテゴリを反映したものであり、{candidate_id} が標的酵素または相互作用"
        "パートナーであることを確定的に同定したものではない。「現在利用できる証拠によって最も強く"
        "支持される」ものとして読むべきであり、確定した結論として読むべきではない。"
    ),
    "color.genomic_context": (
        "クエリ遺伝子とのゲノム上の近接は、共通のオペロンまたは機能的なまとまりと整合的であるが、"
        "それだけでそれを立証するものではない。"
    ),
    "color.functional_domain": (
        "この候補とクエリの間で、共通または相補的な機能ドメイン注釈が見つかっており、関連する、または"
        "相補的な生化学的役割と整合的である。ただし、ドメインの類似だけでは直接の相互作用を立証できない。"
    ),
    "color.interaction": (
        "外部の相互作用証拠(タンパク質間相互作用データベースおよび/または共発現データ)が、この候補と"
        "クエリの間の機能的な関連を支持しており、これは配列の類似性やゲノム上の位置とは独立している。"
    ),
    "color.evolutionary": (
        "任意の PIH ブリッジによる進化・系統プロファイルの証拠が、比較したゲノム全体にわたる、この候補と"
        "クエリの整合性を支持している。"
    ),
    "color.cellular_compatibility": (
        "任意の PIH ブリッジによる細胞内適合性の証拠が、この候補とクエリが同じ細胞内環境に適合することを"
        "支持している。"
    ),
    "negative_evidence_closer": (
        "負の(生物学的な矛盾を示す)証拠は予約されたカテゴリであり、このバージョンのパイプラインには"
        "実装された指標がない(証拠の構成、セクション 5.7 を参照)。したがって、ここに記載がないことは、"
        "この候補について矛盾する証拠が検査され、問題なしと判断されたことを意味しない。"
    ),
    "interpretation.closer": (
        "上に挙げていないカテゴリは、この実行ではこの候補について評価されなかった。空欄のカテゴリは、"
        "評価できる証拠がなかったことを意味し、証拠が調べられて存在しないと判明したことを意味しない。"
        "{evolutionary_closer}{negative_evidence_closer}"
    ),
    "evolutionary_closer.legacy": (
        "この実行では legacy_additive スコアリングモデルを使用しており、このモデルは進化的証拠と細胞内"
        "適合性の証拠をまったく評価しない。これらのカテゴリは v2_evidence_based スコアリングモデルでのみ"
        "存在する。"
    ),
    "evolutionary_closer.pih": (
        "この実行では、提供された ProteinInteractionHunter(PIH)データバンドルを用いて、進化的証拠と"
        "細胞内適合性の証拠を評価した。いずれかが空欄であるのは、この候補に該当する証拠が見つからなかった"
        "ことを意味し、そのカテゴリが省略されたことを意味しない。"
    ),
    "evolutionary_closer.no_pih": (
        "進化的証拠と細胞内適合性の証拠には、任意の外部データバンドル(ProteinInteractionHunter/PIH)が"
        "必要だが、この実行では提供されなかったため、このレポートのどの候補についても、いずれのカテゴリも"
        "評価されていない。"
    ),
    "conserved_query.note": (
        "クエリ {query_id} 自身が、陰性参照に対して {strength} の BLAST ヒットを持つ"
        "(negative_hit_strength={strength_raw})。広く保存されたタンパク質の相互作用パートナーも保存されている"
        "可能性があり、その場合は Negative_hit(陰性ヒット)に分類されることがある。現在の "
        "interaction_scoring.candidate_sources の設定(negative_hit とその強・中・弱のサブバケツがすべて無効)では、"
        "そのバケツに入る候補は、スコアにかかわらずすべての出力シートから除外される。このクエリについて想定される"
        "パートナーが出力に見当たらない場合は、candidate_sources.negative_hit(または強度別のサブバケツの"
        "いずれか)を有効にして再実行することを検討されたい。"
    ),
}

STRINGS: dict[str, dict[str, str]] = {"en": _EN, "ja": _JA}

#: Japanese glosses for classification labels, keyed by the lower-cased English label.
#: Used only in prose, as ``English(日本語)``; tables and Excel keep the English value alone.
LABEL_JA: dict[str, str] = {
    # negative_hit_strength values
    "strong": "強",
    "medium": "中",
    "weak": "弱",
    "none": "なし",
    # tier names (the part after "Tier N —")
    "very strong": "非常に強い",
    "moderate": "中程度",
    # candidate_source values
    "candidates": "候補",
    "candidates_relaxed": "候補(緩和条件)",
    "positive_all_sources": "全陽性ソース一致",
    "negative_unmatched": "陰性未ヒット",
    "no_hit": "ヒットなし",
    "negative_hit": "陰性ヒット",
    "negative_strong_hit": "陰性ヒット(強)",
    "negative_medium_hit": "陰性ヒット(中)",
    "negative_weak_hit": "陰性ヒット(弱)",
    # evidence category display labels (word_report.category_refs_for_scoring_model)
    "sequence/source classification": "配列・ソース分類",
    "genomic context": "ゲノム文脈",
    "functional/domain": "機能・ドメイン",
    "interaction": "相互作用",
    "evolutionary": "進化",
    "cellular compatibility": "細胞内適合性",
    "co-occurrence": "共起",
    "domain complementarity": "ドメイン相補性",
    "interaction (string ppi)": "STRING PPI による相互作用",
}


def normalize_language(language: Any) -> str:
    """``language`` as a supported code; None/"" mean the default. Raises ValueError otherwise."""
    if language is None or language == "":
        return DEFAULT_LANGUAGE
    if language not in SUPPORTED_LANGUAGES:
        raise ValueError(f"unsupported report language {language!r}; expected one of {SUPPORTED_LANGUAGES}")
    return str(language)


def t(key: str, language: str = DEFAULT_LANGUAGE, **values: Any) -> str:
    """Wording for ``key`` in ``language``, ``str.format``-ed with ``values``.

    A missing key raises KeyError on purpose: a typo must fail a test, not
    silently print a key name into a report.
    """
    template = STRINGS[normalize_language(language)][key]
    return template.format(**values) if values else template


def bilingual(label: str, language: str = DEFAULT_LANGUAGE) -> str:
    """``label`` for English; ``label(日本語)`` in Japanese when a gloss exists, else ``label``."""
    if normalize_language(language) == "ja":
        gloss = LABEL_JA.get(label.strip().lower())
        if gloss is not None:
            return f"{label}({gloss})"
    return label


__all__: tuple[str, ...] = (
    "DEFAULT_LANGUAGE",
    "LABEL_JA",
    "STRINGS",
    "SUPPORTED_LANGUAGES",
    "bilingual",
    "normalize_language",
    "t",
)
