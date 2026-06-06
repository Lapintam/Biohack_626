# v2 crux probe — tissue-confound result (R1)

Date: 2026-06-05 · script: `scripts/v2_tissue_confound_probe.py` · branch: `feat/v2-cross-indication`

**Question (the v2 load-bearing assumption):** does a marker→drug signal survive the §8a tissue-confound safeguards — i.e. is it *tissue-agnostic* (holds within tissues + after tissue adjustment), or just pooled-significant lineage?

**Substrate:** GDSC, 30 tissues, ~960 lines for palbociclib. Within-tissue tested at n≥20 (14 tissues); tissue-adjusted = OLS `auc ~ z(meth) + C(tissue)` (rare tissues pooled to "Other").

## Result

| Pair | Pooled ρ | Within-tissue (n≥20) | Between-tissue ρ | Tissue-adj effect retained | Verdict |
|---|---|---|---|---|---|
| **MTAP → Palbociclib** | −0.440 | 13/14 neg, 5 sig-neg (lung, breast, colorectal, H&N, ovary) | −0.626 | **60%** (p=1.4e−9) | **TISSUE-AGNOSTIC → L2/L3** |
| **CDKN2A → Palbociclib** | −0.336 | 12/14 neg, 4 sig-neg | −0.464 | **74%** (p=5.7e−9) | **TISSUE-AGNOSTIC → L2/L3** |
| **DAPK1 → Luminespib** | −0.346 | 12/14 neg, 4 sig-neg | −0.635 | **34%** (p=0.007) | **PARTIAL — mostly lineage** |

## Interpretation

1. **v2 thesis validated for the positive control.** 9p21/MTAP-CDKN2A → palbociclib is genuinely tissue-agnostic: significant *within* multiple unrelated tissues and 60–74% of the effect survives tissue adjustment. This is "molecular state beats tissue taxonomy," recovered cold and confound-controlled (evidence level **L3**).
2. **The safeguards discriminate.** DAPK1→luminespib (a v1 "novel hit") is largely between-tissue (only 34% retained) → correctly demoted to PARTIAL/lineage. The ladder has teeth; the §8a safeguards are not decorative.
3. **The real v2 engine = rank by *tissue-adjusted* association**, not pooled. Pooled `rank_response_markers` will surface lineage hits; ranking by the tissue-adjusted coefficient (or within-tissue consistency) surfaces genuinely taxonomy-crossing states. Next step (R1b): a tissue-agnostic marker scan to find a *novel* L2/L3 cross-indication hit (DAPK1 doesn't qualify cleanly).

## R1b — tissue-adjusted marker scan (the engine attempt)

Script: `scripts/v2_tissue_agnostic_scan.py`. Ranks genes by **within-tissue-residualized** correlation with AUC (fixed-effects on tissue) = the §8a-safeguarded discovery rule.

**Validate (Palbociclib):** MTAP rank **106/14608** (r_adj=−0.186, q=2.4e−6), CDKN2A rank **170** (q=6.1e−6) — both **survive tissue adjustment** as top-1% significant. The positive control is genuinely tissue-agnostic. ✓

**But the discovery side does NOT work naively:**
- **~7,000–8,000 of 14,608 genes are FDR-significant** after adjustment, and the **top single markers per drug are obscure / non-mechanistic** (FAM167B→Venetoclax, DYSFIP1→Sorafenib, KCNAB3→Talazoparib, STXBP3→Olaparib, …).
- Signature of **methylation's pervasive inter-gene co-correlation + a dominant within-tissue axis** (global methylation / proliferation / cell-state). Single-gene tissue-adjusted ranking is **non-specific** — it can't separate specific drug biology from the background axis.
- ⇒ "rank by tissue-adjusted correlation, take top marker" would surface obscure genes riding a global axis, then mechanism-launder them (audit pitfall #8). **Rejected as the discovery rule.**

**Strategic implication (converges with the audit):**
- The defensible v2 discovery is **mechanism-anchored, not top-of-scan**: test whether a drug's *known target / pathway* genes carry a tissue-agnostic methylation signal (marker mechanistically connected by construction). Needs drug-target/pathway annotations (R2 / §9.3).
- Alternatively/additionally: **control for the global within-tissue axis** (e.g., regress out PC1 / mean methylation / a proliferation proxy) before ranking.
- The **positive control (MTAP/CDKN2A→palbociclib) survives and remains the demo spine**; novel hits require the mechanism anchor or global-axis control — not brute-force ranking.

## R3a — mechanism-anchored engine (the corrected, demoable discovery)

Script: `scripts/v2_mechanism_anchored.py` → `outputs/v2_mechanism_anchored.csv`. For each drug, restrict the search to its **own mechanism genes** (GDSC `PUTATIVE_TARGET` + a curated `PATHWAY_NAME`→genes map) and find the one with a genuinely **tissue-agnostic** (within-tissue + adjusted) methylation→AUC signal, evidence-graded L0–L3.

**Results (level counts: L3=16, L2=83, L1=36, L0=64):**

| Role | Drug(s) | Marker (pathway) | Level | r_adj / r_pool | within-tissue sig | Note |
|---|---|---|---|---|---|---|
| **Positive control** | Palbociclib (+ Ribociclib L2) | **MTAP / 9p21** (Cell cycle) | **L3** | −0.196 / −0.461 | 5 tissues | confirmation; **replicated across 2 CDK4/6i** |
| **Novel lead A** | Rapamycin, Dactolisib, AZD8055 | **INPP4B** (PI3K/MTOR) | **L3** ×3 | −0.14…−0.21 | 3 tissues each | **replicated across 3 mTOR-pathway drugs**; PI3K-pathway TS silencing → mTORi dependence |
| **Novel lead B** | SB216763 (+ CHIR-99021) | **SFRP1** (WNT) | **L3** | −0.169 / −0.183 | **6 tissues** | classic pan-cancer WNT-antagonist methylation; most cross-tissue |
| Other L3 | Nutlin-3a→CDKN2A (p53), Entinostat→HDAC1, Nelarabine→RRM1, … | | L3 | | ≥3 | mechanism-anchored hypotheses |
| **Negative control** | Temozolomide | MGMT | **null** | +0.078 (wrong dir) | 0 | engine not spuriously firing |

**Why this is the right engine:** every hit is interpretable by construction (marker ∈ the drug's target/pathway), tissue-agnostic by the §8a safeguards, and **internally replicated** (CDK4/6i→9p21 ×2; mTOR-pathway→INPP4B ×3) — replication across a drug class is real-biology evidence, not noise. The positive control is **confirmation** (CDKN2A→palbo published); the novel leads are **L3 hypotheses**. The ranked table = the portfolio opportunity map.

**Honest caveats specific to R3a:** best-of-mechanism-genes per drug is a mild selection (mitigated by the within-≥3-tissue consistency bar); the curated pathway map is hand-built (auditable, license-clean); INPP4B's TS-vs-oncogene role is context-dependent — present as hypothesis.

## Honest caveats
- There IS a real between-tissue (lineage) component (MTAP between-tissue ρ=−0.626); we claim only the surviving within-tissue core, not the pooled effect.
- A 0/1 median-split on a tissue-agnostic marker is still post-selection; the discovery statistic remains the tissue-adjusted continuous association.
- Cell-line, tumor-compartment, gene-level methylation — all prior rails still hold.
