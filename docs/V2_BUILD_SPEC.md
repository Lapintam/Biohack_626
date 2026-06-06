# Strata v2 — Cross-Indication Opportunity Engine · AS-BUILT SPEC

Date: 2026-06-06. Branch: `feat/v2-cross-indication`. Status: **Saturday slice BUILT** (research pass R1/R1b/R2/R3a complete; demo live additively). This is the concrete spec the external audit asked for — documenting what exists, the exact method, and what's roadmap. Vision/strategy lives in `OPPORTUNITY_GRAPH_PLAN.md`; provenance in `notes/V2_TISSUE_CONFOUND.md` + `notes/V2_PRIOR_ART_ANNOTATION.md`.

## 1. What it is
Start from **molecular state, not disease taxonomy.** For each drug, find a methylation-defined state whose drug-response association **holds across the indication taxonomy** (tissue-agnostic), anchored to the drug's own mechanism, evidence-graded, ranked into an opportunity table. v1 (`strata.html`) = one drug, within-indication. v2 (`strata-v2.html`) = across indications.

## 2. The method (exact)
1. **Tissue-confound safeguards** (`scripts/v2_tissue_confound_probe.py`): a marker→drug signal is *tissue-agnostic* only if it survives **within-tissue** (Spearman per tissue, n≥20) **and tissue adjustment** (within-tissue residualized correlation `r_adj`). Pooled significance alone = possibly lineage (proven: the unsupervised null + R1b).
2. **Mechanism anchor** (`scripts/v2_mechanism_anchored.py`): for each drug, restrict the search to its **own mechanism genes** = GDSC `PUTATIVE_TARGET` gene(s) ∪ a curated `PATHWAY_NAME`→genes map (`PATHWAY_GENES` in the script; tumor-suppressor/regulator-emphasized, license-clean). Pick the gene with the strongest tissue-agnostic, sensitivity-direction (`r_adj < 0`) signal. *Why:* naive all-gene tissue-adjusted ranking is non-specific (methylation co-correlation; R1b).
3. **Evidence ladder (L0–L6):** L0 pooled · L1 interpretable/known-biology · **L2 survives tissue adjustment** · **L3 tissue-agnostic (consistent within ≥3 cancer types)** · L4 mechanism path · L5 external support · L6 clinical. Grading: L3 = `r_adj<0`, `p_adj<0.05`, `|r_adj|≥0.12`, within-tissue-significant in ≥3 tissues.

## 3. Results (the demo content)
- **Positive control (CONFIRMATION):** 9p21/**MTAP** → **Palbociclib**, L3 (r_adj=−0.196, sig in 5 tissues); **replicated** by Ribociclib (L2). The CDKN2A→palbociclib link is published (PLOS One 2019) — recovered cold + confound-controlled, *not a discovery claim*.
- **Novel leads (HYPOTHESES, L3):** **INPP4B → mTOR/PI3K inhibitors** (replicated across Rapamycin/Dactolisib/AZD8055); **SFRP1 → WNT modulators** (most cross-tissue, 6 tissues).
- **Negative control:** **MGMT → Temozolomide** = null (r_adj +0.08, 0 tissues) — engine doesn't force canonical results.
- **Opportunity table:** top-20 ranked L3/L2 mechanism-anchored candidates = the portfolio-triage view.

## 4. Artifacts & UI
- **Engine:** `scripts/v2_mechanism_anchored.py` → `outputs/v2_mechanism_anchored.csv`.
- **Demo data:** `scripts/generate_v2_demo.py` → `site/demo-v2/opportunities.json` (featured cards w/ real per-tissue recurrence panels, opportunity table, evidence ladder, scripted agent transcript).
- **Page:** `site/strata-v2.html` + `strata-v2.js` + `strata-v2.css` (additive; v1 untouched, linked from v1 nav). §01 thesis/novelty · §02 agent transcript (recorded) · §03 featured cards (level badge + tissue-recurrence bars + mechanism card + replication chips) · §04 negative control · §05 opportunity table · §06 evidence ladder + honesty.
- **Regenerate:** `PYTHONPATH=. uv run python scripts/v2_mechanism_anchored.py && PYTHONPATH=. uv run python scripts/generate_v2_demo.py`.

## 5. Honesty rails (must hold)
- Novelty = **combination/framing only** (cross-taxonomy + mechanism-anchored + provenance); cite **Comms Biol 2023** (cancer-type-stratified) + **CellMinerCDB** (per-query). Methylation-as-feature, the biology, and CDKN2A→palbo are NOT novel.
- Positive control = **confirmation**; novel leads = **screen-derived hypotheses** (L3), licensed by the control. Mechanism cards = annotation, not proof.
- Tumor cell-line, gene-level methylation, research-use — not blood/liver/patient, not clinical decision support (`SYSTEM_MAP.md` Parts 4 + 7).
- "Graph" is roadmap; v2 ships an **opportunity map / ranked table**, not a knowledge graph.

## 6. NOT built (roadmap)
Multi-omic state, a real knowledge graph + traversal, external validation (PRISM/CTRP/CCLE), patient cohorts, autonomous (live-API) agent over v2, the federated claims universe. Pitch as horizon.
