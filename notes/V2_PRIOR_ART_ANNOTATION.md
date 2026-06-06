# v2 research pass — R2: prior art, annotation feasibility, positive controls (2026-06-05)

Deep-research pass (adversarially verified, cited). Full transcript in the task log; this is the synthesis.

## A. Prior art / novelty

**Closest prior art (cite + beat):**
- **Scompairi et al., Comms Biol 2023** ("pharmacoepigenomic landscape of cancer cell lines") — 721 GDSC lines, 22 cancer types, methylation→drug-response biomarkers, **mechanism-anchored via PPI to drug targets**, 19 biomarkers/17 drugs. **BUT cancer-type-STRATIFIED ("analysis performed for each cancer type separately"), one-shot paper, NOT cross-indication, NOT a re-runnable engine.** The single most important paper to differentiate from.
- **CellMinerCDB** (NAR 2021/2026) — the closest *portal*; integrates methylation as a first-class feature, methylation⇄drug correlation + LASSO — but **per-query, interactive, one gene/one drug at a time; no systematic cross-indication ranking.**
- CMap/LINCS = expression-signature repurposing (not methylation). DepMap/PharmacoDB carry methylation as data, not as a cross-indication opportunity product.
- **PMC7906305** ("between-tissue differences on pan-cancer drug-sensitivity predictions") = **published corroboration of our R1b finding** that naive pan-tissue ranking is confounded → cite as justification for mechanism-anchoring.

**NOT novel:** methylation as a drug-response feature; GDSC methylation→drug biomarker discovery; the tissue-agnostic paradigm (all FDA agnostic markers are mutation/fusion/MSI/TMB — none methylation); computational repurposing (crowded); **the CDKN2A-methylation→palbociclib link itself is already published in GDSC (PLOS One 2019, n=522) → our positive control is CONFIRMATION, not discovery.**

**Defensible novelty (the only claim that survives):** the *combination/framing* — (1) methylation-defined state as the **primary unit, ranked ACROSS the indication taxonomy** (vs cancer-type-stratified Comms Bio 2023 / per-query CellMinerCDB); (2) **mechanism-anchored** to the drug's own target/pathway to defeat the co-methylation confound; (3) **provenance-stamped**, re-runnable engine. Combination novelty, NOT modality or biology.

**Sharpest honest pitch sentence:**
> "Existing pharmaco-epigenomic work ranks methylation–drug biomarkers *one cancer type at a time* (Comms Bio 2023) or lets you explore correlations *one query at a time* (CellMinerCDB); Strata treats a methylation-defined state as the primary unit and ranks it *across* the indication taxonomy — anchoring each scan to the drug's own known target/pathway to defeat the between-tissue co-methylation confound, and provenance-stamping every opportunity."

## B. Annotation feasibility — SOLVED for Saturday

**GDSC2's own columns suffice, zero external dependency, no license risk** (verified on the local xlsx):
- `PUTATIVE_TARGET`: populated **249/286 drugs (87%)**, 185 distinct targets (mix of gene symbols like `EGFR`,`PARP1, PARP2`,`MEK1, MEK2` and mechanism classes like `Antimetabolite`).
- `PATHWAY_NAME`: **286/286 (100%)**, 24 coarse categories (`Cell cycle`, `PI3K/MTOR signaling`, `DNA replication`, `Chromatin histone acetylation`, …).
→ **Mechanism anchor = for a drug, test methylation only at its target gene(s) + `PATHWAY_NAME` genes.** Build it from the local file.

**Optional, license-clean enrichment:** DGIdb (MIT, TSV) to map the 37 gaps + mechanism-class drugs to gene symbols; MSigDB Hallmark (50 GMT, gene-symbol) and **Reactome (CC0/CC-BY)** for pathway expansion. **Avoid for a commercial build:** full DrugBank (CC BY-NC), raw KEGG (paid commercial), MSigDB KEGG_LEGACY subset. *(MSigDB exact EULA clause not fully verified — confirm before commercial redistribution.)*

## C. Positive / negative controls (ranked)
1. **CDKN2A/P16-methylation → palbociclib** — cleanest, fully methylation-recoverable in GDSC; **prior-known (PLOS One 2019) → confirmation, not discovery.** The spine.
2. **MTAP-loss → PRMT5/MAT2A inhibitors** (TNG908/462, AMG-193, clinically hot, tissue-agnostic) — strongest *cross-indication* narrative, **but CNV-led** (methylation rides the 9p21 co-deletion) → partial control if methylation-only.
3. **MLH1-methylation → MSI → immunotherapy** — molecular-state control, but the drug readout (checkpoint blockade) is **not recoverable in cell-line viability**; narrative anchor only.
4. **MGMT → temozolomide** — null at gene-level (multifactorial, protein/expression-level signal, glioma-concentrated) → **use as the honest NEGATIVE control** (engine isn't spuriously firing).

## Bottom lines
- Mechanism-anchored discovery is **buildable this week from GDSC's own annotations.**
- The pitch must claim **combination/framing novelty**, cite Comms Bio 2023 + CellMinerCDB, and call the positive control **confirmation**.
- Add **MGMT as a negative control** and **MTAP→PRMT5 as the cross-indication opportunity narrative** (with the CNV caveat).
