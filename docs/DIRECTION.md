# Strata — Current Direction (2026-06-06)

## Where we are
- **v1** (shipped, deployed): single-drug, within-indication methylation→drug-response subgroup agent. `methylkb.vercel.app/strata.html`.
- **v2** (built, deployed, submitted to `Lapintam/Biohack_626`): cross-indication opportunity engine — mechanism-anchored, tissue-agnostic, evidence-graded. `methylkb.vercel.app/strata-v2.html`.
- **Data:** joined master table (`data/master/`) + full gene-level (1,055×14,608) and promoter (614×17,181) methylation matrices, all keyed by `COSMIC_ID`; drug response (286 drugs), genomics (TMB/MSI/ploidy/mutation flags).

## What we've learned (the honest spine)
1. **Gene-level imputed methylation inflates signal** — it encodes copy-number at deletion loci (9p21/MTAP-CDKN2A), so the original positive control was substantially a CNV artifact, not promoter methylation.
2. **Naive unbiased all-gene discovery is non-specific** (R1b) — the methylome is pervasively co-correlated (~30–60% of genes "significant").
3. **Mechanism-anchoring** gives interpretable, specific hits but is **bounded to known biology**.
4. **Promoter resolution + global-axis (latent-factor, SVA/RUV-style) control SOLVES the specificity problem** — unbiased genome-wide discovery collapses to a tractable ~few-hundred-gene candidate list (~2–5% significant), *and* removes the CNV confound. Canonical mutation/CNV-driven genes correctly don't appear (their role isn't promoter methylation).

## Two tracks
- **Track A — the defensible product (demo-ready):** the mechanism-anchored, tissue-agnostic, evidence-graded engine + honesty rails. Credibility = the framework that caught its own CNV artifact. This is the pitch.
- **Track B — the research frontier (GOING FOR IT):** promoter-resolution, global-axis-controlled **unbiased discovery** (now tractable) → **multi-omic validation** (RNA-seq = functional silencing; CNV = exclude deletion) → **external datasets** (PRISM/CTRP) → **patient cohorts + compartment attribution** (clinical-translation north star).

## The pipeline (Track B)
```
promoter methylation (CCLE RRBS, true regulatory signal, no CNV proxy)
  → global-axis-controlled all-gene scan  (tissue + methylome PCs removed)
  → ~few-hundred candidate genes / drug    [specific, but unvalidated]
  → RNA-seq silencing filter                (promoter-meth ↑ ⇒ expression ↓ = functional)
  → handful of VALIDATED novel leads
  → external validation (PRISM/CTRP/literature)
  → (later) patient cohorts + compartment attribution
```

## Immediate next steps
1. **RNA-seq layer** (building now): CCLE expression, same lines/crosswalk → for each promoter candidate, test whether its methylation tracks its own expression *down* (functional silencing). Turns the candidate list into validated leads.
2. **External validation** (#2, queued): re-test surviving leads in PRISM/CTRP.

## The meta-point
The product's value is the **rigorous, honest discovery framework** — tissue-confound safeguards, the evidence ladder (L0–L6), and the CNV / promoter / co-correlation controls — which repeatedly catch confounds that naive methylation-drug-response work reports as findings. Track B is that framework maturing into genuine, validated discovery.
