# Patient-methylation annotation — engine + held-out validation (2026-06-06)

`strata/annotate.py` (scorer) + `scripts/build_annotation_assets.py` (model) + `scripts/validate_annotation.py` (validation).

## How it works
Given a patient promoter-methylation profile {gene: beta}, for each validated lead: percentile the
patient's beta within the cell-line cohort, flip by the lead's direction (meth->sensitive vs ->resistant),
weight by |gdsc_r|*(screens), aggregate per drug -> sensitive / resistant / neutral call + driver genes.
Stateless: needs only `api/reports/annotation_model.json` (lead defs + per-gene cohort distributions,
24 drugs / 212 scorable genes / ~0.9 MB) -> runs inside the API.

## Held-out validation (honest, no leakage)
Split lines train/test, rediscover leads on TRAIN only, build score from train leads + train cohort,
predict on TEST, correlate predicted sensitivity vs MEASURED AUC:
- **Olaparib**: held-out Spearman rho = **-0.251** (p=0.001, n=170)
- **Trametinib**: held-out rho = **-0.373** (p=5e-7, n=170); top-quintile predicted-sensitive median AUC 0.864 vs bottom 0.975
Negative rho = works (higher predicted sensitivity -> lower AUC). Modest but real on unseen lines.

## Honest framing
Research-use; cell-line-derived; weak single-gene signals aggregated -> enrichment/triage, not deterministic.
In-cohort cell lines carry training leakage (illustrative); a real patient tumor/PBMC sample is leakage-free (the honest demo input).

## Sanity check
Scorer surfaces SLFN11 as a Gemcitabine driver — a known biomarker of DNA-damage-agent sensitivity.

## Next
- /annotate endpoint (POST profile -> predictions) on the FastAPI; redeploy Fly (model ships in image).
- Live LLM agent endpoint (full loop): NL -> agent -> report API + /annotate -> narrated annotated report.
- Front-end patient section. Real TCGA tumor / PBMC demo sample (leakage-free).
