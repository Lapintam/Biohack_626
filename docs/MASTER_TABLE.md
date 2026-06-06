# Master table — `data/master/` (gitignored, regenerable)

`scripts/build_master_csv.py` → joins everything keyed by `COSMIC_ID`.

## `master_cellline_drug.csv` (242,036 rows × 140 cols)
**Grain: one row per (cell line × drug).** Column groups:
- **Identity:** cell_line, COSMIC_ID, sanger_model_id, tissue, cancer_type
- **Drug + sensitivity:** drug, drug_target, drug_pathway, auc (lower = sensitive), ln_ic50
- **Genomics:** TMB (=mutational_burden), MSI (MSS/MSI), ploidy (coalesced WGS→WES→SNP6), mmr_status, mlh1_promoter_meth, gender, age, braf/kras/pik3ca/pten_mut
- **Methylation panel (119 cols):** `meth_gene__<GENE>` (gene-level imputed) + `meth_prom__<GENE>` (CCLE RRBS promoter) for ~60 curated cancer genes.

## Full methylation matrices (separate, joinable by COSMIC_ID)
- `methylation_genelevel_bycosmic.csv.gz` — 1,055 × 14,608 (all genes, gene-level)
- `methylation_promoter_bycosmic.csv.gz`  — 614 × 17,181 (all genes, promoter)

## Notes
- Methylation: gene-level is **imputed** (can encode copy-number at deletion loci, e.g. 9p21 — see `notes/V2_PROMOTER_VALIDATION.md`); promoter is **RRBS** (NaN where uncovered/deleted).
- Regenerate: `PYTHONPATH=. uv run python scripts/build_master_csv.py`
- Cell-line × drug × methylation × genomics all share `COSMIC_ID` (100% join across the 3 GDSC files; 557 lines have promoter methylation).
