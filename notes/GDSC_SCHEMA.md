# GDSC Data Schema (discovered 2026-06-04)

## Site migration note (important)
`www.cancerrxgene.org` is **decommissioned** — every path returns **HTTP 410 Gone**, including the
old `bulk_download` page and the GDSC1000 `suppData/TableS2*.xlsx` supplementary files
(`TableS2J.xlsx` methylation etc. — all 410). The site now serves a single migration banner that
redirects GDSC data to **Cell Model Passports** (Sanger DepMap):
- Downloads hub: https://cellmodelpassports.sanger.ac.uk/downloads (JS-rendered; links resolved via headless browser)
- DepMap docs: https://depmap.sanger.ac.uk/documentation/datasets/

All three datasets below were obtained from the Cell Model Passports download hub.

## Common join key
The **canonical join key across all three files is `SANGER_MODEL_ID` (SIDM, e.g. `SIDM00001`)**, NOT COSMIC_ID.
- The drug-response file has **no COSMIC_ID column at all** — only `SANGER_MODEL_ID`.
- The methylation matrix is indexed by SIDM.
- `model_list` provides the SIDM -> COSMIC_ID crosswalk (`model_id` -> `COSMIC_ID`).

The loaders honor the task contract (COSMIC_ID-keyed) by mapping SIDM -> COSMIC_ID via `model_list`
and dropping rows whose SIDM has no COSMIC_ID. COSMIC_ID coverage:
- methylation: 1055 / 1523 SIDMs have a COSMIC_ID
- drug: 969 / 969 SIDMs have a COSMIC_ID
- all 969 drug SIDMs are present in the methylation matrix

---

## 1. Methylation beta matrix
- **URL:** https://figshare.com/ndownloader/files/43146652  (the "Methylomics data" link on the CMP
  downloads page; this is the DepMap/CMP multi-omics Figshare+ bundle).
  - NOTE: `figshare.com/ndownloader/...` returns **HTTP 202 forever** (async prep that never completes).
    The working host is **`https://ndownloader.figshare.com/files/43146652`**, which 302-redirects to a
    short-lived (10 s) presigned S3 URL — `curl -L` in one shot captures it.
- **Local file:** `data/gdsc/methylation_imputed.csv.gz` (61 MB gz, ~175 MB raw)
- **Original name (from gzip header):** `20231023_092657_imputed_methylation.csv`
- **Format:** plain CSV, gzip-compressed. First column is unnamed (the SIDM row id); remaining columns
  are feature columns. Header row 0; data rows follow.
- **Shape:** 1523 rows (cell lines, SIDM) x 14608 feature columns.
- **Values:** beta in [0,1] (imputed). Verified >99% of values fall in [0,1].
- **Feature id -> gene mapping:** **features ARE gene symbols directly.** Columns are HGNC gene symbols
  (`A1BG`, `A2BP1`, `ABCB1`, ...). This is a **gene-level aggregated/imputed methylation** matrix, NOT
  raw CpG-island / 450K-probe IDs. So no separate annotation table is needed: column name == gene symbol.
  - This DIFFERS from the old Iorio-2016 GDSC1000 expectation (~14,725 CpG *islands*). Here the ~14.6k
    features are *genes*, not islands. Feature count is similar by coincidence; semantics differ.
- **MGMT feature (downstream dependency):** **`MGMT`** — exact column name is the bare gene symbol
  `MGMT` (column position 7625 in the raw header). Present and unambiguous. A downstream task can read
  the MGMT promoter-methylation surrogate as `methylation["MGMT"]`.

## 2. Drug response (GDSC2 fitted dose-response)
- **URL:** https://cmp.cog.sanger.ac.uk/download/GDSC2_fitted_dose_response_27Oct23.xlsx
- **Local file:** `data/gdsc/GDSC2_fitted_dose_response_27Oct23.xlsx`
- **Format:** Excel, single sheet `Sheet1`. Long format.
- **Shape:** 242036 rows x 16 columns. 969 cell lines, 286 drugs.
- **Columns (real):** `DATASET, NLME_RESULT_ID, NLME_CURVE_ID, CELL_LINE_NAME, SANGER_MODEL_ID,
  CANCER_TYPE, DRUG_ID, DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME, MIN_CONC, MAX_CONC, LN_IC50, AUC,
  RMSE, Z_SCORE`
- **Loader renames:** `SANGER_MODEL_ID`->(mapped to)`COSMIC_ID`, `DRUG_NAME`->`drug_name`,
  `LN_IC50`->`ln_ic50`, `AUC`->`auc`.
- **Temozolomide spelling:** `Temozolomide` (capital T, single token). Loader lowercases `drug_name`
  comparisons are case-insensitive-safe but the stored value is kept as-is (`Temozolomide`).

## 3. Cell-line annotations (model list)
- **URL:** https://cog.sanger.ac.uk/cmp/download/model_list_20260420.csv
- **Local file:** `data/gdsc/model_list_20260420.csv`
- **Format:** CSV.
- **Shape:** 2266 rows x 98 columns.
- **Key columns used:** `model_id` (SIDM), `COSMIC_ID`, `model_name`, `tissue`, `cancer_type`.
  - `COSMIC_ID` non-null for 1124 / 2266 models.
- **Loader output:** index = COSMIC_ID; columns include `cell_line_name` (from `model_name`),
  `tissue`, plus `sanger_model_id`, `cancer_type`.

---

## Status flags
- Methylation features are **gene-level betas**, not CpG islands. The MGMT downstream dependency is
  satisfied (column `MGMT`). Beta range [0,1] holds. ~14.6k features (> the 1000 the test needs).
