# TCGA-LAML Data Schema (discovered 2026-06-04)

Real-patient anchor for the methylation-subgroup vignette: cluster acute myeloid
leukemia (AML) patients on DNA methylation, test for differential overall survival.

## Source — UCSC Xena, GDC TCGA Acute Myeloid Leukemia (LAML) hub
Cohort: "GDC TCGA Acute Myeloid Leukemia (LAML)". Browser:
https://xenabrowser.net (datapages are JS-rendered). The GDC hub serves flat
matrices from a public S3 bucket with stable, predictable paths
(`https://gdc-hub.s3.us-east-1.amazonaws.com/download/<dataset>`):

| Dataset | URL | Local file | Size (gz) |
|---|---|---|---|
| Methylation (Illumina HumanMethylation450 beta) | https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LAML.methylation450.tsv.gz | `data/tcga_laml/TCGA-LAML.methylation450.tsv.gz` | 663,582,724 B (~663 MB) |
| Survival (curated OS) | https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LAML.survival.tsv.gz | `data/tcga_laml/TCGA-LAML.survival.tsv.gz` | 1,649 B |

Both downloaded with a single `curl -sL` (the S3 bucket serves directly, no async
prep / presign dance). `data/` is gitignored.

## Sample-id convention (the join key)
Both matrices key on the **TCGA aliquot barcode at the sample level**, e.g.
`TCGA-AB-2820-03A` (16 chars). The methylation column headers and the survival
`sample` column use the **same** barcodes, so they join directly — **no trimming
or normalization needed**. Verified exact-string intersection.

- Methylation column ids: `TCGA-AB-2820-03A`, `TCGA-AB-2834-03A`, ... (all 16-char, `-03A` = primary blood-derived cancer aliquot A).
- Survival `sample` ids: same form. The survival table also carries `_PATIENT`
  (12-char patient barcode, e.g. `TCGA-AB-2840`) and contains **multiple aliquots
  per patient** (`-03A`, `-03B`, `-03D` with identical OS). The clinical loader
  keeps the first occurrence per barcode; the join then falls out on the `-03A`
  barcodes shared with the methylation matrix.

## 1. Methylation beta matrix
- **Format:** tab-separated, gzip-compressed. First column `Composite Element REF`
  = 450K probe id (cg/rs/ch...). Header row = sample barcodes.
- **Native shape:** 486,427 probes (rows) x 194 samples (columns).
- **Orientation:** probes x samples. The loader **transposes to samples x probes**
  to match the engine convention (samples as ROWS).
- **Values:** beta in [0, 1]. (>99% of finite values in range; test asserts this.)
- **Missing-value handling (documented contract):** the 450K matrix has many NaNs
  (probes not measured in a given sample). The loader **drops probes missing in
  >20% of samples**, then **fills remaining NaN with the per-probe mean**.
  - After filtering: **417,326 probes** retained, **0 NaN** remaining.
  - PERFORMANCE NOTE: the filter + fill are done in the *native probes-x-samples*
    orientation (probes as rows), vectorized with numpy, and the transpose is done
    LAST. Doing the column-wise `isna().mean()` / `fillna(mean())` *after* the
    transpose (over ~417K columns) is pathologically slow in pandas (>10 min,
    killed); the row-wise/numpy path is ~0.6 s. Result is identical.

## 2. Survival (curated OS)
- **Format:** tab-separated, gzip-compressed.
- **Shape:** 250 rows x 4 columns.
- **Columns:** `sample`, `OS.time` (overall-survival time, **days**, float),
  `OS` (event: 1 = dead, 0 = censored/alive, int), `_PATIENT`.
- **Loader output (`load_tcga_clinical`):** index = `sample` barcode (str);
  columns `os_time` (float days) <- `OS.time`, `os_event` (int 0/1) <- `OS`.
  Rows with missing OS/OS.time, non-0/1 event, or negative time are dropped;
  duplicate barcodes deduped (keep first).
  - BUG AVOIDED: building the output DataFrame must use `.values` for the columns
    and a fresh `pd.Index` for `sample`, or pandas aligns the source RangeIndex
    against the new string index and silently yields all-NaN columns.

## Vignette result (`scripts/tcga_vignette.py`, 2026-06-04)
- **n patients (methylation x survival):** 169
- **n clusters (auto-k via silhouette on top 2000 variable probes):** 2 (sizes 85 / 84)
- **log-rank p = 0.844 — NOT significant.**
- **Honest read:** the *unsupervised top-variance* 2-way split is survival-neutral
  here. AML methylation subtypes *do* often differ in survival, but those classes
  track specific cytogenetic / mutational groups (e.g. CEBPA, IDH, RUNX1-RUNX1T1),
  not the global highest-variance partition that silhouette-driven k-means lands
  on. The vignette is still a real-patient methylation-subgroup demonstration; the
  null survival split is reported as-is, not massaged. A supervised or
  feature-targeted clustering would be the next step to recover prognostic
  subgroups.
- **Figures:** `outputs/tcga_laml_km.png`, `outputs/tcga_laml_umap.png`.
