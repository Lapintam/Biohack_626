# External (orthogonal) drug-response datasets — schema & provenance

Two independent cell-line drug-sensitivity screens acquired to **replicate**
methylation→drug-response leads first found in GDSC2. Both are keyed back to
`COSMIC_ID` via the Cell Model Passports model list
(`data/gdsc/model_list_20260420.csv`), the same crosswalk the GDSC/CCLE loaders
use. External line IDs are DepMap/Broad `ACH-…` IDs → matched to the
`BROAD_ID` column → `COSMIC_ID`.

Loaders live in `strata/data/external_drug.py`; both return the same LONG schema
as `strata.data.gdsc.load_gdsc_drug_response`:
`["COSMIC_ID", "drug_name", "auc"]` (COSMIC_ID str, drug_name lowercased/stripped).

---

## PRISM Repurposing — Secondary screen (must-have)

- **Source:** DepMap download manifest `https://depmap.org/portal/api/download/files`
  (a CSV). Release **"PRISM Secondary Repurposing 20Q2"**, filename
  `secondary-screen-dose-response-curve-parameters.csv`. The manifest `url` is a
  figshare signed link: `https://ndownloader.figshare.com/files/36794595`.
- **Local file:** `data/prism/secondary-screen-dose-response-curve-parameters.csv`
  (~290 MB; verified > 5 MB).
- **Shape:** 753,778 data rows × 20 columns.
- **ID type:** `depmap_id` (`ACH-…`) → `BROAD_ID` → `COSMIC_ID`.
  531 of 737 unique lines map to a COSMIC_ID (rest dropped). 10,137 rows have a
  null `depmap_id` (control barcodes / unresolved) and are dropped.
- **Drug-name column:** `name` (e.g. `palbociclib`, `trametinib`).
- **AUC column:** `auc` = area under the fitted **dose-response viability curve**
  (4-param log-logistic, drc R package, fit to viability where viability is
  2^(log2 fold-change vs DMSO)). Viability ≈ 1 = no effect, ≈ 0 = full kill.
  Min/median/max observed = 0.000 / 0.907 / 4.889 (some curve-fit overshoot > 1).
  **AUC direction: LOWER auc = MORE sensitive.** Same convention as GDSC2 AUC.
  **No transform applied.**
- **Screens present (`screen_id`):** HTS002 (602,495 rows), MTS010 (116,977),
  MTS006 (33,472), MTS005 (834). All retained — `drug_name` is collapsed across
  screens/doses by the loader (one mean AUC per (line, drug)). The readme notes
  MTS010 is the redo of the oncology compounds; we do not privilege it because we
  want maximal compound coverage and a single per-(line,drug) value.

## CTRP — Cancer Therapeutics Response Portal v2 (acquired)

- **Source:** DepMap download manifest, release **"Harmonized CTD^2 25Q2"**.
  - `CTRPAUCMatrix.csv` — wide AUC matrix (rows = `ModelID` `ACH-…`,
    cols = `CompoundID` `DPC-…`), GCS-hosted signed URL.
  - `CTRPResponseCurves.csv` — long curve-parameter table; used only to map
    `CompoundID` → `CompoundName`.
  - (The classic `CTRPv2.0_2015_ctd2_ExpandedDataset.zip` with the legacy
    `v20.data.curves_post_qc.txt` / `v20.meta.per_compound.txt` /
    `v20.meta.per_cell_line.txt` three-file join is also listed in the manifest,
    but its NCI URL `ctd2-data.nci.nih.gov/...` now returns an HTML redirect
    landing page, not the zip, and the figshare mirrors 403. The Harmonized
    25Q2 files are the clean, currently-downloadable equivalent and give the same
    AUC-per-(line,compound) signal, already keyed by `ACH-` like PRISM.)
- **Local files:** `data/ctrp/CTRPAUCMatrix.csv` (~6 MB),
  `data/ctrp/CTRPResponseCurves.csv` (~45 MB).
- **Shape:** AUC matrix 841 lines × 545 compounds (→ ~458k long rows after melt,
  dropping NaNs); 545/545 CompoundIDs resolve to a name.
- **ID type:** `ModelID` (`ACH-…`) → `BROAD_ID` → `COSMIC_ID`. 614 of 841 lines
  map to COSMIC.
- **Drug-name column:** `CompoundName` (joined from `CTRPResponseCurves.csv`).
- **AUC column:** Harmonized normalized `auc`, in [0, 1]
  (min/median/max = 0.000 / 0.892 / 1.000), area under the **viability** curve.
  **AUC direction: LOWER auc = MORE sensitive.** Same as GDSC/PRISM.
  **No transform applied.**

---

## AUC-direction summary (the load-bearing fact for replication)

| Dataset | AUC column | Range | Lower AUC = more sensitive? | Transform |
|---------|-----------|-------|-----------------------------|-----------|
| GDSC2   | `AUC`     | ~[0,1]| yes                         | none      |
| PRISM   | `auc`     | 0–4.9 | **yes**                     | none      |
| CTRP    | `auc`     | 0–1   | **yes**                     | none      |

All three are areas under a **viability** dose-response curve, so the sign of any
methylation↔AUC association is directly comparable across all three — a lead where
hypermethylation predicts *lower* GDSC AUC (more sensitive) should reproduce as
*lower* PRISM/CTRP AUC. No inversion needed for any source.
