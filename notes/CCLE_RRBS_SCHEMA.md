# CCLE RRBS promoter (TSS±1kb) methylation — schema & provenance

## Source

- **File:** `CCLE_RRBS_TSS1kb_20181022.txt.gz` (40,101,337 bytes; md5 `a228ad21f1016fb6508f6252d5187a2a`)
- **Release:** CCLE 2019 (Ghandi et al., *Nature* 2019, "Next-generation characterization of the
  Cancer Cell Line Encyclopedia"). RRBS promoter methylation in the **TSS ± 1 kb** window.
- **How obtained:** The DepMap portal proxies these as **time-limited GCS signed URLs**. The stable
  way to resolve a current link is the DepMap file manifest CSV:
  `https://depmap.org/portal/api/download/files` (a CSV of `release,release_date,filename,url,md5_hash`).
  Grep that for `CCLE_RRBS_TSS1kb_20181022.txt.gz`; the `url` column is a signed
  `https://storage.googleapis.com/depmap-external-downloads/ccle/ccle_2019/CCLE_RRBS_TSS1kb_20181022.txt.gz?GoogleAccessId=...&Expires=...&Signature=...`
  link. (The bare `Expires`-stamped URL stops working after expiry — re-fetch the manifest to refresh.)
  Stored locally at `data/ccle/` (gitignored).
- The legacy `data.broadinstitute.org/ccle/...` and the undated `...20180614...` names now 404 / are not
  the served version; the served promoter file is the `20181022` re-process. We use it.

## Raw shape

- **21,338 data rows** + 1 header row.
- **846 columns:** `locus_id`, `CpG_sites_hg19`, `avg_coverage`, then **843 cell-line columns**.

### Columns

| column            | meaning                                                                    |
|-------------------|----------------------------------------------------------------------------|
| `locus_id`        | `GENE_chr_start_end`, e.g. `MGMT_10_131264949_131265949` (TSS±1kb window). One `NA_NA` placeholder row exists and is dropped. |
| `CpG_sites_hg19`  | `;`-joined hg19 CpG coordinates covered in that window                      |
| `avg_coverage`    | mean RRBS read coverage over the window                                     |
| cell-line cols    | **CCLE_ID** = `NAME_TISSUE` (e.g. `DMS53_LUNG`, `PK59_PANCREAS`). No `ACH-` IDs in this file. |

### Values

- **Methylation ratio (beta) in [0, 1]** — coverage-weighted average of per-CpG methylation in the
  TSS±1kb window. `1.00000` and `0.00000` are common (fully methylated / unmethylated).
- **Missing = `NaN`** (literal token in the file; pandas reads as NaN). These are windows with no/low
  coverage in that line.

## Crosswalk to COSMIC_ID

- Map each cell-line column **CCLE_ID → COSMIC_ID** via `data/gdsc/model_list_20260420.csv`
  (`CCLE_ID` → `COSMIC_ID`).
- Of 843 cell-line columns: **835 match a CCLE_ID** in the model list; **614 map to a non-null,
  unique COSMIC_ID**. Unmapped columns are dropped (cap on the join n).

## Reduction to gene level

- `locus_id` parses as `^(gene)_([0-9XYMT]+)_(start)_(end)$`. The leading token is the **gene symbol**.
- **17,181 unique gene symbols.** When a gene has multiple TSS windows (alternative promoters / isoforms),
  we take the **mean** beta across its windows (NaN-skipping).
- All 5 demo genes present as loci: MGMT, CDKN2A, MTAP, INPP4B, SFRP1.

## Loader output (`strata.data.ccle_methylation.load_ccle_promoter_methylation`)

- **Drop-in for `gdsc.load_gdsc_methylation`:** rows = **COSMIC_ID (str)**, cols = **gene symbols (str)**,
  values = promoter methylation beta in [0, 1]; duplicate COSMIC_IDs deduped (keep first).
- Resulting shape: **614 lines × 17,181 genes**. Values remain `NaN` where the line had no coverage
  for any of a gene's promoter windows (documented; not imputed — unlike the GDSC gene-level matrix
  which is imputed). Downstream code should treat NaN as missing.

## Join coverage vs GDSC drug response

- COSMIC_IDs overlapping `gdsc.load_gdsc_drug_response()`: see `tests/test_ccle.py` / final report.
