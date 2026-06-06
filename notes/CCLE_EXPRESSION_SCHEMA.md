# CCLE RNA-seq expression — schema & provenance

## Source

- **File:** `CCLE_expression_full.csv` (646,185,698 bytes; md5 `f40d4734899efc99c0c40e77c8f79c63`)
- **Release:** DepMap Public **19Q4** (2019-12-19) / CCLE 2019 (Ghandi et al., *Nature* 2019,
  "Next-generation characterization of the Cancer Cell Line Encyclopedia"). RSEM RNA-seq.
- **How obtained:** Same DepMap manifest path as the RRBS file. The manifest CSV
  `https://depmap.org/portal/api/download/files` (`release,release_date,filename,url,md5_hash`)
  lists this file with a **figshare** `url`:
  `https://ndownloader.figshare.com/files/20234343`
  (figshare links are stable — unlike the GCS signed URLs used for RRBS, this one does not expire).
  Grep the manifest for `CCLE_expression_full.csv` + release `DepMap Public 19Q4`.
  Stored locally at `data/ccle/` (gitignored). md5 verified against the manifest.

## Raw shape

- **1,249 rows (cell lines) × 58,676 columns (genes).** Lines × genes orientation (rows = lines).
- **Row index:** `ACH-` DepMap IDs (`BROAD_ID`), e.g. `ACH-001097`. No COSMIC or CCLE_ID in the file.
- **Columns:** `SYMBOL (ENSG...)`, e.g. `TSPAN6 (ENSG00000000003)`. Includes all biotypes (~58k genes).

### Values

- **`log2(TPM + 1)`** — RSEM transcript-per-million, log-transformed (DepMap convention for the
  `*_expression*` files). Range observed ~0 to ~14. `0.0` = not expressed. **Native units kept;
  not rescaled.** (This is the "expression DOWN" axis for functional-silencing tests against beta.)

## Crosswalk to COSMIC_ID

- Map each line's **`ACH-` ID → COSMIC_ID** via `data/gdsc/model_list_20260420.csv` **`BROAD_ID`** column
  (`BROAD_ID` → `COSMIC_ID`).
- Of 1,249 `ACH-` lines: **733 map to a non-null, unique COSMIC_ID.** Unmapped lines are dropped.

## Reduction to gene level

- Strip the `(ENSG...)` suffix from each column to a **bare gene symbol** (regex
  `\s*\(ENSG\d+(?:\.\d+)?\)\s*$`).
- 58,676 columns → **57,138 unique symbols**; **172 symbols map to >1 Ensembl ID** and are collapsed by
  the **mean** across those columns.
- All 5 demo genes present: MGMT, CDKN2A, MTAP, INPP4B, SFRP1.

## Loader output (`strata.data.ccle_expression.load_ccle_expression`)

- **Drop-in companion to `ccle_methylation.load_ccle_promoter_methylation` / `gdsc.load_gdsc_methylation`:**
  rows = **COSMIC_ID (str)**, cols = **gene symbols (str)**; same orientation and key contract.
  Values = `log2(TPM + 1)` (expression, not beta). Duplicate COSMIC_IDs deduped (keep first).
- Resulting shape: **733 lines × 57,138 genes.**

## Join coverage (the n for silencing validation)

- **Expression ∩ CCLE promoter methylation** (`ccle_methylation.load_ccle_promoter_methylation`, by
  COSMIC_ID): **613 lines.** This is the n for the methylation↔expression functional-silencing test.
- **Expression ∩ GDSC drug response** (`gdsc.load_gdsc_drug_response`): **663 lines.**
- (Methylation has 614 lines total; 613 of them also have expression — near-complete overlap.)
