"""Loaders for TCGA-LAML (Acute Myeloid Leukemia) DNA methylation + survival.

Data provenance and schema are documented in ``notes/TCGA_SCHEMA.md``.

Source: UCSC Xena, GDC TCGA Acute Myeloid Leukemia (LAML) hub
(https://xenabrowser.net / https://xena.ucsc.edu). Two flat matrices:

  * ``TCGA-LAML.methylation450.tsv.gz`` -- Illumina HumanMethylation450 beta
    matrix, probes (rows) x samples (columns). Transposed here to samples x probes
    to match the engine convention (samples as ROWS).
  * ``TCGA-LAML.survival.tsv.gz`` -- curated survival table with columns
    ``sample``, ``OS.time`` (days), ``OS`` (1=dead, 0=censored), ``_PATIENT``.

Sample-id convention
--------------------
Both matrices key on the TCGA aliquot barcode at the sample level, e.g.
``TCGA-AB-2820-03A``. The methylation column headers and the survival ``sample``
column use the SAME barcodes, so they join directly with no trimming. The survival
table contains multiple aliquots per patient (``-03A``, ``-03B``, ``-03D`` with
identical OS); we keep the first occurrence per barcode and the intersection of
barcodes present in BOTH matrices drives the join (~190+ patients).

Missing-value handling (methylation)
------------------------------------
The 450K matrix has many NaNs (probes not measured in a given sample). In
``load_tcga_methylation`` we drop probes missing in >20% of samples, then fill any
remaining NaN with that probe's mean across samples, so downstream clustering /
variable-feature selection sees a dense matrix.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd

from strata.config import DATA_DIR

TCGA_DIR = DATA_DIR / "tcga_laml"

METHYLATION_FILE = TCGA_DIR / "TCGA-LAML.methylation450.tsv.gz"
SURVIVAL_FILE = TCGA_DIR / "TCGA-LAML.survival.tsv.gz"

# Drop a probe if more than this fraction of samples are missing it.
MAX_PROBE_MISSING = 0.20


@lru_cache(maxsize=1)
def load_tcga_methylation() -> pd.DataFrame:
    """Illumina 450K beta matrix, rows = sample barcode (str), cols = probe id (str).

    The source file is probes x samples; it is transposed to samples x probes.
    Probes missing in >20% of samples are dropped, then remaining NaNs are filled
    with the per-probe mean. Values are beta in [0, 1].
    """
    # Source: rows = probe (Composite Element REF), columns = sample barcodes.
    raw = pd.read_csv(METHYLATION_FILE, sep="\t", index_col=0)
    raw.index = raw.index.astype(str)
    raw.columns = raw.columns.astype(str)

    # Do all NaN handling in the native probes-x-samples orientation (probes are
    # ROWS here). pandas column-wise ops over hundreds of thousands of probes are
    # pathologically slow, so we operate row-wise / with numpy, then transpose
    # last. The result is identical to filtering after the transpose.

    # Drop probes (rows) missing in > MAX_PROBE_MISSING of samples (columns).
    keep = raw.isna().mean(axis=1) <= MAX_PROBE_MISSING
    raw = raw.loc[keep]

    # Fill remaining NaNs with the per-probe (row) mean.
    row_mean = raw.mean(axis=1).values
    vals = raw.values
    nan_idx = np.where(np.isnan(vals))
    vals[nan_idx] = np.take(row_mean, nan_idx[0])

    # Transpose -> samples (rows) x probes (cols), matching the engine convention.
    meth = pd.DataFrame(vals.T, index=raw.columns, columns=raw.index)
    meth.index.name = "sample"

    return meth


@lru_cache(maxsize=1)
def load_tcga_clinical() -> pd.DataFrame:
    """Overall-survival table, index = sample barcode (str), cols os_time, os_event.

    ``os_time`` = OS.time in days (float). ``os_event`` = OS (1=dead, 0=censored,
    int). Indexed by the same sample barcode used by the methylation matrix; the
    first aliquot per barcode is kept.
    """
    surv = pd.read_csv(SURVIVAL_FILE, sep="\t")
    surv["sample"] = surv["sample"].astype(str)

    # Use .values so the new sample-string index is not aligned against the
    # source RangeIndex (which would silently produce all-NaN columns).
    out = pd.DataFrame(
        {
            "os_time": pd.to_numeric(surv["OS.time"], errors="coerce").values,
            "os_event": pd.to_numeric(surv["OS"], errors="coerce").values,
        },
        index=pd.Index(surv["sample"].values, name="sample"),
    )

    # Need both fields; OS must be a clean 0/1 event flag and time non-negative.
    out = out.dropna(subset=["os_time", "os_event"])
    out = out[out["os_event"].isin([0, 1])]
    out = out[out["os_time"] >= 0]
    out["os_event"] = out["os_event"].astype(int)

    # Keep the first aliquot per barcode (barcodes are already sample-level).
    out = out[~out.index.duplicated(keep="first")]
    return out
