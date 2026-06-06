"""Loader for CCLE RRBS promoter (TSS±1kb) DNA-methylation.

Data provenance and schema are documented in ``notes/CCLE_RRBS_SCHEMA.md``.

This is a **drop-in replacement** for ``strata.data.gdsc.load_gdsc_methylation``:
both return a beta matrix whose rows are ``COSMIC_ID`` (str) and whose columns are
gene symbols (str), values in ``[0, 1]``. The difference is resolution -- this one
is promoter/TSS-region methylation from CCLE RRBS rather than gene-level imputed
values, so promoter-specific signal (e.g. MGMT) is preserved.

Source: ``CCLE_RRBS_TSS1kb_20181022.txt.gz`` (CCLE 2019 / Ghandi et al. 2019),
methylation in the TSS +/- 1 kb window. Rows are ``GENE_chr_start_end`` loci, columns
are CCLE_IDs (``NAME_TISSUE``). We map CCLE_ID -> COSMIC_ID via the Cell Model
Passports model list, reduce multiple promoter windows per gene by the mean, and
transpose to lines x genes.
"""

from __future__ import annotations

import re
from functools import lru_cache

import pandas as pd

from strata.config import DATA_DIR

CCLE_DIR = DATA_DIR / "ccle"
RRBS_TSS1KB_FILE = CCLE_DIR / "CCLE_RRBS_TSS1kb_20181022.txt.gz"
MODEL_LIST_FILE = DATA_DIR / "gdsc" / "model_list_20260420.csv"

# locus_id looks like "MGMT_10_131264949_131265949"; the leading token is the gene
# symbol, followed by chrom (1-22/X/Y/MT) and start/end. Gene symbols may themselves
# contain underscores, so anchor on the trailing _chrom_start_end.
_LOCUS_RE = re.compile(r"^(?P<gene>.+)_(?P<chrom>[0-9XYMT]+)_(?P<start>\d+)_(?P<end>\d+)$")


@lru_cache(maxsize=1)
def _ccle_to_cosmic() -> "pd.Series":
    """Series mapping CCLE_ID (NAME_TISSUE) -> COSMIC_ID (str), non-null & unique."""
    ml = pd.read_csv(MODEL_LIST_FILE, dtype=str, low_memory=False)
    ml = ml.dropna(subset=["CCLE_ID", "COSMIC_ID"])
    m = ml.set_index("CCLE_ID")["COSMIC_ID"]
    m = m[~m.index.duplicated(keep="first")]
    return m


def _gene_of(locus_id: str) -> "str | None":
    m = _LOCUS_RE.match(locus_id)
    return m.group("gene") if m else None


def load_ccle_promoter_methylation() -> pd.DataFrame:
    """Promoter (TSS+/-1kb) methylation beta matrix keyed by COSMIC_ID.

    Returns a DataFrame whose rows are ``COSMIC_ID`` (str) and whose columns are
    gene symbols (str); values are promoter methylation beta in ``[0, 1]``.

    Missing windows are kept as ``NaN`` (the RRBS file leaves uncovered windows
    empty; values are *not* imputed). Genes with multiple promoter windows are
    collapsed by the NaN-skipping mean across windows.
    """
    raw = pd.read_csv(RRBS_TSS1KB_FILE, sep="\t", low_memory=False)

    # Drop the annotation columns; keep locus_id as the row index.
    raw = raw.set_index("locus_id").drop(columns=["CpG_sites_hg19", "avg_coverage"])

    # Coerce all cell-line columns to numeric (NaN tokens -> NaN).
    raw = raw.apply(pd.to_numeric, errors="coerce")

    # Map locus -> gene symbol; drop loci that don't parse (e.g. the "NA_NA" row).
    genes = pd.Index(raw.index).map(_gene_of)
    raw = raw[pd.notna(genes)]
    raw.index = pd.Index([g for g in genes if g is not None], name="gene")

    # Collapse multiple promoter windows per gene -> mean beta (skip NaN).
    by_gene = raw.groupby(level=0).mean()  # genes x cell-lines

    # Transpose to lines x genes.
    mat = by_gene.T  # rows = CCLE_ID, cols = gene

    # Map CCLE_ID -> COSMIC_ID; drop unmapped lines.
    mapping = _ccle_to_cosmic()
    keep = mat.index.intersection(mapping.index)
    mat = mat.loc[keep]
    mat.index = mapping.loc[keep].astype(str).values
    mat.index.name = "COSMIC_ID"

    # Dedupe any accidental duplicate COSMIC_IDs (keep first).
    mat = mat[~mat.index.duplicated(keep="first")]

    # Stable column order; ensure column labels are str.
    mat.columns = mat.columns.astype(str)
    mat = mat.sort_index(axis=1)
    return mat
