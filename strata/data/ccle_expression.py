"""Loader for CCLE RNA-seq gene expression (log2(TPM+1)).

Data provenance and schema are documented in ``notes/CCLE_EXPRESSION_SCHEMA.md``.

This is a **drop-in companion** to ``strata.data.ccle_methylation.load_ccle_promoter_methylation``
and ``strata.data.gdsc.load_gdsc_methylation``: it returns a matrix whose rows are
``COSMIC_ID`` (str) and whose columns are gene symbols (str), same orientation/keys.
Where those return methylation beta, this returns RNA-seq expression so methylation can
be tested against expression (functional silencing: promoter methylation up -> expression down).

Source: ``CCLE_expression_full.csv`` (DepMap Public 19Q4 / CCLE 2019, Ghandi et al. 2019),
RSEM ``log2(TPM + 1)`` over all genes. Rows are ``ACH-`` (DepMap) cell-line IDs, columns are
``SYMBOL (ENSG...)``. We strip the Ensembl suffix to bare gene symbols, collapse genes that map
to multiple Ensembl IDs by the mean, and map ``ACH-`` -> ``COSMIC_ID`` via the Cell Model
Passports model list (``BROAD_ID`` column).
"""

from __future__ import annotations

import re
from functools import lru_cache

import pandas as pd

from strata.config import DATA_DIR

CCLE_DIR = DATA_DIR / "ccle"
EXPRESSION_FILE = CCLE_DIR / "CCLE_expression_full.csv"
MODEL_LIST_FILE = DATA_DIR / "gdsc" / "model_list_20260420.csv"

# Columns look like "TSPAN6 (ENSG00000000003)"; strip the trailing "(ENSG...)" to a bare symbol.
_ENSG_SUFFIX_RE = re.compile(r"\s*\(ENSG\d+(?:\.\d+)?\)\s*$")


@lru_cache(maxsize=1)
def _broad_to_cosmic() -> "pd.Series":
    """Series mapping DepMap BROAD_ID (ACH-...) -> COSMIC_ID (str), non-null & unique."""
    ml = pd.read_csv(MODEL_LIST_FILE, dtype=str, low_memory=False)
    ml = ml.dropna(subset=["BROAD_ID", "COSMIC_ID"])
    m = ml.set_index("BROAD_ID")["COSMIC_ID"]
    m = m[~m.index.duplicated(keep="first")]
    return m


def _strip_ensg(col: str) -> str:
    return _ENSG_SUFFIX_RE.sub("", str(col)).strip()


def load_ccle_expression() -> pd.DataFrame:
    """CCLE RNA-seq expression matrix keyed by COSMIC_ID.

    Returns a DataFrame whose rows are ``COSMIC_ID`` (str) and whose columns are gene
    symbols (str); values are ``log2(TPM + 1)`` (native CCLE units; *not* rescaled).

    Cell lines are keyed by ``ACH-`` (DepMap) IDs in the source and mapped to ``COSMIC_ID``
    via the model list ``BROAD_ID`` column; unmapped lines are dropped. Genes that map to
    multiple Ensembl IDs (same bare symbol) are collapsed by the mean across those columns.
    Duplicate COSMIC_IDs are deduped (keep first).
    """
    df = pd.read_csv(EXPRESSION_FILE, index_col=0)
    df.index = df.index.astype(str)

    # Strip "(ENSG...)" suffixes to bare gene symbols.
    df.columns = [_strip_ensg(c) for c in df.columns]

    # Collapse genes with multiple Ensembl IDs (duplicate bare symbols) -> mean expression.
    if df.columns.duplicated().any():
        df = df.T.groupby(level=0).mean().T

    # Map ACH- (BROAD_ID) -> COSMIC_ID; drop unmapped lines.
    mapping = _broad_to_cosmic()
    keep = df.index.intersection(mapping.index)
    df = df.loc[keep]
    df.index = mapping.loc[keep].astype(str).values
    df.index.name = "COSMIC_ID"

    # Dedupe any accidental duplicate COSMIC_IDs (keep first).
    df = df[~df.index.duplicated(keep="first")]

    # Stable column order; ensure column labels are str.
    df.columns = df.columns.astype(str)
    df = df.sort_index(axis=1)
    return df
