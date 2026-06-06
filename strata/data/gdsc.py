"""Loaders for GDSC / Cell Model Passports data.

Data provenance and schema are documented in ``notes/GDSC_SCHEMA.md``.

Key fact: the native join key across the GDSC2 drug file, the methylation matrix,
and the model list is the **Sanger model id** (``SANGER_MODEL_ID`` / ``model_id``,
e.g. ``SIDM00001``) -- the drug file has no COSMIC_ID column at all. These loaders
honor the COSMIC_ID-keyed contract by mapping SIDM -> COSMIC_ID via the model list
and dropping rows whose SIDM lacks a COSMIC_ID.
"""

from __future__ import annotations

from functools import lru_cache

import pandas as pd

from strata.config import DATA_DIR

GDSC_DIR = DATA_DIR / "gdsc"

METHYLATION_FILE = GDSC_DIR / "methylation_imputed.csv.gz"
METHYLATION_GENELEVEL_FILE = GDSC_DIR / "methylation_genelevel_bycosmic.csv.gz"
METHYLATION_PROMOTER_FILE = GDSC_DIR / "methylation_promoter_bycosmic.csv.gz"
DRUG_RESPONSE_FILE = GDSC_DIR / "GDSC2_fitted_dose_response_27Oct23.xlsx"
MODEL_LIST_FILE = GDSC_DIR / "model_list_20260420.csv"


@lru_cache(maxsize=1)
def _sidm_to_cosmic() -> "pd.Series":
    """Series mapping Sanger model id (SIDM) -> COSMIC_ID (str), non-null only."""
    ml = pd.read_csv(MODEL_LIST_FILE, dtype=str, low_memory=False)
    m = ml.set_index("model_id")["COSMIC_ID"]
    return m.dropna()


def load_gdsc_methylation() -> pd.DataFrame:
    """Gene-level methylation beta matrix keyed by COSMIC_ID.

    Returns a DataFrame whose rows are COSMIC_ID (str) and whose columns are
    gene symbols (str); values are beta in [0, 1].

    The source matrix is indexed by Sanger model id (SIDM); it is reindexed to
    COSMIC_ID via the model list. SIDMs without a COSMIC_ID are dropped.
    """
    df = pd.read_csv(METHYLATION_FILE, index_col=0)
    df.index = df.index.astype(str)

    mapping = _sidm_to_cosmic()
    keep = df.index.intersection(mapping.index)
    df = df.loc[keep]
    df.index = mapping.loc[keep].astype(str).values
    df.index.name = "COSMIC_ID"

    # Guard against accidental duplicate COSMIC_IDs (keep first).
    df = df[~df.index.duplicated(keep="first")]
    return df


def _load_cosmic_keyed_betas(path, dropna_cols: bool = False) -> pd.DataFrame:
    """Load a COSMIC_ID-keyed beta matrix (gene-symbol columns) from `path`.

    These matrices are already indexed by COSMIC_ID, so no SIDM crosswalk is
    needed. Rows are COSMIC_ID (str), columns are gene symbols. The promoter
    matrix contains NaNs (un-imputed); set ``dropna_cols`` to drop all-NaN
    columns. Values are beta (~[0, 1]).
    """
    df = pd.read_csv(path, index_col=0, low_memory=False)
    df.index = df.index.astype(str)
    df.index.name = "COSMIC_ID"
    df = df[~df.index.duplicated(keep="first")]
    if dropna_cols:
        df = df.dropna(axis=1, how="all")
    return df


def load_gdsc_methylation_genelevel() -> pd.DataFrame:
    """Gene-level methylation beta matrix, natively keyed by COSMIC_ID.

    Rows = COSMIC_ID (str), columns = gene symbols (~14.6k), values = beta.
    Equivalent resolution to :func:`load_gdsc_methylation` but distributed
    already COSMIC-keyed (no SIDM->COSMIC reindex).
    """
    return _load_cosmic_keyed_betas(METHYLATION_GENELEVEL_FILE)


def load_gdsc_methylation_promoter() -> pd.DataFrame:
    """Promoter-region methylation beta matrix, keyed by COSMIC_ID.

    Rows = COSMIC_ID (str), columns = gene symbols (~17.2k promoters), values =
    beta in [0, 1] with NaNs where a promoter was not measured (NOT imputed).

    Promoter-CpG resolution preserves silencing signals that gene-level
    averaging washes out (e.g. MGMT promoter methylation). Callers should
    handle NaNs (e.g. impute per-column or restrict to covered genes).
    """
    return _load_cosmic_keyed_betas(METHYLATION_PROMOTER_FILE)


def load_gdsc_drug_response() -> pd.DataFrame:
    """Long-format GDSC2 fitted dose-response.

    Columns: ``["COSMIC_ID", "drug_name", "ln_ic50", "auc"]``. COSMIC_ID is a str,
    mapped from the file's ``SANGER_MODEL_ID``; rows with no COSMIC_ID are dropped.
    """
    raw = pd.read_excel(
        DRUG_RESPONSE_FILE,
        usecols=["SANGER_MODEL_ID", "DRUG_NAME", "LN_IC50", "AUC"],
    )

    mapping = _sidm_to_cosmic()
    sidm = raw["SANGER_MODEL_ID"].astype(str)
    cosmic = sidm.map(mapping)

    out = pd.DataFrame(
        {
            "COSMIC_ID": cosmic,
            "drug_name": raw["DRUG_NAME"].astype(str),
            "ln_ic50": pd.to_numeric(raw["LN_IC50"], errors="coerce"),
            "auc": pd.to_numeric(raw["AUC"], errors="coerce"),
        }
    )
    out = out.dropna(subset=["COSMIC_ID"]).reset_index(drop=True)
    out["COSMIC_ID"] = out["COSMIC_ID"].astype(str)
    return out


def load_gdsc_annotations() -> pd.DataFrame:
    """Cell-line annotations indexed by COSMIC_ID.

    Index = COSMIC_ID (str). Columns include ``cell_line_name`` and ``tissue``
    (plus ``sanger_model_id`` and ``cancer_type``).
    """
    ml = pd.read_csv(MODEL_LIST_FILE, dtype=str, low_memory=False)
    ml = ml[ml["COSMIC_ID"].notna()].copy()

    out = pd.DataFrame(
        {
            "COSMIC_ID": ml["COSMIC_ID"].astype(str),
            "sanger_model_id": ml["model_id"].astype(str),
            "cell_line_name": ml["model_name"],
            "tissue": ml["tissue"],
            "cancer_type": ml["cancer_type"],
        }
    )
    out = out[~out["COSMIC_ID"].duplicated(keep="first")]
    out = out.set_index("COSMIC_ID")
    return out
