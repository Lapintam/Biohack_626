"""Loaders for orthogonal (non-GDSC) cell-line drug-response screens.

These exist to **replicate** methylation->drug-response leads first found in GDSC2
against independent screens. Both loaders return the same LONG schema as
``strata.data.gdsc.load_gdsc_drug_response``:

    columns = ["COSMIC_ID", "drug_name", "auc"]

with ``COSMIC_ID`` as str and ``drug_name`` lowercased/stripped.

Provenance, file URLs, shapes, and AUC semantics are documented in
``notes/EXTERNAL_DRUGRESP_SCHEMA.md``.

AUC direction (critical for replication): in **every** source here -- GDSC2,
PRISM, and CTRP -- ``auc`` is the area under a *viability* dose-response curve, so
**lower auc = more sensitive**. No source is inverted, so no transform is applied;
the sign of a methylation<->auc association is directly comparable to GDSC.

Join strategy: external IDs are DepMap/Broad ``ACH-...`` ids, mapped to
``COSMIC_ID`` via the ``BROAD_ID`` column of the Cell Model Passports model list
(``data/gdsc/model_list_20260420.csv``). Lines that do not map are dropped.
"""

from __future__ import annotations

from functools import lru_cache

import pandas as pd

from strata.config import DATA_DIR

PRISM_DIR = DATA_DIR / "prism"
CTRP_DIR = DATA_DIR / "ctrp"
MODEL_LIST_FILE = DATA_DIR / "gdsc" / "model_list_20260420.csv"

PRISM_DOSE_RESPONSE_FILE = (
    PRISM_DIR / "secondary-screen-dose-response-curve-parameters.csv"
)
CTRP_AUC_MATRIX_FILE = CTRP_DIR / "CTRPAUCMatrix.csv"
CTRP_RESPONSE_CURVES_FILE = CTRP_DIR / "CTRPResponseCurves.csv"


@lru_cache(maxsize=1)
def _broad_to_cosmic() -> "pd.Series":
    """Series mapping DepMap/Broad id (``ACH-...``) -> COSMIC_ID (str).

    Non-null and de-duplicated on ``BROAD_ID`` (keep first).
    """
    ml = pd.read_csv(MODEL_LIST_FILE, dtype=str, low_memory=False)
    ml = ml.dropna(subset=["BROAD_ID", "COSMIC_ID"])
    m = ml.set_index("BROAD_ID")["COSMIC_ID"]
    m = m[~m.index.duplicated(keep="first")]
    return m


def _to_long_cosmic(
    df: pd.DataFrame, *, broad_col: str, drug_col: str, auc_col: str
) -> pd.DataFrame:
    """Map a (broad_id, drug, auc) frame to the canonical LONG COSMIC schema.

    - maps ``broad_col`` (``ACH-...``) -> COSMIC_ID, dropping unmapped rows;
    - lowercases/strips ``drug_col``;
    - coerces ``auc`` to numeric and drops NaN;
    - collapses duplicate (COSMIC_ID, drug_name) pairs by mean AUC.
    """
    mapping = _broad_to_cosmic()
    out = pd.DataFrame(
        {
            "COSMIC_ID": df[broad_col].astype(str).map(mapping),
            "drug_name": df[drug_col].astype(str).str.strip().str.lower(),
            "auc": pd.to_numeric(df[auc_col], errors="coerce"),
        }
    )
    out = out.dropna(subset=["COSMIC_ID", "auc"])
    out = out[out["drug_name"].str.len() > 0]
    out["COSMIC_ID"] = out["COSMIC_ID"].astype(str)
    # One value per (line, drug): mean across screens/doses/replicates.
    out = (
        out.groupby(["COSMIC_ID", "drug_name"], as_index=False)["auc"]
        .mean()
        .reset_index(drop=True)
    )
    return out


def load_prism_drug_response() -> pd.DataFrame:
    """PRISM Repurposing secondary-screen dose-response, LONG format.

    Columns: ``["COSMIC_ID", "drug_name", "auc"]`` (COSMIC_ID str).

    Source: ``secondary-screen-dose-response-curve-parameters.csv`` (DepMap
    "PRISM Secondary Repurposing 20Q2"). The native key is ``depmap_id``
    (``ACH-...``), mapped to COSMIC_ID via ``BROAD_ID``; unmapped lines and
    control-barcode rows (null ``depmap_id``) are dropped.

    ``auc`` is the area under the fitted viability dose-response curve;
    **lower auc = more sensitive** (same direction as GDSC2). No transform.
    Multiple screens (HTS002/MTS010/MTS006/MTS005) are collapsed to one mean AUC
    per (line, drug).
    """
    raw = pd.read_csv(
        PRISM_DOSE_RESPONSE_FILE,
        usecols=["depmap_id", "name", "auc"],
        low_memory=False,
    )
    return _to_long_cosmic(raw, broad_col="depmap_id", drug_col="name", auc_col="auc")


def load_ctrp_drug_response() -> pd.DataFrame:
    """CTRP (CTD^2) v2 AUC, LONG format.

    Columns: ``["COSMIC_ID", "drug_name", "auc"]`` (COSMIC_ID str).

    Source: DepMap "Harmonized CTD^2 25Q2" -- ``CTRPAUCMatrix.csv`` (wide:
    rows = ``ModelID`` ``ACH-...``, cols = ``CompoundID`` ``DPC-...``) melted to
    long; ``CompoundID`` -> ``CompoundName`` via ``CTRPResponseCurves.csv``.
    ModelID (``ACH-...``) is mapped to COSMIC_ID via ``BROAD_ID``; unmapped lines
    are dropped.

    ``auc`` is the harmonized normalized area under the viability curve, in [0, 1];
    **lower auc = more sensitive** (same direction as GDSC2/PRISM). No transform.

    Raises ``NotImplementedError`` if the CTRP files are not present locally.
    """
    if not (CTRP_AUC_MATRIX_FILE.exists() and CTRP_RESPONSE_CURVES_FILE.exists()):
        raise NotImplementedError(
            "CTRP files not found. Expected "
            f"{CTRP_AUC_MATRIX_FILE} and {CTRP_RESPONSE_CURVES_FILE}. "
            "Download the DepMap 'Harmonized CTD^2 25Q2' CTRPAUCMatrix.csv and "
            "CTRPResponseCurves.csv (see notes/EXTERNAL_DRUGRESP_SCHEMA.md)."
        )

    auc_wide = pd.read_csv(CTRP_AUC_MATRIX_FILE, index_col=0)
    auc_wide.index.name = "ModelID"

    # CompoundID -> CompoundName.
    curves = pd.read_csv(
        CTRP_RESPONSE_CURVES_FILE, usecols=["CompoundName", "CompoundID"]
    )
    cid_to_name = (
        curves.dropna(subset=["CompoundID"])
        .drop_duplicates("CompoundID")
        .set_index("CompoundID")["CompoundName"]
    )

    long = (
        auc_wide.reset_index()
        .melt(id_vars="ModelID", var_name="CompoundID", value_name="auc")
        .dropna(subset=["auc"])
    )
    long["name"] = long["CompoundID"].map(cid_to_name)
    long = long.dropna(subset=["name"])

    return _to_long_cosmic(long, broad_col="ModelID", drug_col="name", auc_col="auc")
