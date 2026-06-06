"""
Link GDSC2 drug-response labels to cell-line metadata and emit a clean,
NN-ready dataset.

The two files share a cell-line identifier:
    GDSC2_fitted_dose_response_27Oct23.xlsx :: SANGER_MODEL_ID
    model_list_20260420.csv                 :: model_id

The actual row-to-row LINK is done in numpy with `np.searchsorted` (a vectorized
key join), not pandas.merge -- pandas is used only to read the spreadsheets.
Each response row (one cell line x one drug) is matched to exactly one metadata
row, so the join is one-to-many and lossless on the response side.

Outputs (written next to this script):
    gdsc2_merged.csv   -- merged + cleaned long table (labels + curated metadata)
    gdsc2_merged.npz   -- numpy arrays ready for a neural net:
                            X_num   (float32)  standardized numeric features
                            X_cat   (int64)    integer-encoded categorical features
                            y       (float32)  LN_IC50 label
                            plus the encoders/column names for reconstruction
"""

import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RESPONSE_XLSX = os.path.join(HERE, "GDSC2_fitted_dose_response_27Oct23.xlsx")
MODEL_CSV = os.path.join(HERE, "model_list_20260420.csv")
OUT_CSV = os.path.join(HERE, "gdsc2_merged.csv")
OUT_NPZ = os.path.join(HERE, "gdsc2_merged.npz")

# Shared key.
RESP_KEY = "SANGER_MODEL_ID"
MODEL_KEY = "model_id"

# Response columns we carry through (labels + the join key + drug identity).
RESP_COLS = ["SANGER_MODEL_ID", "CELL_LINE_NAME", "DRUG_ID", "DRUG_NAME",
             "PUTATIVE_TARGET", "PATHWAY_NAME", "LN_IC50", "AUC", "Z_SCORE"]

# Curated, well-populated metadata features (see fill-rate analysis). Sparse
# clinical fields are dropped so every model sees the same columns.
META_CATEGORICAL = ["tissue", "cancer_type", "tissue_status", "growth_properties",
                    "gender", "ethnicity", "msi_status", "smoking_status"]
META_NUMERIC = ["mutational_burden", "ploidy_wes", "age_at_sampling"]
# Identifiers kept for reference only (not used as model features).
META_IDS = ["COSMIC_ID", "model_name"]


def numpy_link(resp_key_values, model_key_values):
    """Map each response key to the row index of its metadata record, in numpy.

    Returns an int array `idx` such that metadata_row[idx[i]] corresponds to
    response_row[i]. Raises if any response key is missing from the metadata.
    """
    resp_key_values = np.asarray(resp_key_values)
    model_key_values = np.asarray(model_key_values)

    # Sort the metadata keys so we can binary-search into them.
    order = np.argsort(model_key_values, kind="stable")
    sorted_keys = model_key_values[order]

    # Vectorized lookup: position of each response key within the sorted metadata.
    pos = np.searchsorted(sorted_keys, resp_key_values)
    pos = np.clip(pos, 0, len(sorted_keys) - 1)

    matched = sorted_keys[pos] == resp_key_values
    if not matched.all():
        missing = np.unique(resp_key_values[~matched])
        raise KeyError(f"{len(missing)} response key(s) absent from metadata, "
                       f"e.g. {missing[:5].tolist()}")

    # Translate positions in the sorted view back to original metadata rows.
    return order[pos]


def main():
    print("Reading files ...")
    resp = pd.read_excel(RESPONSE_XLSX)
    models = pd.read_csv(MODEL_CSV, low_memory=False)

    resp = resp[[c for c in RESP_COLS if c in resp.columns]].copy()
    meta_cols = [MODEL_KEY] + META_CATEGORICAL + META_NUMERIC + META_IDS
    models = models[[c for c in meta_cols if c in models.columns]].copy()

    # --- clean each side before linking -------------------------------------
    n0 = len(resp)
    resp = resp.dropna(subset=["LN_IC50", RESP_KEY, "DRUG_ID"])
    resp = resp.drop_duplicates(subset=[RESP_KEY, "DRUG_ID"], keep="first")
    print(f"Response rows: {n0:,} -> {len(resp):,} after dropping null labels / dup pairs")

    # metadata key must be unique for a clean one-to-many link
    models = models.drop_duplicates(subset=[MODEL_KEY], keep="first")

    # --- THE LINK (numpy) ---------------------------------------------------
    idx = numpy_link(resp[RESP_KEY].to_numpy(), models[MODEL_KEY].to_numpy())
    print(f"Linked {len(idx):,} response rows to metadata via np.searchsorted")

    # Gather metadata aligned to every response row.
    meta_aligned = models.iloc[idx].reset_index(drop=True)
    resp = resp.reset_index(drop=True)

    merged = pd.concat(
        [resp, meta_aligned.drop(columns=[MODEL_KEY])], axis=1
    )

    # --- consistent missing-value handling (apples-to-apples) ---------------
    for c in META_CATEGORICAL:
        if c in merged:
            merged[c] = merged[c].fillna("Unknown").astype(str).str.strip()
    for c in META_NUMERIC:
        if c in merged:
            merged[c] = pd.to_numeric(merged[c], errors="coerce")
            merged[c] = merged[c].fillna(merged[c].median())

    merged.to_csv(OUT_CSV, index=False)
    print(f"Wrote {OUT_CSV}  ({merged.shape[0]:,} rows x {merged.shape[1]} cols)")

    # --- build numpy tensors for the neural net -----------------------------
    cat_present = [c for c in META_CATEGORICAL if c in merged]
    num_present = [c for c in META_NUMERIC if c in merged]

    # integer-encode categoricals (store the category lists for reconstruction)
    cat_codes = np.zeros((len(merged), len(cat_present)), dtype=np.int64)
    cat_vocab = {}
    for j, c in enumerate(cat_present):
        cats, codes = np.unique(merged[c].to_numpy(), return_inverse=True)
        cat_codes[:, j] = codes
        cat_vocab[c] = cats

    # standardize numerics (z-score over the whole table; re-fit on train in NN)
    X_num = merged[num_present].to_numpy(dtype=np.float32)
    mu = X_num.mean(axis=0)
    sd = X_num.std(axis=0)
    sd[sd == 0] = 1.0
    X_num = (X_num - mu) / sd

    y = merged["LN_IC50"].to_numpy(dtype=np.float32)

    np.savez_compressed(
        OUT_NPZ,
        X_num=X_num,
        X_cat=cat_codes,
        y=y,
        num_cols=np.array(num_present),
        cat_cols=np.array(cat_present),
        sanger_model_id=merged[RESP_KEY].to_numpy().astype(str),
        drug_id=merged["DRUG_ID"].to_numpy(),
        **{f"vocab__{c}": cat_vocab[c].astype(str) for c in cat_present},
    )
    print(f"Wrote {OUT_NPZ}")
    print(f"  X_num {X_num.shape} (cols: {num_present})")
    print(f"  X_cat {cat_codes.shape} (cols: {cat_present})")
    print(f"  y     {y.shape}  mean={y.mean():.3f} std={y.std():.3f}")
    print(f"  cell lines: {merged[RESP_KEY].nunique()} | drugs: {merged['DRUG_ID'].nunique()}")


if __name__ == "__main__":
    main()
