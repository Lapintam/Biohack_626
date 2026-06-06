"""
Scan all drugs: how methylation-predictable is each drug's response (AUC),
cell-line-blind, vs the cancer-type-mean baseline?

For speed and rigor across ~286 drugs we use a FIXED, unsupervised gene panel:
the top-`N_GENES` most-variable genes across the whole methylation matrix,
selected once (label-blind, drug-blind -> no leakage into any drug's CV). Each
drug then gets a 5-fold cell-line-blind RidgeCV on that panel.

Outputs outputs/all_drugs_methylation_scan.csv ranked by Pearson r gain over the
cancer-type-mean baseline.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.model_selection import KFold
from sklearn.linear_model import RidgeCV

from strata.data import gdsc
from strata.config import OUTPUT_DIR

SEED = 1337
N_GENES = 3000
MIN_LINES = 40
FOLDS = 5
ALPHAS = np.logspace(-1, 4, 12)


def metrics(y_true, y_pred):
    rmse = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
    if len(y_true) <= 2 or np.ptp(y_pred) == 0:
        return rmse, float("nan")
    return rmse, float(pearsonr(y_true, y_pred)[0])


def main():
    meth = gdsc.load_gdsc_methylation()
    drug = gdsc.load_gdsc_drug_response()
    anno = gdsc.load_gdsc_annotations()

    # fixed, unsupervised gene panel (selected once across all lines)
    var = meth.var(axis=0).to_numpy()
    panel = np.argsort(var)[::-1][:N_GENES]
    Mpanel = meth.iloc[:, panel]
    panel_idx = {cid: i for i, cid in enumerate(Mpanel.index)}
    Mvals = Mpanel.to_numpy(dtype=np.float64)
    # z-score the panel once (global); Ridge is scale-tolerant and this is feature-only
    Mz = (Mvals - Mvals.mean(0)) / np.where(Mvals.std(0) == 0, 1.0, Mvals.std(0))

    cancer_map = anno["cancer_type"].fillna("Unknown")

    kf = KFold(n_splits=FOLDS, shuffle=True, random_state=SEED)
    rows = []
    drugs = sorted(drug["drug_name"].unique())
    for di, dn in enumerate(drugs):
        d = drug[drug["drug_name"] == dn]
        auc = d.groupby("COSMIC_ID")["auc"].mean()
        ids = [c for c in auc.index if c in panel_idx]
        if len(ids) < MIN_LINES:
            continue
        ridx = np.array([panel_idx[c] for c in ids])
        X = Mz[ridx]
        y = auc.loc[ids].to_numpy(dtype=np.float64)
        cancer = cancer_map.reindex(ids).fillna("Unknown").to_numpy()

        base_rmse, base_r, meth_rmse, meth_r = [], [], [], []
        for tr, te in kf.split(X):
            s = pd.Series(y[tr]).groupby(pd.Series(cancer[tr])).mean()
            pct = np.array([s.get(c, y[tr].mean()) for c in cancer[te]])
            br = metrics(y[te], pct); base_rmse.append(br[0]); base_r.append(br[1])
            m = RidgeCV(alphas=ALPHAS).fit(X[tr], y[tr])
            mr = metrics(y[te], m.predict(X[te])); meth_rmse.append(mr[0]); meth_r.append(mr[1])

        rows.append({
            "drug": dn, "n_lines": len(ids),
            "auc_std": float(np.std(y)),
            "base_rmse": float(np.mean(base_rmse)),
            "base_r": float(np.nanmean(base_r)),
            "meth_rmse": float(np.mean(meth_rmse)),
            "meth_r": float(np.nanmean(meth_r)),
        })
        if (di + 1) % 40 == 0:
            print(f"  ...{di+1}/{len(drugs)} drugs")

    res = pd.DataFrame(rows)
    res["r_gain"] = res["meth_r"] - res["base_r"]
    res["rmse_gain_pct"] = (res["base_rmse"] - res["meth_rmse"]) / res["base_rmse"] * 100
    res = res.sort_values("r_gain", ascending=False).reset_index(drop=True)

    out = os.path.join(OUTPUT_DIR, "all_drugs_methylation_scan.csv")
    res.to_csv(out, index=False)

    fmt = {"auc_std": "{:.3f}".format, "base_rmse": "{:.4f}".format,
           "base_r": "{:+.3f}".format, "meth_rmse": "{:.4f}".format,
           "meth_r": "{:+.3f}".format, "r_gain": "{:+.3f}".format,
           "rmse_gain_pct": "{:+.1f}".format}
    print(f"\nScanned {len(res)} drugs (>= {MIN_LINES} lines). Saved {out}\n")
    print("TOP 15 most methylation-predictable (by Pearson r gain over tissue baseline):")
    print(res.head(15).to_string(index=False, formatters=fmt))
    print(f"\nmean meth_r={res.meth_r.mean():.3f}  mean base_r={res.base_r.mean():.3f}  "
          f"drugs where methylation beats baseline: {(res.r_gain>0).sum()}/{len(res)}")
    # where is palbociclib?
    if (res.drug == "Palbociclib").any():
        rank = res.index[res.drug == "Palbociclib"][0]
        print(f"Palbociclib rank: {rank+1}/{len(res)}  "
              f"(meth_r={res.loc[rank,'meth_r']:+.3f}, base_r={res.loc[rank,'base_r']:+.3f})")


if __name__ == "__main__":
    main()
