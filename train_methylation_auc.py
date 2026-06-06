"""
Predict drug response (AUC) from gene methylation (+ cancer type) for one drug,
evaluated cell-line-blind.

Setup
-----
One row per cell line (Beta over ~14.6k genes) -> AUC for a chosen drug. Because
each cell line appears once, a plain K-fold over rows IS a cell-line-blind split:
every test line is unseen in training.

p >> n (≈14.6k genes, ≈1k lines), so we use regularized linear models with
train-only feature selection + standardization (no leakage). The honest bar is
the cancer-type-mean baseline, since methylation strongly tracks tissue. We
report four nested models so each added ingredient is auditable:

    1. global-mean AUC          (trivial floor)
    2. cancer-type-mean AUC     (the real bar -- pure tissue effect)
    3. methylation only         (Ridge on top-variance genes)
    4. methylation + cancer_type (genes + tissue one-hot)

A model only "counts" if it beats the cancer-type-mean baseline.

Usage:
    python train_methylation_auc.py --drug Palbociclib
    python train_methylation_auc.py --drug Temozolomide --n-genes 3000
"""

import argparse
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.model_selection import KFold
from sklearn.linear_model import RidgeCV

from strata.data import gdsc

SEED = 1337
ALPHAS = np.logspace(-1, 4, 12)   # Ridge strength grid (p>>n -> expect strong reg)


def metrics(y_true, y_pred):
    rmse = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
    # Pearson undefined when the prediction is constant (e.g. global-mean baseline)
    if len(y_true) <= 2 or np.ptp(y_pred) == 0:
        return rmse, float("nan")
    return rmse, float(pearsonr(y_true, y_pred)[0])


def build_matrix(drug_name, n_genes_cap=None):
    meth = gdsc.load_gdsc_methylation()
    drug = gdsc.load_gdsc_drug_response()
    anno = gdsc.load_gdsc_annotations()

    d = drug[drug["drug_name"].str.lower() == drug_name.lower()]
    if d.empty:
        raise SystemExit(f"Drug '{drug_name}' not found.")
    auc = d.groupby("COSMIC_ID")["auc"].mean()

    common = meth.index.intersection(auc.index)
    M = meth.loc[common]
    y = auc.loc[common].to_numpy(dtype=np.float64)
    cancer = anno.reindex(common)["cancer_type"].fillna("Unknown").to_numpy()
    return M, y, cancer, list(M.columns)


def select_and_scale(Mtr, Mte, n_genes):
    """Top-variance genes on TRAIN; z-score with TRAIN stats. Returns (Xtr, Xte, idx)."""
    var = Mtr.var(axis=0)
    idx = np.argsort(var)[::-1][:n_genes]
    Xtr, Xte = Mtr[:, idx], Mte[:, idx]
    mu = Xtr.mean(axis=0)
    sd = Xtr.std(axis=0); sd[sd == 0] = 1.0
    return (Xtr - mu) / sd, (Xte - mu) / sd, idx


def onehot(cancer_tr, cancer_te):
    cats = sorted(pd.unique(cancer_tr))
    index = {c: i for i, c in enumerate(cats)}
    def enc(arr):
        Z = np.zeros((len(arr), len(cats)), dtype=np.float64)
        for i, c in enumerate(arr):
            j = index.get(c)
            if j is not None:
                Z[i, j] = 1.0
        return Z
    return enc(cancer_tr), enc(cancer_te)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drug", default="Palbociclib")
    ap.add_argument("--n-genes", type=int, default=2000)
    ap.add_argument("--folds", type=int, default=5)
    args = ap.parse_args()

    M, y, cancer, genes = build_matrix(args.drug)
    Mv = M.to_numpy(dtype=np.float64)
    print(f"=== {args.drug}: {len(y)} cell lines | {Mv.shape[1]} genes | "
          f"{len(set(cancer))} cancer types ===")
    print(f"AUC mean={y.mean():.3f} std={y.std():.3f}  (lower AUC = more sensitive)\n")

    kf = KFold(n_splits=args.folds, shuffle=True, random_state=SEED)
    acc = {k: {"rmse": [], "r": []} for k in
           ["global_mean", "cancer_mean", "meth_only", "meth_cancer"]}

    for tr, te in kf.split(Mv):
        ytr, yte = y[tr], y[te]

        # 1. global mean
        acc["global_mean"]["rmse"].append(metrics(yte, np.full_like(yte, ytr.mean()))[0])
        acc["global_mean"]["r"].append(float("nan"))

        # 2. cancer-type mean (fallback to global for unseen types)
        s = pd.Series(ytr).groupby(pd.Series(cancer[tr])).mean()
        pred_ct = np.array([s.get(c, ytr.mean()) for c in cancer[te]])
        r = metrics(yte, pred_ct); acc["cancer_mean"]["rmse"].append(r[0]); acc["cancer_mean"]["r"].append(r[1])

        # features
        Xtr, Xte, _ = select_and_scale(Mv[tr], Mv[te], args.n_genes)
        Ctr, Cte = onehot(cancer[tr], cancer[te])

        # 3. methylation only
        m3 = RidgeCV(alphas=ALPHAS).fit(Xtr, ytr)
        r = metrics(yte, m3.predict(Xte)); acc["meth_only"]["rmse"].append(r[0]); acc["meth_only"]["r"].append(r[1])

        # 4. methylation + cancer one-hot
        m4 = RidgeCV(alphas=ALPHAS).fit(np.hstack([Xtr, Ctr]), ytr)
        r = metrics(yte, m4.predict(np.hstack([Xte, Cte]))); acc["meth_cancer"]["rmse"].append(r[0]); acc["meth_cancer"]["r"].append(r[1])

    print(f"{args.folds}-fold cell-line-blind CV (mean +/- std):")
    print(f"{'model':<22}{'RMSE':>16}{'Pearson r':>16}")
    for k, name in [("global_mean", "global mean"), ("cancer_mean", "cancer-type mean"),
                    ("meth_only", "methylation"), ("meth_cancer", "methylation+cancer")]:
        rm = np.array(acc[k]["rmse"]); rr = np.array(acc[k]["r"])
        rstr = "n/a" if np.isnan(rr).all() else f"{np.nanmean(rr):.3f}+/-{np.nanstd(rr):.3f}"
        print(f"{name:<22}{np.mean(rm):.4f}+/-{np.std(rm):.4f}{rstr:>16}")

    base = np.mean(acc["cancer_mean"]["rmse"])
    best = np.mean(acc["meth_cancer"]["rmse"])
    print(f"\nmethylation+cancer vs cancer-type baseline: "
          f"{(base-best)/base*100:+.1f}% RMSE "
          f"({'BEATS' if best < base else 'does NOT beat'} the tissue bar)")

    # ---- interpretability: refit on all data, top genes by |coef| ----
    Xall, _, idx = select_and_scale(Mv, Mv, args.n_genes)
    Call, _ = onehot(cancer, cancer)
    full = RidgeCV(alphas=ALPHAS).fit(np.hstack([Xall, Call]), y)
    coef_genes = full.coef_[:len(idx)]
    order = np.argsort(np.abs(coef_genes))[::-1][:15]
    print(f"\nTop 15 methylation drivers (Ridge coef; negative = methylation->sensitivity):")
    for j in order:
        print(f"  {genes[idx[j]]:<12} coef={coef_genes[j]:+.4f}")


if __name__ == "__main__":
    main()
