"""
Polymer Benchmark — a fixed, model-agnostic test bed for "methylome -> drug response".

The point: the naive structure of this task (some drugs are potent on *everyone*,
e.g. Sepantronium) makes raw R^2 look high without any personalization. This
benchmark exposes a SECOND, much harder score -- the PERSONALIZED R^2, computed on
each drug's mean-removed residuals -- where the per-drug-mean baseline scores 0 by
construction. That residual is the real signal; driving it up is the game.

Artifacts written to outputs/benchmark/:
    manifest.json          - seed, fractions, drug list, split sizes, the contract
    splits.json            - train / val / test COSMIC_IDs (test = never train on these)
    test_methylomes.csv.gz - INPUT: test lines x gene-symbol betas (the methylome inputs)
    truth_test.csv.gz      - ANSWER KEY: (COSMIC_ID, drug_name, true_auc, cancer_type)
    train_responses.csv.gz - (COSMIC_ID, drug_name, auc) for train+val (fit any model on this)

Usage:
    python polymer_benchmark.py build         # create the benchmark files
    python polymer_benchmark.py leaderboard   # score reference models + Polymer

To score your own model, produce a DataFrame with columns
[COSMIC_ID, drug_name, pred_auc] over the test pairs and call:
    from polymer_benchmark import score, load_truth
    print(score(my_pred_df))
"""

import os, json, sys
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from strata.data import gdsc
from strata.config import OUTPUT_DIR

SEED = 1337
VAL_FRAC = TEST_FRAC = 0.15
MIN_DRUGS_PER_LINE = 10     # lines need >=10 measured drugs to score a ranking

SOURCES = {
    "imputed": gdsc.load_gdsc_methylation,
    "genelevel": gdsc.load_gdsc_methylation_genelevel,
    "promoter": gdsc.load_gdsc_methylation_promoter,
}


def bench_dir(source="imputed"):
    sub = "benchmark" if source == "imputed" else f"benchmark_{source}"
    return os.path.join(OUTPUT_DIR, sub)


def artifact_path(source="imputed"):
    suffix = "" if source == "imputed" else f"_{source}"
    return os.path.join(OUTPUT_DIR, f"polymer_model{suffix}.pt")


# ============================================================ build ==========
def _split_lines(meth_index):
    """Reproduce EXACTLY the recommender's cell-line split (seed 1337, 15/15)."""
    rng = np.random.default_rng(SEED)
    lines = np.asarray(meth_index)
    n = len(lines)
    order = rng.permutation(n)
    n_te, n_va = int(n * TEST_FRAC), int(n * VAL_FRAC)
    test = set(lines[order[:n_te]]); val = set(lines[order[n_te:n_te + n_va]])
    train = [c for c in lines if c not in test and c not in val]
    return train, sorted(val), sorted(test)


def build(source="imputed"):
    BENCH_DIR = bench_dir(source)
    os.makedirs(BENCH_DIR, exist_ok=True)
    meth = SOURCES[source]()
    meth.index = meth.index.astype(str)
    drug = gdsc.load_gdsc_drug_response()
    anno = gdsc.load_gdsc_annotations()

    train, val, test = _split_lines(meth.index)
    drug = drug.dropna(subset=["auc"])
    drug = drug[drug["COSMIC_ID"].astype(str).isin(set(meth.index))]
    drug = drug.groupby(["COSMIC_ID", "drug_name"], as_index=False)["auc"].mean()
    drug["cancer_type"] = drug["COSMIC_ID"].astype(str).map(anno["cancer_type"])

    test_set = set(test)
    truth = drug[drug["COSMIC_ID"].astype(str).isin(test_set)].rename(columns={"auc": "true_auc"})
    train_resp = drug[~drug["COSMIC_ID"].astype(str).isin(test_set)][["COSMIC_ID", "drug_name", "auc"]]

    # methylome inputs for the test lines (full gene-symbol beta matrix)
    test_meth = meth.loc[[c for c in test if c in meth.index]]

    truth.to_csv(os.path.join(BENCH_DIR, "truth_test.csv.gz"), index=False, compression="gzip")
    train_resp.to_csv(os.path.join(BENCH_DIR, "train_responses.csv.gz"), index=False, compression="gzip")
    test_meth.to_csv(os.path.join(BENCH_DIR, "test_methylomes.csv.gz"), compression="gzip")
    json.dump({"train": list(train), "val": list(val), "test": list(test)},
              open(os.path.join(BENCH_DIR, "splits.json"), "w"))
    manifest = {
        "name": "Polymer-GDSC methylome->drug-response benchmark",
        "seed": SEED, "val_frac": VAL_FRAC, "test_frac": TEST_FRAC,
        "n_train_lines": len(train), "n_val_lines": len(val), "n_test_lines": len(test),
        "n_test_pairs": int(len(truth)), "n_drugs": int(drug["drug_name"].nunique()),
        "n_genes": int(meth.shape[1]),
        "contract": "Models MUST NOT train on test COSMIC_IDs. Fit on train_responses + "
                    "any methylation for train/val lines; predict test pairs.",
        "headline_metric": "personalized_r2 (R^2 on per-drug-mean residuals; "
                           "per-drug-mean baseline = 0).",
    }
    json.dump(manifest, open(os.path.join(BENCH_DIR, "manifest.json"), "w"), indent=2)
    print(f"built benchmark in {BENCH_DIR}")
    for k, v in manifest.items():
        if k not in ("contract", "headline_metric"):
            print(f"  {k}: {v}")
    return manifest


# ============================================================ scoring ========
def load_truth(source="imputed"):
    BENCH_DIR = bench_dir(source)
    truth = pd.read_csv(os.path.join(BENCH_DIR, "truth_test.csv.gz"))
    train_resp = pd.read_csv(os.path.join(BENCH_DIR, "train_responses.csv.gz"))
    drug_mean = train_resp.groupby("drug_name")["auc"].mean()   # the "everyone" prior
    return truth, drug_mean


def _r2(y, p):
    ss_res = np.sum((y - p) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")


def _safe_pearson(y, p):
    if len(y) < 2 or np.ptp(p) == 0 or np.ptp(y) == 0:
        return float("nan")
    return float(pearsonr(y, p)[0])


def _safe_spearman(a, b):
    if len(a) < 2 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return np.nan
    return spearmanr(a, b).correlation


def score(pred_df, truth=None, drug_mean=None):
    """pred_df: columns [COSMIC_ID, drug_name, pred_auc]. Returns a metrics dict."""
    if truth is None:
        truth, drug_mean = load_truth()
    truth = truth.copy(); truth["COSMIC_ID"] = truth["COSMIC_ID"].astype(str)
    pred_df = pred_df.copy(); pred_df["COSMIC_ID"] = pred_df["COSMIC_ID"].astype(str)
    df = truth.merge(pred_df, on=["COSMIC_ID", "drug_name"], how="inner")
    if len(df) < len(truth):
        print(f"  [warn] {len(truth)-len(df)} truth pairs missing from predictions")
    if len(df) < 2:
        raise ValueError("No overlapping (COSMIC_ID, drug_name) pairs between truth and predictions.")
    y = df["true_auc"].to_numpy(float)
    p = df["pred_auc"].to_numpy(float)

    # global metrics (the "easy" R^2 -- dominated by per-drug potency)
    rmse = float(np.sqrt(np.mean((y - p) ** 2)))
    r = _safe_pearson(y, p)
    r2 = _r2(y, p)

    # PERSONALIZED metrics: remove each drug's training mean from BOTH sides
    dm = df["drug_name"].map(drug_mean).to_numpy(float)
    yr, pr = y - dm, p - dm
    pers_r2 = _r2(yr, pr)                       # <-- the honest, much-lower R^2
    pers_r = _safe_pearson(yr, pr)

    # within-line ranking + top-5 recovery
    raw_rho, dem_rho, hits5 = [], [], []
    for cid, g in df.groupby("COSMIC_ID"):
        if len(g) < MIN_DRUGS_PER_LINE:
            continue
        yy, pp = g["true_auc"].to_numpy(float), g["pred_auc"].to_numpy(float)
        raw_rho.append(_safe_spearman(pp, yy))
        d = g["drug_name"].map(drug_mean).to_numpy(float)
        dem_rho.append(_safe_spearman(pp - d, yy - d))
        if len(g) >= 11:
            true_top = set(g.sort_values("true_auc")["drug_name"].head(10))
            pred_top = set(g.sort_values("pred_auc")["drug_name"].head(5))
            hits5.append(len(true_top & pred_top))

    def _mean(a):
        a = [x for x in a if x == x]      # drop NaNs
        return float(np.mean(a)) if a else float("nan")

    return {
        "n_pairs": len(df),
        "RMSE": rmse, "Pearson_r": r, "R2": r2,
        "personalized_R2": pers_r2, "personalized_r": pers_r,
        "within_line_rho_raw": _mean(raw_rho),
        "within_line_rho_demeaned": _mean(dem_rho),
        "top5_recovery": _mean(hits5),
    }


# ====================================================== reference models =====
def _predict_global_mean(truth, train_resp):
    g = train_resp["auc"].mean()
    return truth.assign(pred_auc=g)[["COSMIC_ID", "drug_name", "pred_auc"]]


def _predict_per_drug_mean(truth, train_resp):
    dm = train_resp.groupby("drug_name")["auc"].mean()
    glob = train_resp["auc"].mean()
    pred = truth["drug_name"].map(dm).fillna(glob)
    return truth.assign(pred_auc=pred.values)[["COSMIC_ID", "drug_name", "pred_auc"]]


def _predict_per_drug_tissue(truth, train_resp, anno):
    ct_map = anno["cancer_type"]
    ct_map.index = ct_map.index.astype(str)
    tr = train_resp.copy()
    tr["cancer_type"] = tr["COSMIC_ID"].astype(str).map(ct_map)
    dt = tr.groupby(["drug_name", "cancer_type"])["auc"].mean()
    dm = tr.groupby("drug_name")["auc"].mean()
    glob = tr["auc"].mean()
    ct = truth["COSMIC_ID"].astype(str).map(ct_map)
    pred = [dt.get((d, c), dm.get(d, glob))
            for d, c in zip(truth["drug_name"], ct)]
    return truth.assign(pred_auc=pred)[["COSMIC_ID", "drug_name", "pred_auc"]]


def _predict_polymer(truth, source="imputed"):
    from polymer_infer import Polymer
    meth = SOURCES[source]()
    meth.index = meth.index.astype(str)
    poly = Polymer(path=artifact_path(source))
    out = []
    meth_index = set(meth.index.astype(str))
    for cid, g in truth.groupby("COSMIC_ID"):
        cid_s = str(cid)
        if cid_s not in meth_index:
            continue
        betas = meth.loc[cid_s].to_dict()
        rows, _ = poly.recommend(betas)              # all drugs
        pred = {r["drug"]: r["pred_auc"] for r in rows}
        for dn in g["drug_name"]:
            if dn in pred:
                out.append((cid_s, dn, pred[dn]))
    return pd.DataFrame(out, columns=["COSMIC_ID", "drug_name", "pred_auc"])


def leaderboard(source="imputed"):
    BENCH_DIR = bench_dir(source)
    truth, drug_mean = load_truth(source)
    train_resp = pd.read_csv(os.path.join(BENCH_DIR, "train_responses.csv.gz"))
    anno = gdsc.load_gdsc_annotations()

    models = {
        "GlobalMean": _predict_global_mean(truth, train_resp),
        "PerDrugMean": _predict_per_drug_mean(truth, train_resp),
        "PerDrug+Tissue": _predict_per_drug_tissue(truth, train_resp, anno),
        "Polymer (NN)": _predict_polymer(truth, source),
    }
    rows = []
    for name, pred in models.items():
        s = score(pred, truth, drug_mean)
        s["model"] = name
        rows.append(s)
    lb = pd.DataFrame(rows).set_index("model")
    cols = ["RMSE", "R2", "Pearson_r", "personalized_R2", "personalized_r",
            "within_line_rho_demeaned", "top5_recovery", "n_pairs"]
    lb = lb[cols].sort_values("personalized_R2", ascending=False)
    out = os.path.join(BENCH_DIR, "leaderboard.csv")
    lb.to_csv(out)
    pd.set_option("display.width", 160)
    print(f"\n============== POLYMER BENCHMARK LEADERBOARD [source={source}] ==============")
    print("(ranked by personalized_R2 — the hard, mean-removed signal)\n")
    print(lb.round(4).to_string())
    print(f"\nNote: 'R2' is the easy global score (per-drug potency dominates).")
    print(f"      'personalized_R2' strips each drug's mean — PerDrugMean = 0 there by design.")
    print(f"saved {out}")
    return lb


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", nargs="?", default="build", choices=["build", "leaderboard"])
    ap.add_argument("--source", default="imputed", choices=list(SOURCES))
    args = ap.parse_args()
    if args.cmd == "build":
        build(args.source)
    elif args.cmd == "leaderboard":
        if not os.path.exists(os.path.join(bench_dir(args.source), "truth_test.csv.gz")):
            build(args.source)
        leaderboard(args.source)
