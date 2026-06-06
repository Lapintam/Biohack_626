"""Held-out validation of the patient-annotation scorer (honest, no leakage).

Split cell lines into train/test. Rediscover leads on TRAIN only (global-axis scan
+ silencing). Build the per-line sensitivity score from train leads + the TRAIN
cohort distribution. Predict on TEST lines and correlate the predicted sensitivity
score against the MEASURED drug AUC on those held-out lines.
Expect NEGATIVE correlation (higher sensitivity score -> lower AUC = more sensitive).

Run: PYTHONPATH=. uv run python scripts/validate_annotation.py [DrugName]
"""

from __future__ import annotations

import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from strata.data import gdsc
from strata.data import ccle_methylation as ccle
from strata.data import ccle_expression as cexpr
from scripts.v2_global_axis_scan import scan
from scripts.v2_silencing_validation import silencing


def _auc(drug, name):
    s = (drug[drug["drug_name"].str.lower() == name.lower()]
         .drop_duplicates("COSMIC_ID").set_index("COSMIC_ID")["auc"])
    return s[~s.index.duplicated(keep="first")].astype(float)


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "Olaparib"
    prom = ccle.load_ccle_promoter_methylation()
    expr = cexpr.load_ccle_expression()
    drug = gdsc.load_gdsc_drug_response()
    tissue = gdsc.load_gdsc_annotations()["tissue"]
    auc = _auc(drug, name)

    lines = [c for c in auc.index if c in prom.index]
    test = [c for c in lines if int(c) % 10 < 3]      # ~30% deterministic held-out
    train = [c for c in lines if c not in set(test)]
    print(f"{name}: {len(lines)} lines -> train {len(train)} / test {len(test)}")

    # rediscover leads on TRAIN only
    res, n = scan(prom, auc.loc[train], tissue, control="tissue+pcs")
    cand = res[res["qvalue"] < 0.1].head(40)
    leads = []
    for c in cand.itertuples():
        g = c.gene
        if g in expr.columns and g in prom.columns:
            sr, sp, _ = silencing(prom[g], expr[g])
            if not np.isnan(sr) and sr < 0 and sp < 0.01:
                leads.append((g, -1.0 if c.r < 0 else 1.0, abs(c.r)))   # gene, dir_sign(sensitive=-r), weight
    print(f"  rediscovered {len(leads)} train leads")
    if len(leads) < 3:
        print("  too few leads to score"); return

    # TRAIN cohort distributions for percentile
    cohort = {g: np.sort(prom.loc[train, g].dropna().values) for g, _, _ in leads}

    def score_line(cosmic):
        contribs, wsum = 0.0, 0.0
        for g, dsign, w in leads:
            b = prom.loc[cosmic, g] if g in prom.columns else np.nan
            arr = cohort[g]
            if b == b and len(arr):                      # not NaN
                pct = np.searchsorted(arr, b) / len(arr)
                sens = (pct - 0.5) * 2.0 * (1 if dsign < 0 else -1)  # dsign<0 == sensitive lead
                contribs += sens * w; wsum += w
        return contribs / wsum if wsum else np.nan

    pred = pd.Series({c: score_line(c) for c in test}).dropna()
    meas = auc.loc[pred.index]
    rho, p = spearmanr(pred.values, meas.values)
    print(f"\n  HELD-OUT: predicted sensitivity score vs measured AUC")
    print(f"  Spearman rho = {rho:+.3f}  (p={p:.2g}, n={len(pred)} held-out lines)")
    print(f"  {'(negative = WORKS: higher predicted sensitivity -> lower AUC)' if rho < 0 else '(positive = does NOT validate)'}")
    # sensitivity of the top-decile predicted-sensitive vs rest
    hi = pred.sort_values(ascending=False)
    top = meas.loc[hi.head(max(5, len(hi)//5)).index].median()
    bot = meas.loc[hi.tail(max(5, len(hi)//5)).index].median()
    print(f"  median AUC: top-quintile predicted-sensitive={top:.3f} vs bottom-quintile={bot:.3f}  (lower=better)")


if __name__ == "__main__":
    main()
