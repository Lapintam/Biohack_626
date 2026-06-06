"""Track B payoff — silencing-validation of promoter discovery candidates.

Pipeline per drug:
  1. global-axis-controlled promoter scan -> ~few-hundred response-associated
     candidate genes (tissue + methylome-PCs removed).
  2. SILENCING FILTER: keep only candidates whose promoter methylation tracks
     their OWN expression DOWN (Spearman(meth, expr) < 0, FDR-significant) =
     functional epigenetic silencing, not a passenger correlation.
  -> VALIDATED LEADS: genes that are (a) tissue-agnostically response-associated
     AND (b) functionally silenced by promoter methylation.

Run: PYTHONPATH=. uv run python scripts/v2_silencing_validation.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

from strata.data import gdsc
from strata.data import ccle_methylation as ccle
from strata.data import ccle_expression as cexpr
from scripts.v2_global_axis_scan import scan, _auc

MAX_CANDIDATES = 500   # cap response-candidates fed to the silencing test


def silencing(meth_col, expr_col):
    common = meth_col.dropna().index.intersection(expr_col.dropna().index)
    if len(common) < 50:
        return np.nan, np.nan, len(common)
    r, p = spearmanr(meth_col.loc[common], expr_col.loc[common])
    return float(r), float(p), len(common)


def run(prom, expr, drug, tissue, name):
    print("#" * 78)
    print(f"# {name}")
    print("#" * 78)
    auc = _auc(drug, name)
    res, n = scan(prom, auc, tissue, control="tissue+pcs")
    cand = res[res["qvalue"] < 0.1].head(MAX_CANDIDATES)
    rows = []
    for c in cand.itertuples():
        g = c.gene
        if g in expr.columns and g in prom.columns:
            sr, sp, sn = silencing(prom[g], expr[g])
            if not np.isnan(sr):
                rows.append((g, c.r, c.qvalue, sr, sp, sn))
    s = pd.DataFrame(rows, columns=["gene", "resp_r", "resp_q", "silencing_r", "silencing_p", "sil_n"])
    if len(s) == 0:
        print("  (no candidates with expression data)\n"); return None
    s["silencing_q"] = multipletests(s["silencing_p"], method="fdr_bh")[1]
    validated = s[(s["silencing_r"] < 0) & (s["silencing_q"] < 0.05)].copy()
    validated = validated.reindex(validated["resp_r"].abs().sort_values(ascending=False).index)
    print(f"  scan n={n} lines | response-candidates(q<0.1): {len(cand)} | "
          f"with expression: {len(s)} | FUNCTIONALLY SILENCED & validated: {len(validated)}")
    print(f"  {'gene':<12}{'resp_r':>8}{'resp_q':>9}{'sil_r':>8}{'sil_q':>9}   direction")
    for v in validated.head(15).itertuples():
        direction = "meth->SENSITIVE" if v.resp_r < 0 else "meth->resistant"
        print(f"  {v.gene:<12}{v.resp_r:>8.3f}{v.resp_q:>9.2g}{v.silencing_r:>8.3f}{v.silencing_q:>9.2g}   {direction}")
    print()
    return validated


def main():
    prom = ccle.load_ccle_promoter_methylation()
    expr = cexpr.load_ccle_expression()
    drug = gdsc.load_gdsc_drug_response()
    tissue = gdsc.load_gdsc_annotations()["tissue"]
    for dn in ("Palbociclib", "Trametinib", "Nutlin-3a (-)", "Temozolomide", "Olaparib"):
        run(prom, expr, drug, tissue, dn)


if __name__ == "__main__":
    main()
