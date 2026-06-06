"""v2 engine (R1b) — rank markers by TISSUE-ADJUSTED association, not pooled.

The pooled scan (`rank_response_markers`) surfaces lineage-confounded hits. The
taxonomy-crossing engine ranks each gene by its WITHIN-TISSUE (tissue-adjusted)
correlation with drug AUC: residualize both AUC and each gene's methylation on
tissue (subtract per-tissue means), then correlate the residuals. This is the
fixed-effects / §8a-safeguarded version of the discovery rule.

Two uses:
  1. validate: for Palbociclib, do MTAP/CDKN2A still rank high once tissue is removed?
  2. discover: scan many drugs, rank (drug, top tissue-adjusted marker) pairs to find
     a NOVEL tissue-agnostic cross-indication hit (DAPK1 was mostly lineage).

Run: PYTHONPATH=. uv run python scripts/v2_tissue_agnostic_scan.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from strata.data import gdsc

MIN_TISSUE_N = 15  # tissues smaller than this are dropped (within-tissue centering needs n)


def _auc(drug, drug_name):
    s = (drug[drug["drug_name"].str.lower() == drug_name.lower()]
         .drop_duplicates("COSMIC_ID").set_index("COSMIC_ID")["auc"])
    return s[~s.index.duplicated(keep="first")].astype(float)


def tissue_adjusted_scan(meth, auc, tissue) -> pd.DataFrame:
    """Return DataFrame[gene, r_adj, pvalue, qvalue, n, n_tissues] ranked by |r_adj|.
    r_adj = Pearson correlation of within-tissue-centered methylation vs within-tissue-
    centered AUC (= fixed-effects-on-tissue association)."""
    common = meth.index.intersection(auc.index).intersection(tissue.index)
    t = tissue.loc[common]
    keep_t = t.map(t.value_counts()) >= MIN_TISSUE_N
    common = common[keep_t.values]
    M = meth.loc[common]
    a = auc.loc[common]
    t = tissue.loc[common]
    n, n_tissues = len(common), t.nunique()

    # within-tissue center AUC and each gene (subtract per-tissue means)
    a_c = (a - a.groupby(t).transform("mean")).values
    tmeans = M.groupby(t.values).mean()            # n_tissues x genes
    M_c = M.values - tmeans.loc[t.values].values   # broadcast-subtract per-row tissue mean

    # vectorized Pearson of each centered gene column vs centered AUC
    a_c = a_c - a_c.mean()
    num = M_c.T @ a_c
    den = np.linalg.norm(M_c, axis=0) * np.linalg.norm(a_c)
    den[den == 0] = np.nan
    r = num / den

    df_resid = n - n_tissues - 1
    t_stat = r * np.sqrt(df_resid / (1 - r**2))
    pvals = 2 * stats.t.sf(np.abs(t_stat), df_resid)

    out = pd.DataFrame({"gene": M.columns, "r_adj": r, "pvalue": pvals})
    out = out.dropna(subset=["r_adj", "pvalue"])
    out["qvalue"] = multipletests(out["pvalue"], method="fdr_bh")[1]
    out["n"], out["n_tissues"] = n, n_tissues
    return out.reindex(out["r_adj"].abs().sort_values(ascending=False).index).reset_index(drop=True)


def _rank_of(df, gene):
    hit = df.index[df["gene"] == gene]
    return (int(hit[0]) + 1, df.loc[hit[0]]) if len(hit) else (None, None)


def main():
    meth = gdsc.load_gdsc_methylation()
    drug = gdsc.load_gdsc_drug_response()
    ann = gdsc.load_gdsc_annotations()
    tissue = ann["tissue"]

    # ---- 1) VALIDATE on Palbociclib ----
    print("#" * 74)
    print("# VALIDATE — tissue-adjusted scan: Palbociclib")
    print("#" * 74)
    pal = tissue_adjusted_scan(meth, _auc(drug, "Palbociclib"), tissue)
    print(f"n={pal['n'].iloc[0]} lines across {pal['n_tissues'].iloc[0]} tissues; "
          f"{int((pal['qvalue']<0.1).sum())} markers q<0.1")
    print("top 12 tissue-adjusted markers:")
    for _, row in pal.head(12).iterrows():
        print(f"   {row['gene']:<12} r_adj={row['r_adj']:+.3f}  q={row['qvalue']:.2g}")
    for g in ("MTAP", "CDKN2A"):
        rk, row = _rank_of(pal, g)
        print(f"   [{g}] tissue-adjusted rank {rk}/{len(pal)}  r_adj={row['r_adj']:+.3f} q={row['qvalue']:.2g}")

    # ---- 2) DISCOVER a novel tissue-agnostic hit across drugs ----
    print("\n" + "#" * 74)
    print("# DISCOVER — top tissue-adjusted marker per drug (novel cross-indication leads)")
    print("#" * 74)
    candidates = ["Trametinib", "Selumetinib", "Olaparib", "Talazoparib", "Dabrafenib",
                  "Alpelisib", "Dasatinib", "Nutlin-3a (-)", "Venetoclax", "Erlotinib",
                  "Afatinib", "Crizotinib", "Luminespib", "Vorinostat", "Nilotinib",
                  "AZD7762", "Sorafenib", "Lapatinib"]
    present = set(drug["drug_name"].unique())
    rows = []
    for dn in candidates:
        if dn not in present:
            continue
        sc = tissue_adjusted_scan(meth, _auc(drug, dn), tissue)
        if len(sc) == 0:
            continue
        top = sc.iloc[0]
        rows.append((dn, top["gene"], float(top["r_adj"]), float(top["qvalue"]),
                     int(top["n"]), int(top["n_tissues"]), int((sc["qvalue"] < 0.1).sum())))
    res = pd.DataFrame(rows, columns=["drug", "top_marker", "r_adj", "q", "n", "n_tissues", "n_sig"])
    res = res.reindex(res["r_adj"].abs().sort_values(ascending=False).index).reset_index(drop=True)
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
