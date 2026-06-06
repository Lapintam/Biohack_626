"""
Associate gene methylation (Beta) with drug response (AUC) for a given drug,
with cancer type as an added dimension.

This builds the tidy table the analysis needs -- one row per cell line:

    COSMIC_ID | Beta(gene) | Drug_Name | AUC | Cancer_Type

...and answers two questions:

  1. Which genes' methylation correlates with AUC for this drug?
     (Spearman per gene, BH-FDR corrected -- strata.engine.rank_response_markers)
  2. Does the lead association survive adding cancer type?
     Tissue is a known confounder here, so we report the raw association AND a
     cancer-type-stratified view + a tissue-controlled partial correlation.

Lower AUC = more sensitive. A NEGATIVE rho means "more methylation -> lower AUC
-> more sensitive" (methylation marks a responder subgroup).

Usage:
    python methylation_auc_association.py --drug Palbociclib --gene CDKN2A
    python methylation_auc_association.py --drug Temozolomide          # top markers only
"""

import argparse
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

from strata.data import gdsc
from strata.engine.associate import rank_response_markers, discover_responder_subgroup


def build_table(drug_name, gene=None):
    """Tidy per-cell-line table: Beta(gene) | Drug_Name | AUC | Cancer_Type.
    If `gene` is None, Beta is omitted (use the full matrix for marker ranking).
    Returns (table, meth_aligned, auc_series, cancer_series)."""
    meth = gdsc.load_gdsc_methylation()
    drug = gdsc.load_gdsc_drug_response()
    anno = gdsc.load_gdsc_annotations()

    d = drug[drug["drug_name"].str.lower() == drug_name.lower()].copy()
    if d.empty:
        raise SystemExit(f"Drug '{drug_name}' not found. "
                         f"e.g. {sorted(drug['drug_name'].unique())[:8]}")
    # one AUC per cell line for this drug (collapse rare dup curves by mean)
    auc = d.groupby("COSMIC_ID")["auc"].mean()

    common = meth.index.intersection(auc.index)
    meth_a = meth.loc[common]
    auc = auc.loc[common]
    cancer = anno.reindex(common)["cancer_type"]

    table = pd.DataFrame({
        "COSMIC_ID": common,
        "Drug_Name": drug_name,
        "AUC": auc.values,
        "Cancer_Type": cancer.values,
    })
    if gene is not None:
        if gene not in meth.columns:
            raise SystemExit(f"Gene '{gene}' not in methylation matrix.")
        table.insert(1, "Beta", meth_a[gene].values)
    return table, meth_a, auc, cancer


def partial_spearman_controlling_tissue(beta, auc, cancer):
    """Spearman(beta, auc) after removing the cancer-type (tissue) main effect.
    Rank-transform both, subtract per-cancer-type rank means (group-mean center),
    then correlate residuals -> a within-cancer-type association."""
    df = pd.DataFrame({"b": beta, "a": auc, "c": cancer.values}).dropna()
    rb = rankdata(df["b"]); ra = rankdata(df["a"])
    df = df.assign(rb=rb, ra=ra)
    rb_res = df["rb"] - df.groupby("c")["rb"].transform("mean")
    ra_res = df["ra"] - df.groupby("c")["ra"].transform("mean")
    rho, p = spearmanr(rb_res, ra_res)
    return float(rho), float(p), len(df)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drug", default="Palbociclib")
    ap.add_argument("--gene", default=None, help="focus gene (e.g. CDKN2A)")
    ap.add_argument("--top", type=int, default=15)
    args = ap.parse_args()

    table, meth_a, auc, cancer = build_table(args.drug, args.gene)
    print(f"\n=== {args.drug}: {len(table)} cell lines with methylation + AUC ===")
    print(f"AUC  mean={table.AUC.mean():.3f}  std={table.AUC.std():.3f}")
    print(f"cancer types: {cancer.nunique()}  (top: "
          f"{', '.join(cancer.value_counts().head(3).index)})")

    # ---- 1. genome-wide marker ranking: methylation vs AUC ----
    print(f"\n--- Top {args.top} methylation markers of {args.drug} response "
          f"(Spearman beta vs AUC, q<0.1) ---")
    markers = rank_response_markers(meth_a, auc, min_n=30)
    if markers.empty:
        print("  (no genes pass q<0.1)")
    else:
        with pd.option_context("display.width", 120):
            print(markers.head(args.top).to_string(
                index=False,
                formatters={"rho": "{:+.3f}".format, "pvalue": "{:.2e}".format,
                            "qvalue": "{:.2e}".format}))

    # ---- 2. focus gene + cancer-type control ----
    if args.gene is not None:
        beta = table["Beta"].to_numpy()
        a = table["AUC"].to_numpy()
        rho, p = spearmanr(beta, a)
        print(f"\n--- Focus: {args.gene} methylation vs {args.drug} AUC ---")
        print(f"raw Spearman rho={rho:+.3f}  p={p:.2e}  n={len(beta)}  "
              f"({'sensitivity' if rho < 0 else 'resistance'} marker)")

        prho, pp, pn = partial_spearman_controlling_tissue(beta, a, cancer)
        print(f"tissue-controlled (within cancer type) rho={prho:+.3f}  p={pp:.2e}  n={pn}")

        # median-split responder subgroup (reuses the engine)
        clusters, assoc = discover_responder_subgroup(meth_a, auc, args.gene)
        print(f"median-split subgroup: effect(AUC drop)={assoc.effect_size:+.3f}  "
              f"p={assoc.pvalue:.2e}  grade={assoc.grade}  "
              f"(responder = {'high' if assoc.responder_label==1 else 'low'} methylation)")

        # per-cancer-type breakdown of the focus gene association
        print(f"\n--- {args.gene} vs AUC within the largest cancer types ---")
        tmp = table.dropna(subset=["Beta", "AUC", "Cancer_Type"])
        rows = []
        for ct, grp in tmp.groupby("Cancer_Type"):
            if len(grp) >= 15:
                r, pv = spearmanr(grp["Beta"], grp["AUC"])
                rows.append((ct, len(grp), r, pv))
        ct_df = (pd.DataFrame(rows, columns=["Cancer_Type", "n", "rho", "p"])
                 .sort_values("n", ascending=False).head(8))
        print(ct_df.to_string(index=False,
              formatters={"rho": "{:+.3f}".format, "p": "{:.2e}".format}))


if __name__ == "__main__":
    main()
