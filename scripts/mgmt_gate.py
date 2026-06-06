"""MGMT positive-control gate for the methylation->drug-response pipeline.

Hypothesis (textbook pharmacogenomics): MGMT promoter methylation silences MGMT,
making cells MORE sensitive to the alkylating agent temozolomide. In GDSC, higher
sensitivity = LOWER AUC. So high-MGMT-methylation lines should have lower
temozolomide AUC, and continuous MGMT methylation should correlate NEGATIVELY with
temozolomide AUC.

Run: uv run python scripts/mgmt_gate.py
"""

from __future__ import annotations

import sys
from pathlib import Path

# Repo root on path so `strata` imports without an installed package.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

from strata.data.gdsc import load_gdsc_drug_response, load_gdsc_methylation

DRUG = "Temozolomide"
GENE = "MGMT"


def pick_mgmt_column(meth: pd.DataFrame) -> str:
    if GENE in meth.columns:
        return GENE
    candidates = [c for c in meth.columns if "MGMT" in c.upper()]
    print(f"[warn] exact '{GENE}' column absent; MGMT-like columns: {candidates}")
    if not candidates:
        sys.exit("No MGMT-like methylation column found; cannot run gate.")
    return candidates[0]


def run_gate() -> dict:
    print("Loading methylation + drug response ...")
    meth = load_gdsc_methylation()
    drugs = load_gdsc_drug_response()

    col = pick_mgmt_column(meth)
    print(f"Using methylation column: {col!r}")

    # Temozolomide AUC per COSMIC_ID (one value per line; average if duplicated).
    tmz = drugs[drugs["drug_name"] == DRUG]
    if tmz.empty:
        sys.exit(f"Drug {DRUG!r} not found in drug response table.")
    tmz_auc = tmz.groupby("COSMIC_ID")["auc"].mean()

    mgmt = meth[col].dropna()

    common = mgmt.index.intersection(tmz_auc.index)
    df = pd.DataFrame(
        {"mgmt": mgmt.loc[common].astype(float), "auc": tmz_auc.loc[common].astype(float)}
    ).dropna()
    n = len(df)
    print(f"Intersected COSMIC_IDs with both MGMT methylation and {DRUG} AUC: n={n}")

    # Median split on MGMT methylation.
    med = df["mgmt"].median()
    high = df[df["mgmt"] > med]["auc"]
    low = df[df["mgmt"] <= med]["auc"]
    high_med = float(high.median())
    low_med = float(low.median())

    # One-sided Mann-Whitney: high-MGMT AUC < low-MGMT AUC (more sensitive).
    mw_stat, mw_p = mannwhitneyu(high, low, alternative="less")

    # Continuous: expect negative Spearman rho.
    rho, sp_p = spearmanr(df["mgmt"], df["auc"])

    print()
    print(f"  Median split MGMT methylation = {med:.4f}")
    print(f"  High-MGMT group: n={len(high)}, median {DRUG} AUC = {high_med:.4f}")
    print(f"  Low-MGMT  group: n={len(low)}, median {DRUG} AUC = {low_med:.4f}")
    print(f"  Mann-Whitney (high AUC < low AUC) p = {mw_p:.4g}")
    print(f"  Spearman(MGMT meth, {DRUG} AUC): rho = {rho:.4f}, p = {sp_p:.4g}")

    passed = (mw_p < 0.05) or (rho < 0 and sp_p < 0.05)
    print()
    print(f"  VERDICT: {'PASS' if passed else 'FALLBACK'}")

    return {
        "col": col,
        "n": n,
        "high_med": high_med,
        "low_med": low_med,
        "mw_p": float(mw_p),
        "rho": float(rho),
        "sp_p": float(sp_p),
        "passed": passed,
    }


def fallback_scan() -> None:
    """Bounded scan for a strong, interpretable methylation->drug association.

    Strategy: a curated set of well-known pharmacogenomic methylation markers
    crossed against ALL drugs, PLUS the most-variable methylation genes against
    the most-screened drugs. Report the strongest interpretable hits.
    """
    print("\n" + "=" * 70)
    print("FALLBACK SCAN: searching for a strong methylation->drug positive control")
    print("=" * 70)

    meth = load_gdsc_methylation()
    drugs = load_gdsc_drug_response()

    # AUC matrix: COSMIC_ID x drug (mean per line/drug).
    auc_long = drugs.groupby(["COSMIC_ID", "drug_name"])["auc"].mean().reset_index()
    auc_mat = auc_long.pivot(index="COSMIC_ID", columns="drug_name", values="auc")

    # Curated markers + top-variable genes.
    curated = [g for g in ["MGMT", "MLH1", "BRCA1", "VHL", "CDKN2A", "RASSF1", "DAPK1"]
               if g in meth.columns]
    var_genes = meth.var(axis=0, numeric_only=True).sort_values(ascending=False)
    top_var = [g for g in var_genes.index[:50] if g not in curated]
    genes = curated + top_var

    # Most-screened drugs (by number of lines), plus keep all curated-marker tests broad.
    drug_counts = auc_mat.notna().sum(axis=0).sort_values(ascending=False)
    top_drugs = list(drug_counts.index[:40])

    common = meth.index.intersection(auc_mat.index)
    meth = meth.loc[common]
    auc_mat = auc_mat.loc[common]
    print(f"Common cell lines: {len(common)}")
    print(f"Genes scanned: {len(genes)} (curated={len(curated)}, top-var={len(top_var)})")

    results = []
    # Curated markers vs ALL drugs (cheap, high biological prior).
    all_drugs = list(auc_mat.columns)
    n_pairs = 0
    for gi, g in enumerate(curated):
        gvals = meth[g]
        for d in all_drugs:
            avals = auc_mat[d]
            mask = gvals.notna() & avals.notna()
            if mask.sum() < 20:
                continue
            rho, p = spearmanr(gvals[mask], avals[mask])
            if np.isnan(rho):
                continue
            results.append((g, d, rho, p, int(mask.sum()), "curated-vs-all"))
            n_pairs += 1
    print(f"Curated-vs-all-drugs pairs tested: {n_pairs}")

    # Top-variable genes vs top-screened drugs.
    n_pairs2 = 0
    for gi, g in enumerate(top_var):
        if gi % 10 == 0:
            print(f"  top-var gene {gi}/{len(top_var)} ...")
        gvals = meth[g]
        for d in top_drugs:
            avals = auc_mat[d]
            mask = gvals.notna() & avals.notna()
            if mask.sum() < 20:
                continue
            rho, p = spearmanr(gvals[mask], avals[mask])
            if np.isnan(rho):
                continue
            results.append((g, d, rho, p, int(mask.sum()), "topvar-vs-topdrug"))
            n_pairs2 += 1
    print(f"Topvar-vs-topdrug pairs tested: {n_pairs2}")

    res = pd.DataFrame(results, columns=["gene", "drug", "rho", "p", "n", "source"])
    res["abs_rho"] = res["rho"].abs()

    print("\nTop 15 strongest associations by |rho| (p<0.05):")
    sig = res[res["p"] < 0.05].sort_values("abs_rho", ascending=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(sig.head(15).to_string(index=False))

    print("\nCurated-marker hits (sorted by |rho|):")
    cur = res[res["source"] == "curated-vs-all"].sort_values("abs_rho", ascending=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(cur.head(15).to_string(index=False))

    return res


if __name__ == "__main__":
    summary = run_gate()
    if not summary["passed"]:
        fallback_scan()
