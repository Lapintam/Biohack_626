"""Global-axis-controlled ALL-GENE discovery scan (unbiased, beyond known genes).

R1b showed a naive all-gene tissue-adjusted scan is non-specific (~8k/14.6k genes
"significant") because the methylome is pervasively co-correlated. Fix (SVA/RUV-style
latent-factor adjustment): residualize BOTH the drug AUC and every gene's methylation
against [tissue fixed-effects + the top-K principal components of the methylome] — the
shared "global axis" — then rank ALL genes by the partial correlation. The global axis
that made everything co-correlate is removed, leaving gene-specific signal.

Run: PYTHONPATH=. uv run python scripts/v2_global_axis_scan.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA

from strata.data import gdsc

MIN_TISSUE_N = 20


def _resid(Y, Z):
    """Residuals of Y on design Z (FWL). Y: (n,) or (n,g)."""
    B, *_ = np.linalg.lstsq(Z, Y, rcond=None)
    return Y - Z @ B


def scan(meth, auc, tissue, control="tissue+pcs", n_pcs=3):
    """Return ranked DataFrame[gene, r, pvalue, qvalue] after adjusting for `control`.
    control: 'none' (pooled), 'tissue', or 'tissue+pcs'."""
    common = meth.index.intersection(auc.index).intersection(tissue.index)
    t = tissue.loc[common]
    keep = t.map(t.value_counts()) >= MIN_TISSUE_N
    common = common[keep.values]
    Mdf = meth.loc[common]
    # coverage filter (promoter RRBS has NaN); then variance filter; then mean-fill
    cov = Mdf.notna().mean(axis=0)
    Mdf = Mdf.loc[:, cov >= 0.60]
    gstd = Mdf.std(axis=0)
    Mdf = Mdf.loc[:, gstd >= 0.02]
    Mdf = Mdf.fillna(Mdf.mean(axis=0))   # per-gene mean-fill (no-op for dense gene-level)
    genes = Mdf.columns
    M = Mdf.values.astype(float)
    a = auc.loc[common].values.astype(float)
    t = tissue.loc[common]
    n, g = M.shape

    cols = [np.ones(n)]
    if control in ("tissue", "tissue+pcs"):
        cols.append(pd.get_dummies(t, drop_first=True).values.astype(float))
    if control == "tissue+pcs" and n_pcs > 0:
        pcs = PCA(n_components=n_pcs, random_state=0).fit_transform(M - M.mean(0))
        cols.append(pcs)
    Z = np.column_stack(cols)

    a_r = _resid(a, Z); M_r = _resid(M, Z)
    a_r = a_r - a_r.mean(); M_r = M_r - M_r.mean(0)
    num = M_r.T @ a_r
    den = np.linalg.norm(M_r, axis=0) * np.linalg.norm(a_r)
    den[den == 0] = np.nan
    r = num / den
    df = n - Z.shape[1] - 1
    tstat = r * np.sqrt(df / (1 - r**2))
    pv = 2 * stats.t.sf(np.abs(tstat), df)
    out = pd.DataFrame({"gene": genes, "r": r, "pvalue": pv}).dropna()
    out["qvalue"] = multipletests(out["pvalue"], method="fdr_bh")[1]
    out = out.reindex(out["r"].abs().sort_values(ascending=False).index).reset_index(drop=True)
    return out, n


def _auc(drug, name):
    s = (drug[drug["drug_name"].str.lower() == name.lower()]
         .drop_duplicates("COSMIC_ID").set_index("COSMIC_ID")["auc"])
    return s[~s.index.duplicated(keep="first")].astype(float)


def report(meth, drug, tissue, name, known=()):
    print("#" * 80)
    print(f"# {name}  — naive vs tissue vs tissue+PCs (global-axis control)")
    print("#" * 80)
    auc = _auc(drug, name)
    for ctrl in ("none", "tissue", "tissue+pcs"):
        res, n = scan(meth, auc, tissue, control=ctrl)
        nsig = int((res["qvalue"] < 0.1).sum())
        top = res.head(8)
        tops = ", ".join(f"{r.gene}({r.r:+.2f})" for r in top.itertuples())
        print(f"  [{ctrl:<11}] n={n}  q<0.1: {nsig:>5}/{len(res)}   top: {tops}")
        for kg in known:
            hit = res.index[res['gene'] == kg]
            if len(hit):
                rk = int(hit[0]) + 1
                print(f"        {kg}: rank {rk}/{len(res)}  r={res.loc[hit[0],'r']:+.3f} q={res.loc[hit[0],'qvalue']:.2g}")
    print()


def main():
    import sys
    from strata.data import ccle_methylation as ccle
    drug = gdsc.load_gdsc_drug_response()
    tissue = gdsc.load_gdsc_annotations()["tissue"]
    src = sys.argv[1] if len(sys.argv) > 1 else "promoter"
    meth = ccle.load_ccle_promoter_methylation() if src == "promoter" else gdsc.load_gdsc_methylation()
    print(f">>> SOURCE: {src}  ({meth.shape[0]} lines x {meth.shape[1]} genes)\n")
    report(meth, drug, tissue, "Palbociclib", known=("MTAP", "CDKN2A"))
    report(meth, drug, tissue, "Trametinib", known=("BRAF", "KRAS", "DUSP6"))
    report(meth, drug, tissue, "Nutlin-3a (-)", known=("TP53", "MDM2", "CDKN2A"))
    report(meth, drug, tissue, "Temozolomide", known=("MGMT",))


if __name__ == "__main__":
    main()
