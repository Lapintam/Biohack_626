"""Phase 1 payoff — external replication of silencing-validated leads.

For each demo drug: take the top silencing-validated leads (response-associated +
functionally silenced, from v2_silencing_validation), then re-test each in the
INDEPENDENT screens (PRISM, CTRP) by the same tissue-adjusted partial correlation.
A lead REPLICATES if the external association has the same sign as GDSC and p<0.05.
Replicated leads = bulletproof cross-indication cards (evidence ladder -> L5).

Run: PYTHONPATH=. uv run python scripts/v2_external_replication.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from strata.config import ROOT
from strata.data import gdsc
from strata.data import ccle_methylation as ccle
from strata.data import ccle_expression as cexpr
from strata.data import external_drug as ext
from scripts.v2_silencing_validation import run as silencing_run

TOP_LEADS = 25      # top validated leads per drug to test for replication
MIN_TISSUE = 15

# demo drug -> matched name in each screen (None = absent)
DRUGS = {
    "Palbociclib":  {"gdsc": "Palbociclib",    "prism": "palbociclib",  "ctrp": None},
    "Trametinib":   {"gdsc": "Trametinib",     "prism": "trametinib",   "ctrp": "trametinib"},
    "Temozolomide": {"gdsc": "Temozolomide",   "prism": "temozolomide", "ctrp": None},
    "Olaparib":     {"gdsc": "Olaparib",       "prism": "olaparib",     "ctrp": "olaparib"},
}


def _series(df, name):
    if name is None:
        return None
    s = (df[df["drug_name"].str.lower() == name.lower()]
         .drop_duplicates("COSMIC_ID").set_index("COSMIC_ID")["auc"])
    return s[~s.index.duplicated(keep="first")].astype(float)


def tissue_adj_corr(meth_g, auc, tissue):
    if auc is None:
        return np.nan, np.nan, 0
    d = pd.DataFrame({"m": meth_g, "a": auc, "t": tissue}).dropna()
    vc = d["t"].value_counts()
    d = d[d["t"].isin(vc[vc >= MIN_TISSUE].index)]
    if len(d) < 30 or d["m"].std() == 0:
        return np.nan, np.nan, len(d)
    Z = np.column_stack([np.ones(len(d)), pd.get_dummies(d["t"], drop_first=True).values.astype(float)])
    mr = d["m"].values - Z @ np.linalg.lstsq(Z, d["m"].values, rcond=None)[0]
    ar = d["a"].values - Z @ np.linalg.lstsq(Z, d["a"].values, rcond=None)[0]
    if np.std(mr) == 0 or np.std(ar) == 0:
        return np.nan, np.nan, len(d)
    r, p = pearsonr(mr, ar)
    return float(r), float(p), len(d)


def main():
    prom = ccle.load_ccle_promoter_methylation()
    expr = cexpr.load_ccle_expression()
    drug = gdsc.load_gdsc_drug_response()
    tissue = gdsc.load_gdsc_annotations()["tissue"]
    prism = ext.load_prism_drug_response()
    try:
        ctrp = ext.load_ctrp_drug_response()
    except Exception:
        ctrp = None

    all_rows = []
    for label, names in DRUGS.items():
        validated = silencing_run(prom, expr, drug, tissue, names["gdsc"])
        if validated is None or len(validated) == 0:
            continue
        leads = validated.head(TOP_LEADS)
        g_auc = _series(drug, names["gdsc"])
        p_auc = _series(prism, names["prism"])
        c_auc = _series(ctrp, names["ctrp"]) if ctrp is not None else None
        print(f"\n===== EXTERNAL REPLICATION: {label} =====")
        print(f"  {'gene':<13}{'gdsc_r':>8}{'prism_r':>9}{'prism_p':>9}{'ctrp_r':>8}{'ctrp_p':>9}  verdict")
        for v in leads.itertuples():
            g = v.gene
            if g not in prom.columns:
                continue
            gr, gp, gn = tissue_adj_corr(prom[g], g_auc, tissue)
            pr, pp, pn = tissue_adj_corr(prom[g], p_auc, tissue)
            cr, cp, cn = tissue_adj_corr(prom[g], c_auc, tissue)
            prism_ok = (not np.isnan(pr)) and (np.sign(pr) == np.sign(gr)) and pp < 0.05
            ctrp_ok = (not np.isnan(cr)) and (np.sign(cr) == np.sign(gr)) and cp < 0.05
            verdict = "REPLICATES" if prism_ok else ("ctrp-only" if ctrp_ok else "no")
            crs = f"{cr:>8.3f}" if not np.isnan(cr) else f"{'--':>8}"
            cps = f"{cp:>9.2g}" if not np.isnan(cp) else f"{'--':>9}"
            print(f"  {g:<13}{gr:>8.3f}{pr:>9.3f}{pp:>9.2g}{crs}{cps}  {verdict}")
            all_rows.append({"drug": label, "gene": g, "gdsc_r": gr, "prism_r": pr,
                             "prism_p": pp, "ctrp_r": cr, "ctrp_p": cp,
                             "prism_replicates": prism_ok, "ctrp_replicates": ctrp_ok,
                             "resp_r": v.resp_r, "silencing_r": v.silencing_r})

    out = pd.DataFrame(all_rows)
    (ROOT / "outputs").mkdir(exist_ok=True)
    out.to_csv(ROOT / "outputs" / "external_replication.csv", index=False)
    bol = out[out["prism_replicates"]]
    print(f"\n\n######## BULLETPROOF (PRISM-replicated) LEADS: {len(bol)} ########")
    for r in bol.sort_values("resp_r").itertuples():
        d = "meth->SENSITIVE" if r.resp_r < 0 else "meth->resistant"
        extra = " +CTRP" if r.ctrp_replicates else ""
        print(f"  {r.drug:<13} {r.gene:<12} gdsc_r={r.gdsc_r:+.2f} prism_r={r.prism_r:+.2f} ({d}){extra}")
    print(f"\n-> outputs/external_replication.csv")


if __name__ == "__main__":
    main()
