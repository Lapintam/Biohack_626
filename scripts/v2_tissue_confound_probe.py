"""v2 crux probe — does a marker->drug signal SURVIVE the tissue-confound safeguards?

The whole cross-indication / molecular-state thesis rests on one empirical claim:
a state's drug-response association is **tissue-agnostic** — it holds *within* tissues
and *after adjusting for* tissue, not just pooled (which could be lineage, per the
unsupervised null). This script applies the SYSTEM_MAP §8a safeguards to candidate
marker->drug pairs and prints a verdict: tissue-agnostic (L2+) vs lineage-confounded.

Run: PYTHONPATH=. uv run python scripts/v2_tissue_confound_probe.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import statsmodels.formula.api as smf

from strata.data import gdsc

PAIRS = [
    ("MTAP", "Palbociclib"),    # the 9p21 positive control (our top marker)
    ("CDKN2A", "Palbociclib"),  # the canonical 9p21 partner
    ("DAPK1", "Luminespib"),    # a novel hit — is it tissue-agnostic or lineage?
]
MIN_TISSUE_N = 20   # min lines in a tissue to trust a within-tissue estimate
RARE_N = 15         # tissues smaller than this are pooled into "Other" for the model


def _frame(meth, drug, ann, marker, drug_name):
    auc = (
        drug[drug["drug_name"].str.lower() == drug_name.lower()]
        .drop_duplicates("COSMIC_ID").set_index("COSMIC_ID")["auc"]
    )
    common = meth.index.intersection(auc.index).intersection(ann.index)
    df = pd.DataFrame({
        "meth": meth.loc[common, marker].astype(float),
        "auc": auc.loc[common].astype(float),
        "tissue": ann.loc[common, "tissue"].astype(str),
    }).dropna()
    return df


def probe(df, marker, drug_name):
    print("=" * 74)
    print(f"{marker} -> {drug_name}   (n={len(df)} lines, {df['tissue'].nunique()} tissues)")
    print("=" * 74)

    # (1) pooled
    rho_p, p_p = spearmanr(df["meth"], df["auc"])
    print(f"  [pooled]        rho={rho_p:+.3f}  p={p_p:.2g}")

    # (2) within-tissue
    within = []
    for t, g in df.groupby("tissue"):
        if len(g) >= MIN_TISSUE_N:
            r, p = spearmanr(g["meth"], g["auc"])
            within.append((t, len(g), r, p))
    within.sort(key=lambda x: -x[1])
    n_neg = sum(1 for _, _, r, _ in within if r < 0)
    n_sig_neg = sum(1 for _, _, r, p in within if r < 0 and p < 0.05)
    med_within = np.median([r for _, _, r, _ in within]) if within else float("nan")
    print(f"  [within-tissue] {len(within)} tissues n>={MIN_TISSUE_N}: "
          f"{n_neg} negative, {n_sig_neg} sig-negative, median rho={med_within:+.3f}")
    for t, n, r, p in within:
        flag = "  <-- sig" if (r < 0 and p < 0.05) else ""
        print(f"        {t:<32} n={n:<4} rho={r:+.3f} p={p:.2g}{flag}")

    # (3) between-tissue component (could the pooled signal be only lineage?)
    tm = df.groupby("tissue").agg(meth=("meth", "mean"), auc=("auc", "mean"),
                                  n=("auc", "size"))
    tm = tm[tm["n"] >= MIN_TISSUE_N]
    if len(tm) >= 3:
        rb, pb = spearmanr(tm["meth"], tm["auc"])
        print(f"  [between-tissue] tissue-mean rho={rb:+.3f} p={pb:.2g}  "
              f"(high => pooled signal partly lineage)")

    # (4) tissue-adjusted model: AUC ~ z(meth) + C(tissue)
    d = df.copy()
    d["methz"] = (d["meth"] - d["meth"].mean()) / d["meth"].std()
    counts = d["tissue"].value_counts()
    d["tissue_m"] = np.where(d["tissue"].map(counts) >= RARE_N, d["tissue"], "Other")
    m_only = smf.ols("auc ~ methz", data=d).fit()
    m_adj = smf.ols("auc ~ methz + C(tissue_m)", data=d).fit()
    b0, p0 = m_only.params["methz"], m_only.pvalues["methz"]
    b1, p1 = m_adj.params["methz"], m_adj.pvalues["methz"]
    retain = abs(b1) / abs(b0) if b0 != 0 else float("nan")
    print(f"  [tissue-adj]    marker-only coef={b0:+.4f} (p={p0:.2g})  ->  "
          f"tissue-adjusted coef={b1:+.4f} (p={p1:.2g})")
    print(f"                  effect retained after adjustment: {retain*100:.0f}%")

    # verdict
    survives_model = (np.sign(b1) == np.sign(b0)) and (p1 < 0.05) and (retain > 0.4)
    survives_within = len(within) >= 3 and n_neg >= max(3, 0.6 * len(within))
    if survives_model and survives_within:
        verdict = "TISSUE-AGNOSTIC (survives within-tissue + adjustment) -> L2/L3"
    elif survives_model or survives_within:
        verdict = "PARTIAL — some tissue-agnostic signal, not clean (investigate)"
    else:
        verdict = "LINEAGE-CONFOUNDED — pooled signal does not survive (NOT L2)"
    print(f"  VERDICT: {verdict}\n")
    return verdict


def main():
    meth = gdsc.load_gdsc_methylation()
    drug = gdsc.load_gdsc_drug_response()
    ann = gdsc.load_gdsc_annotations()
    for marker, drug_name in PAIRS:
        if marker not in meth.columns:
            print(f"skip {marker}: not in methylation matrix\n")
            continue
        df = _frame(meth, drug, ann, marker, drug_name)
        probe(df, marker, drug_name)


if __name__ == "__main__":
    main()
