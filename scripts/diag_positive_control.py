"""Honest positive-control diagnostic: does unsupervised whole-methylome clustering
recover the known CDKN2A-methylation -> Palbociclib-sensitivity signal, or does it
only fall out of a targeted biomarker split?

Run: uv run python scripts/diag_positive_control.py
"""
import pandas as pd
from scipy.stats import mannwhitneyu

from strata.data import gdsc
from strata.engine.features import select_variable_features
from strata.engine.cluster import cluster_samples
from strata.engine.associate import test_drug_association
from strata.engine.annotate import differential_cpgs
from strata.config import N_TOP_FEATURES

DRUG = "Palbociclib"
GENE = "CDKN2A"


def _load():
    meth = gdsc.load_gdsc_methylation()
    drug = gdsc.load_gdsc_drug_response()
    auc = (
        drug[drug["drug_name"].str.lower() == DRUG.lower()]
        .set_index("COSMIC_ID")["auc"]
    )
    auc = auc[~auc.index.duplicated(keep="first")]
    common = meth.index.intersection(auc.index)
    return meth.loc[common], auc.loc[common]


def path_a_unsupervised(meth, auc):
    print("=" * 70)
    print("(a) UNSUPERVISED PATH — whole-methylome clustering")
    print("=" * 70)
    feats = select_variable_features(meth, N_TOP_FEATURES)
    clusters = cluster_samples(meth, feats, k=None)
    assoc = test_drug_association(clusters, auc)
    sig = assoc.pvalue < 0.05
    print(f"  responder = cluster {assoc.responder_label} (grade {assoc.grade})")
    print(f"  association p = {assoc.pvalue:.3g}  -> significant? {'YES' if sig else 'NO'}")
    print(f"  AUC effect (baseline - responder) = {assoc.effect_size:.3f}")

    diff = differential_cpgs(meth, clusters, assoc.responder_label).reset_index(drop=True)
    top20 = diff.head(20)["cpg"].tolist()
    in_top20 = GENE in top20
    # full rank (1-based) of CDKN2A among all differential features
    rank = None
    if GENE in diff["cpg"].values:
        rank = int(diff.index[diff["cpg"] == GENE][0]) + 1
    print(f"  {GENE} in top-20 drivers? {'YES' if in_top20 else 'NO'}")
    print(f"  {GENE} differential rank = {rank} of {len(diff)}")
    print(f"  top-20 drivers: {top20}")
    return {
        "sig": sig,
        "p": assoc.pvalue,
        "cdkn2a_top20": in_top20,
        "cdkn2a_rank": rank,
        "n_features": len(diff),
    }


def path_b_targeted(meth, auc):
    print("=" * 70)
    print("(b) TARGETED / BIOMARKER PATH — median split on CDKN2A methylation")
    print("=" * 70)
    beta = meth[GENE]
    common = beta.index.intersection(auc.index)
    beta, a = beta.loc[common], auc.loc[common]
    med = beta.median()
    high = a[beta >= med]
    low = a[beta < med]
    # one-sided: high-CDKN2A-methylation lines have LOWER Palbociclib AUC (more sensitive)
    stat, p = mannwhitneyu(high, low, alternative="less")
    print(f"  n total = {len(common)} (high={len(high)}, low={len(low)})  CDKN2A median beta = {med:.3f}")
    print(f"  median AUC  high-CDKN2A-methyl = {high.median():.3f}")
    print(f"  median AUC  low-CDKN2A-methyl  = {low.median():.3f}")
    print(f"  Mann-Whitney one-sided (high < low) p = {p:.3g}  -> significant? "
          f"{'YES' if p < 0.05 else 'NO'}")
    return {
        "p": float(p),
        "n": len(common),
        "high_med": float(high.median()),
        "low_med": float(low.median()),
        "sig": p < 0.05,
    }


def main():
    meth, auc = _load()
    a = path_a_unsupervised(meth, auc)
    b = path_b_targeted(meth, auc)

    print("=" * 70)
    print("VERDICT")
    print("=" * 70)
    a_clean = a["sig"] and a["cdkn2a_top20"]
    b_clean = b["sig"]
    paths = []
    if a_clean:
        paths.append("unsupervised")
    if b_clean:
        paths.append("targeted")
    if a["sig"] and not a["cdkn2a_top20"]:
        unsup_note = (f"unsupervised clusters ARE significant for {DRUG} (p={a['p']:.2g}) "
                      f"but {GENE} is NOT a top-20 driver (rank {a['cdkn2a_rank']}/"
                      f"{a['n_features']}) — clustering tracks a different axis (likely lineage)")
    elif not a["sig"]:
        unsup_note = "unsupervised path not significant"
    else:
        unsup_note = "unsupervised path recovers CDKN2A cleanly"

    print(f"  Path(s) that CLEANLY demonstrate {GENE}->{DRUG}: "
          f"{', '.join(paths) if paths else 'NONE'}")
    print(f"  Note: {unsup_note}.")
    print(f"  Targeted split: high-CDKN2A-methyl AUC {b['high_med']:.3f} vs "
          f"low {b['low_med']:.3f}, p={b['p']:.2g} ({'sig' if b['sig'] else 'ns'}).")


if __name__ == "__main__":
    main()
