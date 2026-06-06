"""Scan unsupervised pipeline across drugs for significant responder subgroups.
Run: PYTHONPATH=. uv run python scripts/scan_novel_hits.py
"""
import pandas as pd
from strata.data import gdsc
from strata.engine.features import select_variable_features
from strata.engine.cluster import cluster_samples
from strata.engine.associate import test_drug_association
from strata.engine.annotate import differential_cpgs, annotate_cpgs
from strata.config import N_TOP_FEATURES

DRUGS = ["Trametinib", "Olaparib", "Venetoclax", "Dasatinib", "Nutlin-3a (-)",
         "Alpelisib", "Lapatinib", "Erlotinib", "Sorafenib"]

meth = gdsc.load_gdsc_methylation()
drug = gdsc.load_gdsc_drug_response()

# Cluster ONCE over all methylation samples (drug-independent), then subset per drug.
feats = select_variable_features(meth, N_TOP_FEATURES)
clusters_all = cluster_samples(meth, feats, k=None)
print(f"# global clusters: {dict(clusters_all.value_counts().sort_index())}")

rows = []
for dn in DRUGS:
    auc = drug[drug["drug_name"] == dn].set_index("COSMIC_ID")["auc"]
    auc = auc[~auc.index.duplicated(keep="first")]
    common = meth.index.intersection(auc.index)
    if len(common) < 30:
        print(f"{dn}: SKIP n={len(common)}")
        continue
    cl = clusters_all.loc[common]
    auc_c = auc.loc[common]
    assoc = test_drug_association(cl, auc_c)
    diff = differential_cpgs(meth.loc[common], cl, assoc.responder_label).head(20)
    genes = annotate_cpgs(diff["cpg"].tolist())["gene"].unique().tolist()
    rsize = int((cl == assoc.responder_label).sum())
    rows.append((dn, len(common), rsize, assoc.responder_label, assoc.grade,
                 assoc.pvalue, assoc.effect_size, genes[:8]))
    print(f"{dn}: n={len(common)} responder=c{assoc.responder_label}(n={rsize}) "
          f"grade={assoc.grade} p={assoc.pvalue:.3g} effect={assoc.effect_size:.3f} "
          f"drivers={genes[:6]}")

print("\n# === ranked by effect among p<0.05 & grade in {+,++} ===")
good = [r for r in rows if r[5] < 0.05 and r[4] in ("+", "++")]
for r in sorted(good, key=lambda x: -x[6]):
    print(f"{r[0]}: grade={r[4]} p={r[5]:.2g} effect={r[6]:.3f} "
          f"n_resp={r[2]} drivers={r[7][:6]}")
