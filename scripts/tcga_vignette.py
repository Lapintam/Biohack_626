"""TCGA-LAML patient vignette.

Cluster acute myeloid leukemia patients on DNA methylation and test whether the
methylation subgroups differ in overall survival (log-rank). Real-patient anchor.

Run: ``uv run python scripts/tcga_vignette.py``
"""

from __future__ import annotations

import sys
from pathlib import Path

# Repo root on path so `strata` imports without an installed package.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from strata.config import N_TOP_FEATURES, RANDOM_STATE
from strata.data import tcga_laml
from strata.engine.associate import test_survival_association
from strata.engine.cluster import cluster_samples
from strata.engine.features import select_variable_features
from strata.figures import plot_km, plot_umap


def main() -> None:
    meth = tcga_laml.load_tcga_methylation()
    clin = tcga_laml.load_tcga_clinical()

    common = meth.index.intersection(clin.index)
    meth = meth.loc[common]
    clin = clin.loc[common]
    print(f"n patients (methylation x survival): {len(common)}")
    print(f"methylation matrix: {meth.shape[0]} samples x {meth.shape[1]} probes")

    features = select_variable_features(meth, n_top=N_TOP_FEATURES)
    clusters = cluster_samples(meth, features, k=None, random_state=RANDOM_STATE)
    n_clusters = clusters.nunique()
    print(f"n clusters (auto-k via silhouette): {n_clusters}")
    print("per-cluster size:")
    for label, size in clusters.value_counts().sort_index().items():
        print(f"  subgroup {int(label)}: {int(size)} patients")

    surv = test_survival_association(clusters, clin["os_time"], clin["os_event"])
    p = surv["logrank_p"]
    print(f"log-rank p = {p:.4g}")
    verdict = "SIGNIFICANT (p < 0.05)" if p < 0.05 else "NOT significant (p >= 0.05)"
    print(f"survival difference across methylation subgroups: {verdict}")

    km_path = plot_km(surv, fname="tcga_laml_km.png")
    umap_path = plot_umap(meth, features, clusters, fname="tcga_laml_umap.png")
    print(f"wrote {km_path}")
    print(f"wrote {umap_path}")


if __name__ == "__main__":
    main()
