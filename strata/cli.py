"""Rung-A deterministic discovery pipeline on GDSC (marker-driven)."""
import argparse, pandas as pd
from strata.data import gdsc
from strata.engine.features import select_variable_features
from strata.engine.cluster import cluster_samples
from strata.engine.associate import (test_drug_association, rank_response_markers,
                                     discover_responder_subgroup)
from strata.engine.annotate import annotate_cpgs
from strata import figures
from strata.config import N_TOP_FEATURES, OUTPUT_DIR

def run(drug_name: str, highlight_gene: str | None = None):
    meth = gdsc.load_gdsc_methylation()
    drug = gdsc.load_gdsc_drug_response()
    auc = drug[drug["drug_name"].str.lower() == drug_name.lower()].set_index("COSMIC_ID")["auc"]
    auc = auc[~auc.index.duplicated(keep="first")]
    common = meth.index.intersection(auc.index)
    meth, auc = meth.loc[common], auc.loc[common]

    markers = rank_response_markers(meth, auc)
    top_marker = highlight_gene or markers.iloc[0]["gene"]
    clusters, assoc = discover_responder_subgroup(meth, auc, top_marker)
    top_genes = markers.head(15)["gene"].tolist()

    feats = select_variable_features(meth, N_TOP_FEATURES)
    figures.plot_umap(meth, feats, clusters, response=auc, fname=f"{drug_name}_umap.png")
    figures.plot_auc_box(clusters, auc, fname=f"{drug_name}_auc.png")
    # heatmap over the top discovered markers
    figures.plot_driver_heatmap(meth, clusters, top_genes, fname=f"{drug_name}_heatmap.png")

    from strata.agent.brief import write_brief
    brief = write_brief(drug_name, assoc, [top_marker] + top_genes, len(common))
    (OUTPUT_DIR / f"{drug_name}_brief.md").write_text(brief)
    print(brief)
    print(f"[debug] drug={drug_name} n={len(common)} top_marker={top_marker} "
          f"rho={markers.iloc[0]['rho']:.3f} q={markers.iloc[0]['qvalue']:.2g} "
          f"subgroup_grade={assoc.grade} subgroup_p={assoc.pvalue:.3g} "
          f"top15={top_genes}")
    return markers, assoc

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--drug", default="Palbociclib")
    ap.add_argument("--gene", default=None, help="force highlight a specific marker gene")
    a = ap.parse_args()
    run(a.drug, a.gene)
