"""Build the MASTER table: one row per (cell line x drug), with metadata,
genomics (TMB/MSI/ploidy/...), drug sensitivity, and a curated cancer-gene
methylation panel (gene-level AND promoter beta per gene).

Also exports the FULL methylation matrices as separate COSMIC-keyed joinable
CSVs (gzipped). Outputs -> data/master/.

Run: PYTHONPATH=. uv run python scripts/build_master_csv.py
"""

from __future__ import annotations

import glob
from pathlib import Path

import pandas as pd

from strata.config import ROOT
from strata.data import gdsc
from strata.data import ccle_methylation as ccle

OUT = ROOT / "data" / "master"

# Curated cancer-gene panel (tumor suppressors, oncogenes, DDR, epigenetic,
# + the demo/mechanism-marker genes). Methylation reported gene-level + promoter.
PANEL = [
    # demo / mechanism markers
    "MGMT", "CDKN2A", "CDKN2B", "MTAP", "INPP4B", "SFRP1", "DAPK1", "RASSF1",
    # tumor suppressors
    "TP53", "RB1", "PTEN", "APC", "VHL", "NF1", "NF2", "SMAD4", "STK11",
    "BRCA1", "BRCA2", "ATM", "MLH1", "MSH2", "PMS2", "MSH6", "WT1", "BAP1",
    "SMARCA4", "ARID1A", "PBRM1", "FBXW7", "KEAP1", "SETD2",
    # oncogenes
    "KRAS", "NRAS", "BRAF", "EGFR", "ERBB2", "MET", "ALK", "MYC", "MYCN",
    "CCND1", "CDK4", "CDK6", "MDM2", "MDM4", "PIK3CA", "AKT1", "MTOR",
    "BCL2", "MCL1", "FLT3", "KIT", "JAK2", "IDH1", "IDH2", "EZH2", "KMT2D",
    # DDR / epigenetic
    "PARP1", "HDAC1", "ATR", "CHEK1", "CHEK2",
]

# NOTE: cell_line + sanger_model_id already come from the drug file (fact table),
# so they are NOT remapped here (avoids _x/_y merge collisions).
META_COLS = {
    "tissue": "tissue",
    "mutational_burden": "TMB", "msi_status": "MSI",
    "mismatch_repair_status": "mmr_status",
    "mlh1_promoter_methylation_status": "mlh1_promoter_meth",
    "gender": "gender", "age_at_sampling": "age",
    "braf_mutation_identified": "braf_mut", "kras_mutation_identified": "kras_mut",
    "pik3ca_mutation_identified": "pik3ca_mut", "pten_mutation_identified": "pten_mut",
}


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    # --- fact table: one row per (cell line x drug) ---
    f = glob.glob(str(ROOT / "data/gdsc/GDSC2_fitted_dose_response_*.xlsx"))[0]
    drug = pd.read_excel(f, usecols=["CELL_LINE_NAME", "SANGER_MODEL_ID", "CANCER_TYPE",
                                     "DRUG_NAME", "PUTATIVE_TARGET", "PATHWAY_NAME",
                                     "LN_IC50", "AUC"])
    ml = pd.read_csv(ROOT / "data/gdsc/model_list_20260420.csv", dtype=str, low_memory=False)
    sidm2cosmic = ml.set_index("model_id")["COSMIC_ID"].dropna()
    drug["COSMIC_ID"] = drug["SANGER_MODEL_ID"].map(sidm2cosmic)
    drug = drug.dropna(subset=["COSMIC_ID"])
    drug["COSMIC_ID"] = drug["COSMIC_ID"].astype(str)

    fact = pd.DataFrame({
        "COSMIC_ID": drug["COSMIC_ID"].values,
        "cell_line": drug["CELL_LINE_NAME"].values,
        "sanger_model_id": drug["SANGER_MODEL_ID"].values,
        "cancer_type": drug["CANCER_TYPE"].values,
        "drug": drug["DRUG_NAME"].values,
        "drug_target": drug["PUTATIVE_TARGET"].values,
        "drug_pathway": drug["PATHWAY_NAME"].values,
        "ln_ic50": drug["LN_IC50"].values,
        "auc": drug["AUC"].values,
    })

    # --- metadata + genomics (per COSMIC) ---
    meta = ml[ml["COSMIC_ID"].notna()].copy()
    meta["COSMIC_ID"] = meta["COSMIC_ID"].astype(str)
    # coalesce ploidy across the three assays (WGS -> WES -> SNP6) for coverage
    meta["ploidy"] = (meta.get("ploidy_wgs").fillna(meta.get("ploidy_wes"))
                      .fillna(meta.get("ploidy_snp6")))
    keep = ["COSMIC_ID", "ploidy"] + [c for c in META_COLS if c in meta.columns]
    meta = meta[keep].rename(columns=META_COLS).drop_duplicates("COSMIC_ID")
    # 'tissue' from meta (broad grouping); cancer_type already in fact
    fact = fact.merge(meta, on="COSMIC_ID", how="left")

    # --- methylation panel: gene-level + promoter ---
    mg = gdsc.load_gdsc_methylation()
    mp = ccle.load_ccle_promoter_methylation()
    panel = [g for g in PANEL if g in mg.columns]  # gene-level is the superset
    mg_p = mg[panel].copy(); mg_p.columns = [f"meth_gene__{g}" for g in panel]
    mg_p.index = mg_p.index.astype(str)
    prom_genes = [g for g in panel if g in mp.columns]
    mp_p = mp[prom_genes].copy(); mp_p.columns = [f"meth_prom__{g}" for g in prom_genes]
    mp_p.index = mp_p.index.astype(str)

    fact = fact.merge(mg_p, left_on="COSMIC_ID", right_index=True, how="left")
    fact = fact.merge(mp_p, left_on="COSMIC_ID", right_index=True, how="left")

    master = OUT / "master_cellline_drug.csv"
    fact.to_csv(master, index=False)

    # --- full methylation matrices as separate joinable CSVs (gzipped) ---
    g_out = OUT / "methylation_genelevel_bycosmic.csv.gz"
    p_out = OUT / "methylation_promoter_bycosmic.csv.gz"
    mg.rename_axis("COSMIC_ID").to_csv(g_out, compression="gzip")
    mp.rename_axis("COSMIC_ID").to_csv(p_out, compression="gzip")

    print(f"MASTER: {fact.shape[0]:,} rows x {fact.shape[1]} cols -> {master}")
    print(f"  size: {master.stat().st_size/1e6:.0f} MB")
    print(f"  panel genes: {len(panel)} (gene-level), {len(prom_genes)} (promoter)")
    print(f"  columns: {list(fact.columns[:15])} ... + {len(panel)+len(prom_genes)} methylation cols")
    print(f"FULL gene-level matrix: {mg.shape} -> {g_out} ({g_out.stat().st_size/1e6:.0f} MB)")
    print(f"FULL promoter matrix:   {mp.shape} -> {p_out} ({p_out.stat().st_size/1e6:.0f} MB)")
    print("\n=== master head (key cols) ===")
    show = ["cell_line", "COSMIC_ID", "cancer_type", "drug", "auc", "TMB", "MSI",
            "ploidy", "meth_gene__CDKN2A", "meth_prom__CDKN2A"]
    print(fact[[c for c in show if c in fact.columns]].head(6).to_string(index=False))


if __name__ == "__main__":
    main()
