"""v2 mechanism-anchored discovery (R3a) — the corrected, demoable engine.

Per R1b, brute-force tissue-adjusted ranking is non-specific (methylation
co-correlation). The fix: for each drug, restrict the search to its OWN
mechanism genes — the GDSC `PUTATIVE_TARGET` gene(s) + a curated gene set for
its `PATHWAY_NAME` — and find the one whose methylation carries a genuinely
TISSUE-AGNOSTIC association with the drug's response (survives within-tissue +
tissue adjustment, §8a). Each hit is interpretable by construction and
evidence-graded (L0-L3).

Run: PYTHONPATH=. uv run python scripts/v2_mechanism_anchored.py
"""

from __future__ import annotations

import glob
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import spearmanr

from strata.data import gdsc

MIN_TISSUE_N = 20

# Curated PATHWAY_NAME -> genes whose methylation is mechanistically meaningful
# (tumor suppressors / regulators emphasized; license-clean, auditable).
PATHWAY_GENES = {
    "Cell cycle": ["CDKN2A", "CDKN2B", "MTAP", "RB1", "CCND1", "CDK4", "CDK6",
                   "E2F1", "CDKN1A", "CDKN1B", "CCNE1", "MYC", "CHEK1", "CHEK2"],
    "PI3K/MTOR signaling": ["PTEN", "PIK3CA", "AKT1", "MTOR", "TSC1", "TSC2",
                            "PIK3R1", "STK11", "INPP4B", "RICTOR", "RPTOR"],
    "ERK MAPK signaling": ["KRAS", "NRAS", "BRAF", "MAP2K1", "MAP2K2", "NF1",
                           "DUSP6", "SPRY2", "RAF1", "MAPK1", "MAPK3"],
    "EGFR signaling": ["EGFR", "ERBB2", "ERBB3", "ERBB4", "GRB2"],
    "RTK signaling": ["KIT", "PDGFRA", "PDGFRB", "KDR", "MET", "FGFR1", "FGFR2",
                      "NTRK1", "RET", "ALK", "FLT3"],
    "IGF1R signaling": ["IGF1R", "IRS1", "IRS2", "IGFBP3"],
    "Genome integrity": ["BRCA1", "BRCA2", "ATM", "ATR", "PARP1", "PARP2",
                         "RAD51", "PALB2", "CHEK1", "CHEK2", "MLH1", "MSH2", "TP53BP1"],
    "DNA replication": ["TOP1", "TOP2A", "TYMS", "RRM1", "MGMT", "BRCA1", "ERCC1"],
    "Apoptosis regulation": ["BCL2", "BCL2L1", "MCL1", "BAX", "TP53", "BBC3",
                             "PMAIP1", "BCL2L11", "BIRC5"],
    "p53 pathway": ["TP53", "MDM2", "MDM4", "CDKN2A", "ATM", "CHEK2"],
    "Chromatin histone acetylation": ["HDAC1", "HDAC2", "EP300", "CREBBP",
                                      "KAT2B", "BRD4"],
    "Chromatin histone methylation": ["EZH2", "KMT2A", "KMT2D", "KDM6A",
                                      "SETD2", "DOT1L", "NSD1", "WHSC1"],
    "Chromatin other": ["ARID1A", "SMARCA4", "SMARCB1", "PBRM1", "BAP1"],
    "WNT signaling": ["APC", "CTNNB1", "GSK3B", "AXIN1", "AXIN2", "WIF1",
                      "SFRP1", "SFRP2", "DKK1"],
    "Mitosis": ["AURKA", "AURKB", "PLK1", "TTK", "BUB1", "BUB1B"],
    "Metabolism": ["IDH1", "IDH2", "FH", "SDHB", "LDHA"],
    "Protein stability and degradation": ["HSP90AA1", "VHL", "FBXW7", "CRBN"],
    "Hormone-related": ["ESR1", "AR", "PGR"],
}

# Light target-string -> HGNC symbol aliasing for the clean cases.
TARGET_ALIAS = {
    "MEK1": "MAP2K1", "MEK2": "MAP2K2", "ABL": "ABL1", "BCL-XL": "BCL2L1",
    "BCL-W": "BCL2L2", "VEGFR": "KDR", "PDGFR": "PDGFRA", "HER2": "ERBB2",
    "AKT": "AKT1", "PI3K": "PIK3CA", "JAK2": "JAK2", "SRC": "SRC",
}


def parse_targets(s, valid):
    if not isinstance(s, str):
        return []
    out = []
    for tok in s.replace(";", ",").split(","):
        t = tok.strip()
        t = TARGET_ALIAS.get(t, t)
        if t and t[0].isalpha() and t.upper() == t.replace("-", "").replace(" ", "").upper() and len(t) <= 8:
            if t in valid:
                out.append(t)
    return out


def _auc(drug, name):
    s = (drug[drug["drug_name"].str.lower() == name.lower()]
         .drop_duplicates("COSMIC_ID").set_index("COSMIC_ID")["auc"])
    return s[~s.index.duplicated(keep="first")].astype(float)


def eval_gene(meth, auc, tissue, gene):
    """Tissue-agnostic evaluation of one gene's methylation vs a drug's AUC."""
    common = meth.index.intersection(auc.index).intersection(tissue.index)
    m = meth.loc[common, gene].astype(float)
    a = auc.loc[common].astype(float)
    t = tissue.loc[common]
    keep = m.notna() & a.notna()
    m, a, t = m[keep], a[keep], t[keep]
    # restrict to tissues with enough lines
    big = t.map(t.value_counts()) >= MIN_TISSUE_N
    m, a, t = m[big], a[big], t[big]
    if len(m) < 60 or t.nunique() < 3:
        return None
    n, n_tis = len(m), t.nunique()

    rho_p, _ = spearmanr(m, a)
    # within-tissue centering (fixed effects on tissue)
    a_c = (a - a.groupby(t).transform("mean"))
    m_c = (m - m.groupby(t).transform("mean"))
    a_c, m_c = a_c - a_c.mean(), m_c - m_c.mean()
    denom = np.linalg.norm(m_c) * np.linalg.norm(a_c)
    r_adj = float((m_c @ a_c) / denom) if denom else np.nan
    df = n - n_tis - 1
    p_adj = float(2 * stats.t.sf(abs(r_adj) * np.sqrt(df / (1 - r_adj**2)), df)) if abs(r_adj) < 1 else 0.0
    retained = abs(r_adj) / abs(rho_p) if rho_p else np.nan

    # within-tissue consistency
    neg = sig = 0
    for _, g in pd.DataFrame({"m": m, "a": a, "t": t}).groupby("t"):
        if len(g) >= MIN_TISSUE_N:
            r, p = spearmanr(g["m"], g["a"])
            if r < 0:
                neg += 1
                if p < 0.05:
                    sig += 1

    # evidence level — graded on the TISSUE-ADJUSTED statistic + within-tissue
    # consistency (NOT the pooled/adj ratio, which is ill-defined when pooled~0).
    # negative direction = methylation marks sensitivity.
    direction_ok = r_adj < 0
    if direction_ok and p_adj < 0.05 and abs(r_adj) >= 0.12 and sig >= 3:
        level = "L3"   # tissue-agnostic: survives adjustment + consistent in >=3 tissues
    elif direction_ok and p_adj < 0.05 and (abs(r_adj) >= 0.08 or sig >= 2):
        level = "L2"   # survives adjustment, weaker cross-tissue support
    elif p_adj < 0.05 or (rho_p < 0 and abs(rho_p) > 0.1):
        level = "L1"   # pooled signal that doesn't hold up = likely lineage
    else:
        level = "L0"
    return dict(gene=gene, rho_pooled=float(rho_p), r_adj=r_adj, p_adj=p_adj,
                retained=float(retained), within_neg=neg, within_sig=sig,
                n=n, n_tissues=n_tis, level=level)


def main():
    meth = gdsc.load_gdsc_methylation()
    drug = gdsc.load_gdsc_drug_response()
    ann = gdsc.load_gdsc_annotations()
    tissue = ann["tissue"]
    valid = set(meth.columns)

    f = glob.glob("data/gdsc/GDSC2_fitted_dose_response_*.xlsx")[0]
    meta = (pd.read_excel(f, usecols=["DRUG_NAME", "PUTATIVE_TARGET", "PATHWAY_NAME"])
            .drop_duplicates("DRUG_NAME"))

    order = {"L3": 3, "L2": 2, "L1": 1, "L0": 0}
    rows = []
    for _, mr in meta.iterrows():
        dn, tgt, pw = mr["DRUG_NAME"], mr["PUTATIVE_TARGET"], mr["PATHWAY_NAME"]
        genes = set(parse_targets(tgt, valid)) | set(g for g in PATHWAY_GENES.get(pw, []) if g in valid)
        if not genes:
            continue
        auc = _auc(drug, dn)
        if len(auc) < 100:
            continue
        evals = [e for e in (eval_gene(meth, auc, tissue, g) for g in genes) if e]
        evals = [e for e in evals if e["r_adj"] < 0]  # sensitivity-direction marks
        if not evals:
            continue
        best = max(evals, key=lambda e: (order[e["level"]], abs(e["r_adj"])))
        rows.append((dn, str(pw), best["gene"], best["level"], best["r_adj"], best["rho_pooled"],
                     best["p_adj"], best["retained"], best["within_sig"], best["n_tissues"]))

    res = pd.DataFrame(rows, columns=["drug", "pathway", "marker", "level", "r_adj", "rho_pooled",
                                      "p_adj", "retained", "within_sig", "n_tissues"])
    res = res.sort_values(["level", "r_adj"], key=lambda c: c.map(order) if c.name == "level" else c.abs(),
                          ascending=[False, False]).reset_index(drop=True)

    print("=" * 92)
    print("MECHANISM-ANCHORED TISSUE-AGNOSTIC OPPORTUNITIES (best mechanism gene per drug)")
    print("=" * 92)
    print(f"{'drug':<24}{'pathway':<27}{'marker':<9}{'lvl':<4}{'r_adj':>7}{'r_pool':>8}{'wsig':>5}")
    for _, r in res[res["level"].isin(["L3", "L2"])].head(30).iterrows():
        print(f"{r['drug']:<24}{r['pathway']:<27}{r['marker']:<9}{r['level']:<4}"
              f"{r['r_adj']:>7.3f}{r['rho_pooled']:>8.3f}{r['within_sig']:>5}")

    print("\n--- positive controls (CDK4/6 inhibitors -> 9p21) ---")
    for d in ("Palbociclib", "Ribociclib"):
        pc = res[res["drug"].str.lower() == d.lower()]
        print(pc.to_string(index=False) if len(pc) else f"  ({d} not in results)")

    # explicit MGMT -> temozolomide negative control
    auc_tmz = _auc(drug, "Temozolomide")
    mgmt = eval_gene(meth, auc_tmz, tissue, "MGMT") if "MGMT" in valid else None
    print("\n--- negative control: MGMT -> Temozolomide (expect null) ---")
    print(f"  MGMT r_adj={mgmt['r_adj']:+.3f} p_adj={mgmt['p_adj']:.2g} "
          f"within_sig={mgmt['within_sig']} level={mgmt['level']}" if mgmt else "  (n/a)")

    print(f"\nlevel counts: {res['level'].value_counts().to_dict()}")
    res.to_csv("outputs/v2_mechanism_anchored.csv", index=False)
    print("wrote outputs/v2_mechanism_anchored.csv")


if __name__ == "__main__":
    main()
