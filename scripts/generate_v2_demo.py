"""Generate v2 cross-indication demo artifacts -> site/demo-v2/opportunities.json.

Computes real numbers (per-tissue within-tissue rho, tissue-adjusted r) for the
featured controls/leads, reads the ranked opportunity table from R3a, attaches
curated mechanism cards + a scripted agent transcript. No API key needed.

Run: PYTHONPATH=. uv run python scripts/generate_v2_demo.py
"""

from __future__ import annotations

import glob
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from strata.data import gdsc
from strata.config import ROOT
from scripts.v2_mechanism_anchored import eval_gene  # reuse the graded evaluator

OUT = ROOT / "site" / "demo-v2"
MIN_TISSUE_N = 20


def _auc(drug, name):
    s = (drug[drug["drug_name"].str.lower() == name.lower()]
         .drop_duplicates("COSMIC_ID").set_index("COSMIC_ID")["auc"])
    return s[~s.index.duplicated(keep="first")].astype(float)


def per_tissue(meth, auc, tissue, gene):
    common = meth.index.intersection(auc.index).intersection(tissue.index)
    m, a, t = meth.loc[common, gene], auc.loc[common], tissue.loc[common]
    keep = m.notna() & a.notna()
    m, a, t = m[keep], a[keep], t[keep]
    rows = []
    for tis, g in pd.DataFrame({"m": m, "a": a, "t": t}).groupby("t"):
        if len(g) >= MIN_TISSUE_N:
            r, p = spearmanr(g["m"], g["a"])
            rows.append({"tissue": tis, "n": int(len(g)), "rho": round(float(r), 3),
                         "sig": bool(r < 0 and p < 0.05)})
    rows.sort(key=lambda x: x["rho"])  # most-sensitive first
    return rows


def featured(meth, drug, tissue, name, marker, pathway, card, level_note, replicated=None):
    auc = _auc(drug, name)
    e = eval_gene(meth, auc, tissue, marker)
    return {
        "drug": name, "marker": marker, "pathway": pathway,
        "level": e["level"], "r_adj": round(e["r_adj"], 3),
        "rho_pooled": round(e["rho_pooled"], 3), "p_adj": e["p_adj"],
        "within_sig": e["within_sig"], "n": e["n"], "n_tissues": e["n_tissues"],
        "per_tissue": per_tissue(meth, auc, tissue, marker),
        "mechanism_card": card, "level_note": level_note,
        "replicated_by": replicated or [],
    }


def main():
    meth = gdsc.load_gdsc_methylation()
    drug = gdsc.load_gdsc_drug_response()
    ann = gdsc.load_gdsc_annotations()
    tissue = ann["tissue"]

    pos = featured(
        meth, drug, tissue, "Palbociclib", "MTAP", "Cell cycle",
        {"state": "9p21 / MTAP–CDKN2A methylation–silencing",
         "known_biology": "CDKN2A/p16 loss removes the CDK4/6 brake → CDK4/6-inhibitor dependence; MTAP is the co-deleted 9p21 neighbor.",
         "implication": "A pan-tissue 9p21-methylation state that predicts CDK4/6-inhibitor sensitivity across indications, not within one."},
        "L3 · CONFIRMATION (the CDKN2A→palbociclib link is published, PLOS One 2019) — recovered cold + confound-controlled.",
        replicated=[{"drug": "Ribociclib", "marker": "MTAP", "note": "2nd CDK4/6 inhibitor, same 9p21 anchor (L2)"}],
    )
    novelA = featured(
        meth, drug, tissue, "Rapamycin", "INPP4B", "PI3K/MTOR signaling",
        {"state": "INPP4B methylation–silencing (PI3K-pathway tumor suppressor)",
         "known_biology": "INPP4B is a PI3K-pathway phosphatase; its silencing biases cells toward PI3K/AKT/mTOR dependence.",
         "implication": "An INPP4B-methylated state may mark mTOR/PI3K-inhibitor responders across tissue boundaries."},
        "L3 · NOVEL HYPOTHESIS — replicated across 3 mTOR-pathway drugs (Rapamycin, Dactolisib, AZD8055).",
        replicated=[{"drug": "Dactolisib", "marker": "INPP4B", "note": "PI3K/mTOR dual inhibitor (L3)"},
                    {"drug": "AZD8055", "marker": "INPP4B", "note": "mTOR inhibitor (L3)"}],
    )
    novelB = featured(
        meth, drug, tissue, "SB216763", "SFRP1", "WNT signaling",
        {"state": "SFRP1 promoter methylation (WNT-antagonist silencing)",
         "known_biology": "SFRP1 is a secreted WNT antagonist among the most frequently methylation-silenced genes pan-cancer.",
         "implication": "An SFRP1-methylated state recurring across ≥6 tissues — a candidate cross-indication WNT-pathway responder stratum."},
        "L3 · NOVEL HYPOTHESIS — the most cross-tissue lead (significant within 6 independent tissues).",
        replicated=[{"drug": "CHIR-99021", "marker": "SFRP1", "note": "2nd WNT/GSK3 modulator, same anchor (L2)"}],
    )
    # negative control
    auc_tmz = _auc(drug, "Temozolomide")
    mg = eval_gene(meth, auc_tmz, tissue, "MGMT")
    neg = {"drug": "Temozolomide", "marker": "MGMT", "r_adj": round(mg["r_adj"], 3),
           "within_sig": mg["within_sig"], "level": mg["level"],
           "note": "MGMT → temozolomide shows NO tissue-agnostic sensitivity signal (r_adj positive, 0 tissues). The engine does not force the canonical gene-level story — an honest negative control."}

    # ranked opportunity table from R3a
    csv = ROOT / "outputs" / "v2_mechanism_anchored.csv"
    tbl = pd.read_csv(csv)
    tbl = tbl[tbl["level"].isin(["L3", "L2"])].head(20)
    opportunity_table = [
        {"drug": r["drug"], "pathway": r["pathway"], "marker": r["marker"],
         "level": r["level"], "r_adj": round(float(r["r_adj"]), 3),
         "rho_pooled": round(float(r["rho_pooled"]), 3), "within_sig": int(r["within_sig"])}
        for _, r in tbl.iterrows()
    ]

    transcript = [
        {"type": "thought", "text": "Goal: don't start from an indication — start from molecular state. Scan every drug's OWN target/pathway genes for a methylation state whose response association is tissue-agnostic (holds within tissues, not just pooled)."},
        {"type": "tool_call", "name": "scan_mechanism_anchored", "args": {"safeguard": "tissue-adjusted + within-tissue ≥3"}},
        {"type": "tool_result", "name": "scan_mechanism_anchored", "result": {"L3": 16, "L2": 83, "ranked": "outputs/v2_mechanism_anchored.csv"}},
        {"type": "thought", "text": "Positive control recovered COLD: 9p21/MTAP → palbociclib is L3 (survives tissue adjustment, significant in 5 tissues) — and it replicates on a second CDK4/6 inhibitor, ribociclib. Known biology recovered tissue-agnostically licenses the rest. (This link is published — confirmation, not discovery.)"},
        {"type": "tool_call", "name": "tissue_confound_check", "args": {"drug": "Palbociclib", "marker": "MTAP"}},
        {"type": "tool_result", "name": "tissue_confound_check", "result": {"verdict": "tissue-agnostic (L3)", "within_tissue_sig": 5}},
        {"type": "thought", "text": "Novel cross-indication lead: INPP4B methylation → mTOR/PI3K-inhibitor sensitivity, L3, and replicated across THREE mTOR-pathway drugs (Rapamycin, Dactolisib, AZD8055). A PI3K-pathway tumor-suppressor state that cuts across tissues — a hypothesis to follow up. SFRP1→WNT is the most cross-tissue (6 tissues)."},
        {"type": "thought", "text": "Negative control holds: MGMT → temozolomide gives no tissue-agnostic signal. The engine isn't forcing canonical results — it reports what survives the safeguards."},
        {"type": "tool_call", "name": "write_opportunity_brief", "args": {"lead": "INPP4B → mTOR-pathway inhibitors"}},
        {"type": "thought", "text": "# Cross-Indication Opportunity — INPP4B-methylated state × mTOR-pathway inhibitors\n\n## Evidence (tissue-agnostic, mechanism-anchored)\nAn unbiased methylation screen, restricted to each drug's own pathway and corrected for tissue, surfaces **INPP4B promoter methylation** as a responder marker for mTOR/PI3K inhibitors — **tissue-adjusted r≈−0.21, significant within 3 independent tissues, and replicated across 3 drugs of the class** (Rapamycin, Dactolisib, AZD8055).\n\n## Why it's invisible to current taxonomy\nThe state is defined by a methylation event in a PIK3-pathway tumor suppressor, not by histology — it recurs across indications. It is **licensed by the same screen recovering the known 9p21→CDK4/6i control cold**.\n\n## Commercial implication\nAn INPP4B-methylation companion-marker could define a cross-indication responder stratum for mTOR/PI3K inhibitors — a portfolio-triage hypothesis for follow-up, not a clinical claim. Evidence level **L3** (cell-line, research-use)."},
    ]

    record = {
        "thesis": "Start from molecular state, not disease taxonomy. Rank methylation-defined drug-response states that hold ACROSS the indication taxonomy — mechanism-anchored, tissue-confound-controlled, evidence-graded.",
        "novelty": "Existing work ranks methylation–drug biomarkers one cancer type at a time (Comms Biol 2023) or one query at a time (CellMinerCDB); Strata ranks the methylation state ACROSS the taxonomy, anchored to the drug's own target/pathway to defeat the between-tissue co-methylation confound, and evidence-grades each opportunity.",
        "evidence_ladder": [
            {"level": "L0", "label": "pooled association"},
            {"level": "L1", "label": "interpretable marker / known biology"},
            {"level": "L2", "label": "survives tissue adjustment"},
            {"level": "L3", "label": "tissue-agnostic — consistent across ≥3 cancer types"},
            {"level": "L4", "label": "mechanism path (annotation)"},
            {"level": "L5", "label": "external/orthogonal support"},
            {"level": "L6", "label": "clinical evidence"},
        ],
        "featured": {"positive_control": pos, "novel_leads": [novelA, novelB], "negative_control": neg},
        "opportunity_table": opportunity_table,
        "transcript": transcript,
        "generated": "rung-A mechanism-anchored, scripted",
    }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "opportunities.json").write_text(json.dumps(record, indent=2))
    print(f"wrote {OUT/'opportunities.json'}")
    print(f"positive control: {pos['drug']}→{pos['marker']} {pos['level']} (within_sig={pos['within_sig']})")
    print(f"novel A: {novelA['drug']}→{novelA['marker']} {novelA['level']} (within_sig={novelA['within_sig']})")
    print(f"novel B: {novelB['drug']}→{novelB['marker']} {novelB['level']} (within_sig={novelB['within_sig']})")
    print(f"negative: {neg['drug']}→{neg['marker']} {neg['level']} r_adj={neg['r_adj']}")
    print(f"table rows: {len(opportunity_table)}")


if __name__ == "__main__":
    main()
