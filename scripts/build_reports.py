"""Precompute STANDARDIZED REPORTS (drug + cancer-type) for the report API.

Per drug: global-axis promoter scan -> silencing filter -> external replication
(PRISM/CTRP) -> per-tissue recurrence -> standardized drug report.
Then derive cancer-type reports (which cross-indication leads touch each tissue).
Writes incrementally to api/reports/{drug,cancer_type}/*.json + manifest.json.

The JSON shape here is the contract the FastAPI Pydantic/OpenAPI models mirror.
Designed so a future patient-methylation sample can be scored against the same leads.

Run: PYTHONPATH=. uv run python scripts/build_reports.py
"""

from __future__ import annotations

import json
import re
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

from strata.config import ROOT
from strata.data import gdsc
from strata.data import ccle_methylation as ccle
from strata.data import ccle_expression as cexpr
from strata.data import external_drug as ext
from scripts.v2_global_axis_scan import scan
from scripts.v2_silencing_validation import silencing
from scripts.v2_external_replication import tissue_adj_corr
from scripts.generate_v3_demo import per_tissue_recurrence, GENE_NOTES

OUT = ROOT / "api" / "reports"
TOP_CANDIDATES = 40       # response-candidates fed to silencing
TOP_LEADS = 15            # validated leads kept per drug (replication + recurrence)

CURATED = [
    "Olaparib", "Talazoparib", "Rucaparib", "Niraparib",
    "Trametinib", "Selumetinib", "Refametinib", "PD0325901",
    "Palbociclib", "Temozolomide", "Nutlin-3a (-)", "Rapamycin",
    "AZD8055", "Dabrafenib", "PLX-4720", "Gefitinib", "Erlotinib",
    "Afatinib", "Alpelisib", "MK-2206", "Dasatinib", "Venetoclax",
    "Navitoclax", "Cisplatin", "Gemcitabine", "Vorinostat", "Docetaxel",
]
HONESTY = [
    "Single-gene methylation->response is a weak, diffuse signal; the claim is directional concordance across independent screens + functional silencing, not a large effect.",
    "Cell-line results: defensible novel hypotheses (evidence ladder ~L3-L5), not clinical claims.",
    "The agent expands the hypothesis space the expert adjudicates; it removes the structural gate, it does not replace expert judgment.",
]


def slug(s):
    return re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")


def _auc(drug, name):
    s = (drug[drug["drug_name"].str.lower() == name.lower()]
         .drop_duplicates("COSMIC_ID").set_index("COSMIC_ID")["auc"])
    return s[~s.index.duplicated(keep="first")].astype(float)


def evidence_level(n_screens):
    return {3: "L5", 2: "L4", 1: "L3"}.get(n_screens, "L3")


def build_drug_report(name, prom, expr, drug, tissue, prism, ctrp, drug_meta):
    auc = _auc(drug, name)
    if len(auc) < 60:
        return None, []
    res, n = scan(prom, auc, tissue, control="tissue+pcs")
    cand = res[res["qvalue"] < 0.1].head(TOP_CANDIDATES)

    # silencing filter
    rows = []
    for c in cand.itertuples():
        g = c.gene
        if g in expr.columns and g in prom.columns:
            sr, sp, sn = silencing(prom[g], expr[g])
            if not np.isnan(sr):
                rows.append((g, c.r, c.qvalue, sr, sp))
    if not rows:
        return None, []
    s = pd.DataFrame(rows, columns=["gene", "resp_r", "resp_q", "sil_r", "sil_p"])
    s["sil_q"] = multipletests(s["sil_p"], method="fdr_bh")[1]
    validated = s[(s["sil_r"] < 0) & (s["sil_q"] < 0.05)].copy()
    n_silenced = len(validated)
    validated = validated.reindex(validated["resp_r"].abs().sort_values(ascending=False).index).head(TOP_LEADS)

    p_auc = _auc(prism, name) if prism is not None and (prism["drug_name"].str.lower() == name.lower()).any() else None
    c_auc = _auc(ctrp, name) if ctrp is not None and (ctrp["drug_name"].str.lower() == name.lower()).any() else None

    leads, lead_records = [], []
    n_replicated = 0
    for v in validated.itertuples():
        g = v.gene
        gr, _, _ = tissue_adj_corr(prom[g], auc, tissue)
        pr, pp, _ = tissue_adj_corr(prom[g], p_auc, tissue)
        cr, cp, _ = tissue_adj_corr(prom[g], c_auc, tissue)
        prism_ok = (not np.isnan(pr)) and np.sign(pr) == np.sign(gr) and pp < 0.05
        ctrp_ok = (not np.isnan(cr)) and np.sign(cr) == np.sign(gr) and cp < 0.05
        screens = ["GDSC"] + (["PRISM"] if prism_ok else []) + (["CTRP"] if ctrp_ok else [])
        if prism_ok or ctrp_ok:
            n_replicated += 1
        rec = per_tissue_recurrence(prom[g], auc, tissue)
        same = sum(1 for t in rec if np.sign(t["rho"]) == np.sign(gr))
        sig = sum(1 for t in rec if np.sign(t["rho"]) == np.sign(gr) and t["p"] < 0.05)
        lead = {
            "gene": g, "direction": "sensitive" if v.resp_r < 0 else "resistant",
            "gdsc_r": round(float(gr), 3),
            "prism_r": None if np.isnan(pr) else round(float(pr), 3),
            "prism_p": None if np.isnan(pp) else round(float(pp), 4),
            "ctrp_r": None if np.isnan(cr) else round(float(cr), 3),
            "ctrp_p": None if np.isnan(cp) else round(float(cp), 5),
            "silencing_r": round(float(v.sil_r), 3),
            "screens": screens, "n_screens": len(screens),
            "evidence_level": evidence_level(len(screens)),
            "replicated": bool(prism_ok or ctrp_ok),
            "recurrence": {"tissues_tested": len(rec), "same_direction": same, "significant": sig},
            "gene_note": GENE_NOTES.get(g, ""),
        }
        leads.append(lead)
        lead_records.append({**lead, "drug": name, "by_tissue": rec})

    leads.sort(key=lambda l: (-l["n_screens"], -abs(l["gdsc_r"])))
    dm = drug_meta.get(name.lower(), {})
    report = {
        "type": "drug", "drug": name,
        "target": dm.get("target"),
        "pathway": dm.get("pathway"),
        "funnel": {"genes_scanned": int(len(res)), "response_candidates": int(len(cand)),
                   "functionally_silenced": int(n_silenced), "externally_replicated": int(n_replicated),
                   "lines": int(n)},
        "leads": leads,
        "hero": leads[0] if leads else None,
        "honesty": HONESTY,
    }
    return report, lead_records


def build_cancer_type_reports(all_leads):
    by_tissue = {}
    for lr in all_leads:
        gr = lr["gdsc_r"]
        for t in lr["by_tissue"]:
            if np.sign(t["rho"]) == np.sign(gr) and t["p"] < 0.10:
                by_tissue.setdefault(t["tissue"], []).append({
                    "drug": lr["drug"], "gene": lr["gene"], "direction": lr["direction"],
                    "tissue_rho": t["rho"], "tissue_p": t["p"],
                    "evidence_level": lr["evidence_level"], "replicated": lr["replicated"],
                    "also_in": [x["tissue"] for x in lr["by_tissue"]
                                if x["tissue"] != t["tissue"] and np.sign(x["rho"]) == np.sign(gr) and x["p"] < 0.10],
                })
    reports = {}
    for tis, leads in by_tissue.items():
        leads.sort(key=lambda l: l["tissue_p"])
        connects = sorted({x for l in leads for x in l["also_in"]})
        reports[tis] = {
            "type": "cancer_type", "cancer_type": tis,
            "n_leads": len(leads), "leads": leads[:30],
            "connects_to": connects,
            "honesty": HONESTY,
        }
    return reports


def main():
    (OUT / "drug").mkdir(parents=True, exist_ok=True)
    (OUT / "cancer_type").mkdir(parents=True, exist_ok=True)
    prom = ccle.load_ccle_promoter_methylation()
    expr = cexpr.load_ccle_expression()
    drug = gdsc.load_gdsc_drug_response()
    tissue = gdsc.load_gdsc_annotations()["tissue"]
    prism = ext.load_prism_drug_response()
    try:
        ctrp = ext.load_ctrp_drug_response()
    except Exception:
        ctrp = None

    import glob
    raw = pd.read_excel(glob.glob(str(ROOT / "data/gdsc/GDSC2_fitted_dose_response_*.xlsx"))[0],
                        usecols=["DRUG_NAME", "PUTATIVE_TARGET", "PATHWAY_NAME"]).drop_duplicates("DRUG_NAME")
    drug_meta = {r.DRUG_NAME.lower(): {
        "target": None if pd.isna(r.PUTATIVE_TARGET) else str(r.PUTATIVE_TARGET),
        "pathway": None if pd.isna(r.PATHWAY_NAME) else str(r.PATHWAY_NAME)} for r in raw.itertuples()}

    present = set(drug["drug_name"].str.lower())
    todo = [d for d in CURATED if d.lower() in present]
    print(f"building reports for {len(todo)} drugs: {todo}")

    all_leads, drug_manifest = [], []
    for i, name in enumerate(todo, 1):
        try:
            rep, recs = build_drug_report(name, prom, expr, drug, tissue, prism, ctrp, drug_meta)
        except Exception as e:
            print(f"  [{i}/{len(todo)}] {name}: ERROR {e}")
            continue
        if rep is None:
            print(f"  [{i}/{len(todo)}] {name}: skipped (insufficient)")
            continue
        (OUT / "drug" / f"{slug(name)}.json").write_text(json.dumps(rep, indent=2))
        all_leads.extend(recs)
        f = rep["funnel"]
        drug_manifest.append({"name": name, "slug": slug(name), "target": rep["target"],
                              "n_leads": len(rep["leads"]), "replicated": f["externally_replicated"]})
        print(f"  [{i}/{len(todo)}] {name}: {len(rep['leads'])} leads, {f['externally_replicated']} replicated -> {slug(name)}.json")

    ct = build_cancer_type_reports(all_leads)
    for tis, rep in ct.items():
        (OUT / "cancer_type" / f"{slug(tis)}.json").write_text(json.dumps(rep, indent=2))
    manifest = {
        "drugs": sorted(drug_manifest, key=lambda d: -d["replicated"]),
        "cancer_types": sorted([{"name": t, "slug": slug(t), "n_leads": r["n_leads"]}
                                for t, r in ct.items()], key=lambda c: -c["n_leads"]),
        "gates": json.loads((ROOT / "site" / "demo-v3" / "cards.json").read_text())["gates"],
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"\nDONE: {len(drug_manifest)} drug reports, {len(ct)} cancer-type reports -> {OUT}")


if __name__ == "__main__":
    main()
