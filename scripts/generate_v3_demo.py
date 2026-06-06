"""Generate site/demo-v3/cards.json — data for the three-gates / Track B front end.

Assembles: the three structural gates (narrative), the discovery funnel, and the
externally-replicated hero cards (with per-tissue cross-indication recurrence).

Run: PYTHONPATH=. uv run python scripts/generate_v3_demo.py
"""

from __future__ import annotations

import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from strata.config import ROOT
from strata.data import gdsc
from strata.data import ccle_methylation as ccle

OUT = ROOT / "site" / "demo-v3"

# short, honest gene annotations (why each is "obscure" / not a canonical target)
GENE_NOTES = {
    "ID2": "Inhibitor of DNA-binding 2 — a developmental/differentiation transcription factor. NOT a canonical PARP/DDR gene; no expert would prioritize it for a PARP inhibitor.",
    "TRIM32": "TRIM-family E3 ubiquitin ligase. Not a known PARP-response gene.",
    "OTULINL": "OTULIN-like (deubiquitinase-related). Essentially uncharacterized in drug response.",
    "TSPAN13": "Tetraspanin-13, a membrane scaffold. No prior link to temozolomide.",
    "N4BP3": "NEDD4-binding protein 3. Not a known MEK-pathway gene.",
}
HERO_DRUG_GDSC = {"Olaparib": "Olaparib", "Temozolomide": "Temozolomide", "Trametinib": "Trametinib"}


def per_tissue_recurrence(meth_g, auc, tissue, min_n=10):
    d = pd.DataFrame({"m": meth_g, "a": auc, "t": tissue}).dropna()
    rows = []
    for tis, sub in d.groupby("t"):
        if len(sub) < min_n or sub["m"].std() == 0:
            continue
        rho, p = spearmanr(sub["m"], sub["a"])
        rows.append({"tissue": tis, "n": int(len(sub)), "rho": round(float(rho), 3), "p": round(float(p), 4)})
    rows.sort(key=lambda r: r["rho"])
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    prom = ccle.load_ccle_promoter_methylation()
    drug = gdsc.load_gdsc_drug_response()
    tissue = gdsc.load_gdsc_annotations()["tissue"]
    rep = pd.read_csv(ROOT / "outputs" / "external_replication.csv")

    hero = rep[rep["prism_replicates"]].copy()
    cards = []
    for r in hero.itertuples():
        gdsc_name = HERO_DRUG_GDSC.get(r.drug, r.drug)
        auc = (drug[drug["drug_name"].str.lower() == gdsc_name.lower()]
               .drop_duplicates("COSMIC_ID").set_index("COSMIC_ID")["auc"].astype(float))
        rec = per_tissue_recurrence(prom[r.gene], auc, tissue) if r.gene in prom.columns else []
        same_dir = sum(1 for t in rec if np.sign(t["rho"]) == np.sign(r.gdsc_r))
        sig_dir = sum(1 for t in rec if np.sign(t["rho"]) == np.sign(r.gdsc_r) and t["p"] < 0.05)
        screens = ["GDSC", "PRISM"] + (["CTRP"] if r.ctrp_replicates else [])
        cards.append({
            "drug": r.drug, "gene": r.gene,
            "direction": "sensitive" if r.resp_r < 0 else "resistant",
            "gdsc_r": round(r.gdsc_r, 3), "prism_r": round(r.prism_r, 3), "prism_p": round(r.prism_p, 4),
            "ctrp_r": None if pd.isna(r.ctrp_r) else round(r.ctrp_r, 3),
            "ctrp_p": None if pd.isna(r.ctrp_p) else round(r.ctrp_p, 5),
            "silencing_r": round(r.silencing_r, 3),
            "screens": screens, "n_screens": len(screens),
            "evidence_level": "L5" if len(screens) == 3 else "L4",
            "gene_note": GENE_NOTES.get(r.gene, ""),
            "recurrence": {"tissues_tested": len(rec), "same_direction": same_dir,
                            "significant": sig_dir, "by_tissue": rec[:8]},
        })
    cards.sort(key=lambda c: (-c["n_screens"], c["gdsc_r"]))

    payload = {
        "thesis": "Start from molecular state, not disease taxonomy. Surface opportunities an expert is STRUCTURALLY unable to reach.",
        "gates": [
            {"id": "taxonomy", "title": "The taxonomy prior",
             "expert": "Trained, funded, and regulated by organ/histology. The unit of analysis is the indication.",
             "agent": "Groups cell lines by molecular state across ~30 cancer types; the organ is a variable to adjust away, never a boundary.",
             "evidence": "Hero leads hold tissue-agnostically across unrelated indications."},
            {"id": "prestige", "title": "The known-gene prior",
             "expert": "Evaluates hypotheses about famous genes (KRAS, TP53). Won't fund a hypothesis about a gene they've never heard of.",
             "agent": "Ranks all 17,181 genes with zero prestige prior; surfaces the functionally-validated long tail.",
             "evidence": "Hero gene ID2 is a developmental TF, not a PARP/DDR gene — invisible to a domain prior."},
            {"id": "control", "title": "Multivariate confound control",
             "expert": "Reasons over one or two variables; cannot residualize the whole genome against latent axes in their head.",
             "agent": "Residualizes response + all genes against tissue + latent methylome PCs (SVA/RUV-style) simultaneously.",
             "evidence": "Caught a 9p21 copy-number artifact masquerading as methylation — a confound naive analysis reports as a finding."},
        ],
        "funnel": [
            {"stage": "All promoter genes", "n": 17181, "note": "unbiased, genome-wide"},
            {"stage": "Global-axis-controlled response candidates", "n": "~300–500/drug", "note": "tissue + methylome PCs removed"},
            {"stage": "Functionally silenced (RNA-seq)", "n": "~160/drug", "note": "promoter-meth ↓ own expression"},
            {"stage": "Externally replicated", "n": len(cards), "note": "PRISM (+CTRP); same direction, p<0.05"},
        ],
        "hero_cards": cards,
        "honesty": [
            "Effect sizes are modest (single-gene methylation→response is weak/diffuse). The claim is directional concordance across independent screens + functional silencing — not a large effect.",
            "These are cell-line results: defensible novel HYPOTHESES (evidence ladder ~L5), not clinical claims.",
            "The agent expands the hypothesis space the expert then adjudicates. It removes the structural gate; it does not replace expert judgment.",
        ],
    }
    (OUT / "cards.json").write_text(json.dumps(payload, indent=2))
    print(f"wrote {OUT/'cards.json'}  | hero cards: {len(cards)}")
    for c in cards:
        print(f"  {c['drug']:<13} {c['gene']:<10} {c['n_screens']}-screen {c['evidence_level']}  "
              f"recurs same-dir {c['recurrence']['same_direction']}/{c['recurrence']['tissues_tested']} tissues "
              f"({c['recurrence']['significant']} sig)")


if __name__ == "__main__":
    main()
