"""Build the annotation model: lead definitions + per-gene cohort distributions,
baked into a single JSON so the STATELESS API can score a patient methylation
profile without the raw data.

Output: api/reports/annotation_model.json
  { "drugs": { drug: [ {gene, direction, gdsc_r, n_screens, evidence_level, silencing_r} ] },
    "cohort": { gene: [sorted promoter beta values across cell lines] },
    "n_drugs", "n_genes" }

Run: PYTHONPATH=. uv run python scripts/build_annotation_assets.py
"""

from __future__ import annotations

import glob
import json
import numpy as np

from strata.config import ROOT
from strata.data import ccle_methylation as ccle

OUT = ROOT / "api" / "reports" / "annotation_model.json"


def main():
    drugs = {}
    genes = set()
    for f in glob.glob(str(ROOT / "api" / "reports" / "drug" / "*.json")):
        rep = json.loads(open(f).read())
        leads = []
        for l in rep["leads"]:
            leads.append({
                "gene": l["gene"], "direction": l["direction"],
                "gdsc_r": l["gdsc_r"], "n_screens": l["n_screens"],
                "evidence_level": l["evidence_level"], "silencing_r": l.get("silencing_r"),
            })
            genes.add(l["gene"])
        if leads:
            drugs[rep["drug"]] = leads

    prom = ccle.load_ccle_promoter_methylation()
    cohort = {}
    for g in sorted(genes):
        if g in prom.columns:
            vals = prom[g].dropna().values
            if len(vals) >= 30:
                cohort[g] = [round(float(v), 4) for v in np.sort(vals)]

    # drop leads whose gene has no cohort distribution (can't be scored)
    for d in list(drugs):
        drugs[d] = [l for l in drugs[d] if l["gene"] in cohort]
        if not drugs[d]:
            del drugs[d]

    model = {"drugs": drugs, "cohort": cohort, "n_drugs": len(drugs), "n_genes": len(cohort)}
    OUT.write_text(json.dumps(model))
    print(f"wrote {OUT}")
    print(f"  drugs: {len(drugs)} | scorable genes: {len(cohort)} | size: {OUT.stat().st_size/1e6:.1f} MB")
    print(f"  example: Olaparib has {len(drugs.get('Olaparib', []))} scorable leads")


if __name__ == "__main__":
    main()
