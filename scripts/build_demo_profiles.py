"""Build demo 'patient' methylation profiles for the annotation demo.

Picks a few cell lines spanning the sensitivity spectrum (by measured AUC) and
emits each as a patient-style profile {gene: beta} over the annotation-model genes,
WITH ground-truth measured AUC (so the demo can show prediction vs truth).
Clearly labeled as illustrative in-cohort samples (a real patient tumor/PBMC is leakage-free, added later).

Output: api/reports/demo_profiles/<id>.json + api/reports/demo_profiles/index.json

Run: PYTHONPATH=. uv run python scripts/build_demo_profiles.py
"""

from __future__ import annotations

import json
import numpy as np

from strata.config import ROOT
from strata.data import gdsc
from strata.data import ccle_methylation as ccle

OUT = ROOT / "api" / "reports" / "demo_profiles"
GT_DRUGS = ["Olaparib", "Trametinib", "Temozolomide", "Palbociclib"]


def _auc(drug, name):
    s = (drug[drug["drug_name"].str.lower() == name.lower()]
         .drop_duplicates("COSMIC_ID").set_index("COSMIC_ID")["auc"])
    return s[~s.index.duplicated(keep="first")].astype(float)


def slug(s):
    import re
    return re.sub(r"[^a-z0-9]+", "-", str(s).lower()).strip("-")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    model = json.loads((ROOT / "api" / "reports" / "annotation_model.json").read_text())
    genes = list(model["cohort"].keys())
    prom = ccle.load_ccle_promoter_methylation()
    drug = gdsc.load_gdsc_drug_response()
    ann = gdsc.load_gdsc_annotations()
    ola = _auc(drug, "Olaparib")

    # candidate lines: good model-gene coverage + have olaparib AUC
    cov = prom[genes].notna().mean(axis=1)
    cand = [c for c in ola.index if c in prom.index and cov.get(c, 0) > 0.8]
    ola_c = ola.loc[cand].sort_values()
    picks = list(ola_c.head(2).index) + list(ola_c.tail(2).index)   # 2 most sensitive + 2 most resistant

    aucs = {d: _auc(drug, d) for d in GT_DRUGS}
    index = []
    for c in picks:
        name = ann["cell_line_name"].get(c, c)
        row = prom.loc[c]
        profile = {g: round(float(row[g]), 4) for g in genes if g in prom.columns and row[g] == row[g]}
        gt = {d: (None if c not in aucs[d].index else round(float(aucs[d].loc[c]), 3)) for d in GT_DRUGS}
        pid = slug(name)
        rec = {
            "id": pid, "label": f"{name} (illustrative cell-line sample)",
            "name": str(name), "cosmic_id": str(c),
            "tissue": str(ann["tissue"].get(c, "")),
            "source": "cell_line_in_cohort",
            "note": "Illustrative: this line was in the discovery cohort (training leakage). A real patient tumor/PBMC sample is leakage-free.",
            "n_genes": len(profile),
            "ground_truth_auc": gt,
            "profile": profile,
        }
        (OUT / f"{pid}.json").write_text(json.dumps(rec))
        index.append({"id": pid, "label": rec["label"], "tissue": rec["tissue"],
                      "source": rec["source"], "ground_truth_auc": gt})
        print(f"  {name:<14} cosmic={c} genes={len(profile)} olaparib_AUC={gt['Olaparib']}")
    (OUT / "index.json").write_text(json.dumps({"profiles": index}))
    print(f"wrote {len(index)} demo profiles -> {OUT}")


if __name__ == "__main__":
    main()
