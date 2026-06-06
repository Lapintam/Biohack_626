"""Patient-methylation annotation scorer.

Given a patient's promoter methylation profile {gene: beta}, score it against the
validated leads: for each lead, percentile the patient's beta within the cell-line
cohort, flip by the lead's direction, aggregate per drug -> predicted sensitive /
resistant / neutral + the driver genes.

Pure + stateless: needs only the annotation_model.json (lead defs + cohort
distributions), so it runs inside the API or locally. Research-use only.
"""

from __future__ import annotations

import bisect
import json
from pathlib import Path

DEFAULT_MODEL = Path(__file__).resolve().parent.parent / "api" / "reports" / "annotation_model.json"
THRESHOLD = 0.15      # |score| above which we call sensitive/resistant
MIN_LEADS = 2         # min scorable leads for a drug to be called


def load_model(path=DEFAULT_MODEL):
    return json.loads(Path(path).read_text())


def _percentile(value, sorted_vals):
    """Fraction of the cohort below `value` (0..1)."""
    n = len(sorted_vals)
    if n == 0:
        return None
    return bisect.bisect_left(sorted_vals, value) / n


def annotate_profile(profile: dict, model: dict, threshold: float = THRESHOLD,
                     min_leads: int = MIN_LEADS) -> dict:
    """profile: {gene: beta}. Returns ranked per-drug predictions + drivers."""
    cohort = model["cohort"]
    predictions = []
    for drug, leads in model["drugs"].items():
        contribs = []
        for lead in leads:
            g = lead["gene"]
            if g in profile and g in cohort:
                pct = _percentile(profile[g], cohort[g])
                if pct is None:
                    continue
                signal = (pct - 0.5) * 2.0                       # [-1,1], + = more methylated than cohort
                sens = signal * (1 if lead["direction"] == "sensitive" else -1)
                weight = abs(lead["gdsc_r"]) * (1 + 0.5 * (lead["n_screens"] - 1))
                contribs.append({"gene": g, "sens_signal": round(sens, 3),
                                 "weight": round(weight, 3), "patient_pct": round(pct, 3),
                                 "direction": lead["direction"], "evidence_level": lead["evidence_level"]})
        if len(contribs) < min_leads:
            continue
        wsum = sum(c["weight"] for c in contribs)
        score = sum(c["sens_signal"] * c["weight"] for c in contribs) / wsum if wsum else 0.0
        call = "sensitive" if score > threshold else ("resistant" if score < -threshold else "neutral")
        drivers = sorted(contribs, key=lambda c: -abs(c["sens_signal"] * c["weight"]))[:4]
        predictions.append({
            "drug": drug, "score": round(score, 3), "call": call,
            "n_leads_scored": len(contribs),
            "drivers": [{"gene": d["gene"], "patient_pct": d["patient_pct"],
                         "contribution": round(d["sens_signal"] * d["weight"], 3),
                         "evidence_level": d["evidence_level"]} for d in drivers],
        })
    predictions.sort(key=lambda p: p["score"])   # most sensitive (most negative... no: + = sensitive) -> sort desc
    predictions.sort(key=lambda p: -p["score"])  # highest sensitivity score first
    return {
        "n_drugs_scored": len(predictions),
        "predictions": predictions,
        "honesty": [
            "Research use only. Predictions are cell-line-derived hypotheses, not clinical recommendations.",
            "Each call aggregates weak single-gene methylation signals; treat as enrichment/triage, not a deterministic result.",
            "Scored against the cell-line cohort distribution; a genuinely novel sample (e.g. a patient tumor) carries no training leakage, an in-cohort cell line does.",
        ],
    }


def profile_from_cohort_line(cosmic_id, model_genes, prom):
    """Helper: build a mock-patient profile from a cell line's promoter methylation row."""
    row = prom.loc[str(cosmic_id)]
    return {g: float(row[g]) for g in model_genes if g in prom.columns and not _isnan(row[g])}


def _isnan(x):
    return x != x


if __name__ == "__main__":
    import sys
    from strata.data import ccle_methylation as ccle
    model = load_model()
    prom = ccle.load_ccle_promoter_methylation()
    cosmic = sys.argv[1] if len(sys.argv) > 1 else str(prom.index[0])
    profile = profile_from_cohort_line(cosmic, set(model["cohort"]), prom)
    out = annotate_profile(profile, model)
    print(f"mock patient = cell line {cosmic} | genes in profile: {len(profile)}")
    print(f"drugs scored: {out['n_drugs_scored']}")
    for p in out["predictions"][:8]:
        drv = ", ".join(d["gene"] for d in p["drivers"])
        print(f"  {p['drug']:<14} {p['call']:<10} score={p['score']:+.2f}  drivers: {drv}")
