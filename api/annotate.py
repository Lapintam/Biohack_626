"""Patient-methylation annotation scorer — self-contained copy for the API image.

Pure copy of the dependency-light scorer from ``strata/annotate.py`` (json + bisect
only) so the Fly image (which copies only ``api/``) never has to import the heavier
``strata`` package. Given a patient promoter-methylation profile {gene: beta}, score
it against the validated leads in ``api/reports/annotation_model.json``.

Research use only.
"""

from __future__ import annotations

import bisect
import json
from functools import lru_cache
from pathlib import Path

MODEL_PATH = Path(__file__).parent / "reports" / "annotation_model.json"
THRESHOLD = 0.15      # |score| above which we call sensitive/resistant
MIN_LEADS = 2         # min scorable leads for a drug to be called


def load_model(path=MODEL_PATH):
    return json.loads(Path(path).read_text())


@lru_cache(maxsize=1)
def load_annotation_model() -> dict:
    """Load the annotation model once (lead defs + cohort distributions)."""
    return load_model()


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
