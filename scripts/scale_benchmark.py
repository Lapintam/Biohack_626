"""Honest scale benchmark for the marker->response kernel.

The core Strata operation (`rank_response_markers`) is a Spearman-correlation
scan: every gene's methylation vector vs a drug's response vector. It is
embarrassingly parallel over (gene x drug) and is the kernel that becomes
GPU-bound at biobank scale. This script times it at the current GDSC scale and
writes an HONEST extrapolation (measured small scale -> projected large scale,
linear in the number of correlations) to ``site/demo/scale.json``.

No GPU is used or claimed here. The projection is a back-of-envelope on the
real measured per-correlation cost; it is labelled as such on the page.

Run: ``PYTHONPATH=. uv run python scripts/scale_benchmark.py``
"""

from __future__ import annotations

import json
import time

from strata.config import OUTPUT_DIR
from strata.data import gdsc
from strata.engine.associate import rank_response_markers

DEMO_DIR = (OUTPUT_DIR.parent / "site" / "demo")


def _time_scan(meth, drug, drug_name: str):
    auc = (
        drug[drug["drug_name"].str.lower() == drug_name.lower()]
        .set_index("COSMIC_ID")["auc"]
    )
    auc = auc[~auc.index.duplicated(keep="first")]
    common = meth.index.intersection(auc.index)
    m, a = meth.loc[common], auc.loc[common]
    t0 = time.perf_counter()
    rank_response_markers(m, a)
    dt = time.perf_counter() - t0
    return dt, m.shape[0], m.shape[1]


def main() -> None:
    meth = gdsc.load_gdsc_methylation()
    drug = gdsc.load_gdsc_drug_response()
    n_drugs = int(drug["drug_name"].nunique())

    # Time one representative drug scan (median of a few for stability).
    times = []
    for d in ("Palbociclib", "Luminespib", "Vorinostat"):
        dt, n_lines, n_genes = _time_scan(meth, drug, d)
        times.append(dt)
    per_drug = sorted(times)[len(times) // 2]  # median

    corr_per_drug = n_genes  # one Spearman per gene
    full_gdsc_corr = corr_per_drug * n_drugs
    full_gdsc_secs = per_drug * n_drugs
    per_corr_us = (per_drug / corr_per_drug) * 1e6  # microseconds per correlation

    # Honest biobank-scale extrapolation (linear in #correlations).
    # ~20k genes x ~1000 drugs/outcomes x cohort scaling. The correlation count
    # is what grows; the per-correlation cost also grows ~linearly with cohort n,
    # so we scale by (genes/14608) * (drugs/286) * (cohort_n / current_n).
    current_n = n_lines
    biobank = {
        "genes": 20000,
        "drugs_or_outcomes": 1000,
        "cohort_n": 1_000_000,
    }
    scale_factor = (
        (biobank["genes"] / n_genes)
        * (biobank["drugs_or_outcomes"] / n_drugs)
        * (biobank["cohort_n"] / current_n)
    )
    biobank_corr = full_gdsc_corr * (biobank["genes"] / n_genes) * (
        biobank["drugs_or_outcomes"] / n_drugs
    )
    biobank_cpu_secs = full_gdsc_secs * scale_factor

    record = {
        "kernel": "per-gene Spearman(methylation, drug response), FDR-corrected",
        "current": {
            "genes": int(n_genes),
            "cell_lines": int(n_lines),
            "drugs": int(n_drugs),
            "seconds_per_drug": round(per_drug, 3),
            "microseconds_per_correlation": round(per_corr_us, 2),
            "full_panel_correlations": int(full_gdsc_corr),
            "full_panel_seconds_cpu": round(full_gdsc_secs, 1),
        },
        "biobank_projection": {
            "assumption": biobank,
            "correlations": float(f"{biobank_corr:.3e}"),
            "cpu_seconds": float(f"{biobank_cpu_secs:.3e}"),
            "cpu_hours": round(biobank_cpu_secs / 3600, 1),
            "cpu_days": round(biobank_cpu_secs / 86400, 1),
            "scale_factor_vs_gdsc": float(f"{scale_factor:.3e}"),
        },
        "note": (
            "Projection is linear in the number of correlations on the measured "
            "per-correlation CPU cost; no GPU run is implied. The kernel is "
            "embarrassingly parallel over (gene x drug x cohort), which is exactly "
            "the shape RAPIDS / cuDF / cuML accelerate."
        ),
    }

    DEMO_DIR.mkdir(parents=True, exist_ok=True)
    (DEMO_DIR / "scale.json").write_text(json.dumps(record, indent=2))
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
