"""Generate the Strata demo artifacts the static page reads.

For each demo drug:
  1. ``strata.cli.run`` -> writes 3 figures to ``outputs/`` + returns (markers_df, assoc)
  2. ``strata.agent.loop.run_discovery_agent`` -> transcript + brief (scripted offline)
  3. Copy the 3 figures into ``site/demo/`` and write a per-drug JSON record.

Finally writes ``site/demo/index.json``.

Run:  uv run python scripts/run_demo.py
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

from strata import cli
from strata.agent.loop import run_discovery_agent
from strata.config import OUTPUT_DIR

ROOT = Path(__file__).resolve().parent.parent
SITE_DEMO = ROOT / "site" / "demo"

DRUGS = ["Palbociclib", "Luminespib"]


def _markers_payload(markers_df, n: int = 12) -> list[dict]:
    cols = [c for c in ("gene", "rho", "qvalue") if c in markers_df.columns]
    head = markers_df.head(n)[cols]
    out: list[dict] = []
    for _, row in head.iterrows():
        out.append(
            {
                "gene": str(row["gene"]),
                "rho": float(row["rho"]),
                "qvalue": float(row["qvalue"]),
            }
        )
    return out


def build_run(drug: str) -> dict:
    print(f"\n=== Building demo run: {drug} ===")
    # cli.run returns the FULL ranked significant markers df (all q<0.1 rows)
    markers_df, _assoc = cli.run(drug)
    agent = run_discovery_agent(drug)

    figures = {}
    for kind in ("umap", "auc", "heatmap"):
        fname = f"{drug}_{kind}.png"
        src = OUTPUT_DIR / fname
        dst = SITE_DEMO / fname
        shutil.copyfile(src, dst)
        figures[kind] = fname

    is_positive_control = drug == "Palbociclib"

    # Locate CDKN2A in the full ranked markers df for the positive control
    control_check = None
    if is_positive_control and "CDKN2A" in set(markers_df["gene"]):
        row = markers_df.reset_index(drop=True)
        idx = row.index[row["gene"] == "CDKN2A"][0]
        control_check = {
            "gene": "CDKN2A",
            "rank": int(idx) + 1,
            "total": int(len(row)),
            "rho": float(row.loc[idx, "rho"]),
            "qvalue": float(row.loc[idx, "qvalue"]),
            "neighbor": "MTAP",
        }

    # Cell-line count: stored in the "n" column of markers_df
    n_cell_lines = int(markers_df["n"].iloc[0]) if len(markers_df) else 0
    n_significant = int(len(markers_df))

    return {
        "drug": drug,
        "goal": agent.get("goal", ""),
        "is_positive_control": is_positive_control,
        "mode": agent.get("mode", "scripted"),
        "transcript": agent.get("transcript", []),
        "brief": agent.get("brief", ""),
        "top_markers": _markers_payload(markers_df, 12),
        "figures": figures,
        "control_check": control_check,
        "n_significant": n_significant,
        "n": n_cell_lines,
    }


def main() -> None:
    SITE_DEMO.mkdir(parents=True, exist_ok=True)
    for drug in DRUGS:
        record = build_run(drug)
        out = SITE_DEMO / f"{drug}.json"
        out.write_text(json.dumps(record, indent=2))
        print(f"[wrote] {out}  ({len(record['transcript'])} transcript steps, "
              f"{len(record['top_markers'])} markers)")

    index = {"runs": DRUGS, "generated": "rung-A scripted"}
    (SITE_DEMO / "index.json").write_text(json.dumps(index, indent=2))
    print(f"[wrote] {SITE_DEMO / 'index.json'}")
    print("\nDone. site/demo/ now contains:")
    for f in sorted(SITE_DEMO.iterdir()):
        print(f"  {f.name}  ({f.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
