# MethylGKB · Strata

An AI agent that does **methylation-marker → drug-response subgroup discovery**: it scans a methylome against an outcome, surfaces the marker-defined responder subgroup, recovers known biology as a positive control, and writes it up. Built for the NYC Tech Week *Agents in Biomedical Science* hackathon (Pfizer "Commercial Development Discovery" challenge).

**Live demo:** https://methylkb.vercel.app/strata.html

---

## Read this first
- **`docs/SYSTEM_MAP.md`** — the precise, source-derived map of **what we have**: every data asset (shapes + provenance), every inference rule (exact parameters), the flow, and the epistemic rails. *Start here to understand the system.*
- **`docs/OPPORTUNITY_GRAPH_PLAN.md`** — the **v2 plan under review**: where this goes (cross-indication / molecular-state opportunity graph). Proposal stage, not built.

## Directory map
| Path | What it is |
|---|---|
| `strata/` | the Python package — `data/` loaders, `engine/` (the inference rules), `agent/` (the Claude tool-use loop + scripted fallback), `figures.py`, `cli.py` |
| `site/` | the deployed static demo (vanilla HTML/CSS/JS; replays pre-generated runs in `site/demo/`) |
| `tests/` | pytest suite (24 tests) |
| `scripts/` | run + diagnostic scripts (`run_demo.py`, `scale_benchmark.py`, `tcga_vignette.py`, `mgmt_gate.py`, `diag_positive_control.py`, `scan_novel_hits.py`) |
| `docs/` | canonical docs — `SYSTEM_MAP.md`, `OPPORTUNITY_GRAPH_PLAN.md`, the explainer/system-map HTML (PDF sources), and `superpowers/specs` + `superpowers/plans` (the original build spec + plan) |
| `notes/` | provenance + results — `GDSC_SCHEMA`, `TCGA_SCHEMA`, `MGMT_GATE`, `DEMO_HITS`, `NOVEL_HITS`, `NVIDIA_PITCH` |
| `data/`, `outputs/` | **gitignored** — downloaded datasets + generated figures/briefs |
| `_archive/` | **HISTORICAL, not canonical** — see below |

## Quick start
```bash
uv run pytest -q                                          # tests
uv run python -m strata.cli --drug Palbociclib            # deterministic discovery pipeline
uv run python -m strata.agent.loop                        # the agent (scripted w/o key; live w/ ANTHROPIC_API_KEY in .env)
PYTHONPATH=. uv run python scripts/run_demo.py            # regenerate site/demo artifacts
uv run python -m http.server 8000 --directory site        # view the site locally
vercel deploy --prod --yes --cwd site                     # redeploy to production
```

## Status
- **v1 (shipped + deployed):** marker→drug-response discovery on GDSC; positive control CDKN2A/MTAP → palbociclib; live Claude agent (Palbociclib) with a crash-proof scripted fallback; the demo site.
- **v2 (planned, under review):** cross-indication / molecular-state opportunity graph — see `docs/OPPORTUNITY_GRAPH_PLAN.md`.

## `_archive/` — do not treat as current
Pre-pivot material kept for reference only, **superseded by the Strata work**:
- The original **Polymer claims-universe** pitch (`VISION.md`, `DEMO_DAY_PLAN.md`, `PITCH_FRAME_DISCUSSION.md`, `NVIDIA_INTEGRATION_PLAN.md`) — `VISION.md` is the "where this goes" horizon referenced by the v2 plan.
- The pre-pivot **Methylome Interpreter** lead (`METHYLOME_INTERPRETER*.md`) and the **idea cultivation bed** (`IDEAS_COMPENDIUM.md`).
- Prior **pharmaco-epigenomics research** (`PHARMACO_EPIGENOMICS_RESEARCH.md`, `RESEARCH_AUDIT.md`) — dated, but worth mining for the v2 prior-art pass (PGx landscape, methylation-signal verification).

*Future instances: read `docs/SYSTEM_MAP.md`, not `_archive/`.*
