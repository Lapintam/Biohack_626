# NVIDIA Biotech Ventures — the scale / "why us" framing

This is the part of the Strata pitch aimed at the lunch (the real prize). It is grounded
in the *actual measured* kernel cost (see `scripts/scale_benchmark.py` → `site/demo/scale.json`),
not hand-waving. No GPU is claimed to have run — the projection is an honest linear
extrapolation of the measured per-correlation CPU cost.

## The measured numbers (2026-06-04, this laptop)
- Kernel: one Spearman correlation per gene per drug, FDR-corrected (`rank_response_markers`).
- Current GDSC scale: **14,608 genes × 959 cell lines × 286 drugs = 4.18M correlations**, **~746 s** for a full pass on CPU (~178 µs/correlation; 2.6 s/drug).
- Biobank projection (20,000 genes × 1,000 drugs/outcomes × 1,000,000 patients): **~43 CPU-days for a single pass**, ≈5,000× today's compute.

## The argument (say it as a curve, not an assertion)
1. **Generation is cheap and decentralizing.** Any lab/agent can run the screen on its own methylomes — the per-correlation cost is tiny and LLM-driven hypothesis framing is commoditizing.
2. **Verification — re-screening the whole space as data grows, keeping every association cross-checked — is what explodes.** It is O(genes × drugs × cohort) and grows every time anyone contributes data. That is the layer we own.
3. **The kernel is embarrassingly parallel over (gene × drug × cohort)** — independent correlations, the exact shape cuDF / cuML / RAPIDS accelerate. CPU holds today's 4M-correlation panel in minutes; a continuously-updated biobank-scale universe is GPU-bound by construction.
4. **So our unit economics are GPU-native, not GPU-hopeful.** "We don't *want* GPUs; at the scale this works, we structurally require them." The crossover lands on the layer that is our moat.

## The bridge to the platform (the category, not the tool)
Every licensed marker→response association is a **formal, provenance-stamped claim**. Strata is
the working on-ramp to the **Polymer claims universe** (`VISION.md` in this repo): a federated
graph where contributors generate associations on their own compute (Folding@home-style) and a
GPU-bound center keeps the universe related, contradiction-checked, and coherent. Strata proves
the loop end-to-end on one real domain (methylation × drug response); the universe is where it
compounds.

## The lunch line
> "You build the engines for autonomous science. We're the structured memory they write to —
> and the verification layer that keeps it true is RAPIDS-native by construction. Strata is the
> first working instance; the claims universe is where it scales."

## Honesty rails (keep these if a technical VC pokes)
- The biobank number is a **projection on measured CPU cost**, labeled as such — not a benchmarked GPU run.
- The marker screen **selects on the outcome**, so the median-split subgroup p is **not** independent validation; credibility comes from the screen recovering the known CDKN2A/9p21 → palbociclib positive control. Novel hits are **hypotheses**.
- Methylation here is **gene-level imputed** GDSC (Cell Model Passports), not CpG-island resolution — stated plainly.
