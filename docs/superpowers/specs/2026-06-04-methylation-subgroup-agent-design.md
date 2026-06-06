# Methylation Subgroup-Discovery Agent — Design

**Date:** 2026-06-04
**Event:** NYC Tech Week hackathon — *Agents in Biomedical Science*. **Sponsor challenge:** Pfizer — *Commercial Development Discovery*. Live demo Saturday (2026-06-06). **Team:** solo build.
**Project name:** **MethylGKB** (matches the existing site brand). The subgroup-discovery engine described here is the **Strata** submodule.
**Relationship to the existing demo:** the Pfizer/Strata work is the new *lead*; the abstract Polymer claims-universe plan (`VISION.md`, `DEMO_DAY_PLAN.md`, `PITCH_FRAME_DISCUSSION.md`) and the current `site/` are **kept and built alongside, not removed** — they supply the "where this goes" horizon and the output surface (see §9). We trim later, never up front.

---

## 1. The Pfizer gap we answer

Pfizer's challenge, verbatim intent: commercial development can only evaluate a handful of opportunities over months; the highest-value treatment opportunities are *"the patient subgroup hidden inside a heterogeneous indication"* — and *"by definition the ones the existing taxonomy is least equipped to find."* They want to test whether an **AI-native approach reveals opportunities the expert-gated model is structurally incapable of finding.**

**Our answer:** peripheral-blood / tumor **DNA methylation** is an underexploited axis for resolving hidden **drug-response subgroups**. The paradigm is real — e.g. **CDKN2A/p16 methylation silences the CDK4/6 brake, predicting response to CDK4/6 inhibitors like Pfizer's palbociclib (Ibrance)** — and the bet is that *many more are buried*, undiscoverable by expert review because no human is systematically clustering methylomes against response. An agent can.

---

## 2. The pitch, in one breath

> *"Methylation can define who responds to a drug — silence p16/CDKN2A and a tumor becomes dependent on CDK4/6, the exact target of Pfizer's own Ibrance (palbociclib). We point an agent at an indication's methylation data; it recovers that CDKN2A → palbociclib responder axis cold, unsupervised, as a positive control — then finds a **novel** methylation-defined responder subgroup nobody has named, grounds it in the biology, and writes it up as a commercial-development opportunity. The one dataset that would join blood methylation to drug response sits locked in controlled access — which is exactly why these subgroups stay hidden. Here's the method, proven on what's open, ready the day that data unlocks."*

Two rhetorical jewels: **CDKN2A → palbociclib recovery** (a real, p=6×10⁻²⁷ signal, recovered unsupervised — and it's Pfizer's own drug) and **the locked dataset** (the literal embodiment of Pfizer's "least equipped to find").

> **Gate verdict (2026-06-04):** the originally-planned MGMT → temozolomide control is **null in gene-level GDSC methylation** (rho=−0.018, p=0.58) — promoter-CpG-island silencing is averaged out at gene level. Replaced by the stronger, data-backed **CDKN2A → palbociclib** (rho=−0.336, p=6×10⁻²⁷, n=968). MGMT is dropped from the narrative.

---

## 3. Data plan (verified by data scout, 2026-06-04)

No open dataset pairs peripheral-blood methylation + real drug response. The ideal (BeatAML: AML blood methylation + ex-vivo AUC, 122 drugs) is **dbGaP phs001657 controlled-access — not obtainable this week.** That constraint is turned into the pitch (§2). What we build on:

| Role | Dataset | What it gives | Access |
|---|---|---|---|
| **Engine** | **GDSC** (now via Cell Model Passports / Sanger; `cancerrxgene.org` decommissioned) | ~1,055 cell lines, **gene-level imputed methylation** (~14.6k gene-symbol columns; NOT CpG islands — see gate verdict) **paired with drug IC50/AUC across 286 compounds** incl. palbociclib. Join methylation↔drug via SANGER_MODEL_ID→COSMIC_ID. | Open, processed matrices |
| **Patient vignette** | **TCGA-LAML** (GDC `TCGA-LAML`) | ~194 patients, 450K methylation, **overall survival** + full driver mutations (FLT3/NPM1/DNMT3A/IDH1-2/TET2). Blood/marrow leukemic blasts. | Open via GDC |
| **Honesty slide** | **BeatAML methylation** (`phs001657`) | The perfect substrate (methylation + per-patient drug AUC). | **Locked** — used as the narrative, not built on |

**Positive-control verification was the Day-1 gate (§6) — now resolved.** MGMT→temozolomide came back null at gene-level resolution; the established **CDKN2A→palbociclib** association (rho=−0.336, p=6×10⁻²⁷, n=968) is the data-backed positive control the demo recovers. (See gate verdict in §2.)

---

## 4. Demo run-of-show (~3 min, target rung)

1. **(0:00) The gap.** Pfizer's own words: the best opportunities are the ones the taxonomy can't see. CDKN2A → palbociclib is proof that methylation can define who responds to a drug. Why aren't we mining methylomes for more?
2. **(0:30) Point the agent at GDSC.** Goal: *"find a methylation-defined drug-response subgroup."* It selects features, clusters, evaluates.
3. **(1:00) Positive control — it recovers CDKN2A → palbociclib cold.** The agent surfaces the CDKN2A-methylated subgroup with differential palbociclib response, unsupervised. *"It found that methylation predicts response to your own drug, on its own."*
4. **(1:30) The novel hit.** The agent surfaces a **second**, un-named methylation-defined responder subgroup for another compound, annotates the driver genes/pathway, and proposes the mechanistic hypothesis.
5. **(2:00) The real-patient anchor.** TCGA-LAML: the agent finds a methylation subgroup with differential **survival** + mutation enrichment. "Real patients, real myeloid, real blood-adjacent."
6. **(2:30) The opportunity brief + the close.** The agent emits a one-page commercial-development brief in Pfizer's success language (portfolio ROI, reduced identification time/cost). Close on the locked-data point + scale ("give an agent GPU and the controlled-access methylomes, this runs across every indication").

A pre-recorded backup video at the highest working rung is mandatory (Friday).

---

## 5. Architecture — five components

1. **Data layer** — downloaders + loaders producing tidy matrices: GDSC methylation β (CpG-island × cell line), GDSC drug response (cell line × drug → AUC/IC50), GDSC line annotations; TCGA-LAML β-matrix + clinical/mutation table. Cached locally; deterministic.
2. **Discovery engine** *(deterministic Python — the floor)* — variable-feature selection → dimensionality reduction → consensus/k clustering → subgroups. Pure functions, each independently testable.
3. **Association** — per-subgroup test of the response variable (drug AUC via Mann-Whitney/Kruskal-Wallis; survival via log-rank/KM). Returns effect size + p + the responder stratum, scored against explicit ++/+/neutral/−/−− criteria.
4. **Annotation / interpretation** — driver-CpG differential methylation → gene/pathway mapping → mechanistic hypothesis. *(Stretch C: ground this by querying the existing PolymerGenomicsAPI / methylation knowledge base instead of a bare lookup.)*
5. **Agent + output** — Claude API tool-use loop where components 2–4 are the agent's tools. The agent drives the protocol, judges signal-vs-noise, iterates, and writes the opportunity brief. Output rendered into the existing static `site/` (reuse `report.html` / `explore.html`) plus figures: UMAP colored by response, driver-CpG heatmap, KM curve, differential-AUC boxplot.

**Isolation principle:** the deterministic engine (2–4) is a clean library with no LLM dependency — fully runnable and testable on its own; the agent (5) only *orchestrates* it. This guarantees a working demo even if the agent layer misbehaves, and keeps each unit small enough to reason about.

---

## 6. De-risking ladder (ship at whatever rung we reach)

| Rung | Deliverable | State if we stop here |
|---|---|---|
| **Day-1 gate** | Download GDSC; confirm MGMT island + temozolomide present; reproduce the MGMT→TMZ association by hand. | Go/no-go on the positive-control story; pick fallback if needed. |
| **A (floor)** | Deterministic engine end-to-end on GDSC: cluster → MGMT recovery → one novel responder subgroup → driver annotation → static figures + a templated brief. **No LLM required.** | A complete, reproducible demo with a real result. |
| **B (target)** | Claude agent drives the protocol via tool-use: chooses features/params, judges signal-vs-noise, iterates, writes the brief in natural language. Add the **TCGA-LAML patient vignette**. | The "agent-mediated discovery" wow + real-patient anchor. |
| **C (stretch, only if Fri has room)** | Annotation grounded by querying the existing PolymerGenomicsAPI / methylation KB; light NVIDIA scale framing (cuML clustering / "scales to every methylome"). | Mechanistic grounding + the platform/NVIDIA horizon. |

**Order is strict A → B → (C).** A guarantees a demo exists before any agent risk is taken on.

---

## 7. Tech stack

- **Python.** pandas, numpy, scikit-learn (clustering, feature selection), scipy/statsmodels (association tests), lifelines (survival/KM), umap-learn, matplotlib + plotly (figures).
- **Agent:** Claude API tool-use (prompt-cached), components 2–4 exposed as tools. (Use the `claude-api` skill for the harness.)
- **Output:** existing static `site/` (vanilla HTML/CSS/JS) for the polished brief; a notebook or terminal view shows the live agent run on stage.
- **No new heavy infra.** GPU/RAPIDS is framing only (rung C), not a build dependency.

---

## 8. Scope — IN / OUT

**IN:** GDSC engine; MGMT positive-control recovery; one novel responder subgroup + driver annotation; agent-driven protocol (B); TCGA-LAML survival vignette; opportunity brief; backup video.

**OUT (say no all week):**
- ❌ cfDNA methylation (can't recover reliably in 2 days).
- ❌ CYP450 / liver pharmacology (not assessable from blood; future only).
- ❌ Deconvolution-vs-CBC lineage check (clever, but a side-quest — at most one validation slide).
- ❌ Rebuilding the Polymer claims universe / Evo 2 / v1.3 schema (archived as horizon, §9).
- ❌ Wet-lab, multi-indication breadth, controlled-access data wrangling.
- ❌ Pure on-thesis MPN cohort (GSE277866) as the engine — no response variable to cluster against; usable only as an optional visual cameo if everything else is done.

---

## 9. Disposition of existing material — keep everything, trim later

Nothing currently in the repo is moved or deleted. We build the `strata/` engine and its outputs **alongside** what exists:
- All existing docs (`VISION.md`, `DEMO_DAY_PLAN.md`, `PITCH_FRAME_DISCUSSION.md`, `NVIDIA_INTEGRATION_PLAN.md`, `METHYLOME_INTERPRETER*.md`, etc.) stay in place. They supply the closing "where this goes" vision (federated claims universe, Evo 2 coordinates, GPU-native verification) — referenced in the pitch, not rebuilt.
- The existing `site/` (MethylGKB front end) stays and is **extended** to surface the Strata demo output — not replaced.
- Any trimming/reorganizing happens *after* the demo works, never as a precondition for building.

---

## 10. Risks & mitigations

| Risk | Mitigation |
|---|---|
| MGMT doesn't surface cleanly in GDSC | Day-1 gate; fall back to another known methylation→drug association. Positive-control concept survives. |
| Live agent fails on stage | Deterministic floor (A) always works; record backup video Friday at highest rung. |
| "Just ML + a GPT wrapper" critique | The agent *iterates and judges signal-vs-noise* (B), and recovers a known result then extends it — that's discovery, not a wrapper. Lead with MGMT recovery as proof. |
| "Cell lines, not patients" critique | Owned explicitly: TCGA-LAML real-patient vignette + the locked-BeatAML honesty slide reframe it as *the structural reason these subgroups stay hidden.* |
| Solo time crunch | Strict A→B→(C); freeze features Friday; the brief can be templated (A) before it's agent-written (B). |

---

## 11. Why this wins the room

- **Creative:** an agent that rediscovers a clinical paradigm unsupervised, then finds a new one.
- **Impactful:** answers Pfizer's exact gap — subgroups the taxonomy can't see — with a real, recoverable result.
- **Honest:** the data limitation *is* the argument; nothing overclaimed.
- **Scales:** the method is indication-agnostic and GPU-shaped at scale — the bridge to the platform horizon and the NVIDIA lunch.
