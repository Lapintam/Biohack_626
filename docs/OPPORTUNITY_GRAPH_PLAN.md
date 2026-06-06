# Strata v2 — Cross-Indication Opportunity Engine · Plan for Review

> **Status (2026-06-06): Research pass COMPLETE (R1/R1b/R2/R3a). Saturday slice BUILT** — as-built engine/method/results in **`docs/V2_BUILD_SPEC.md`**; demo at **`site/strata-v2.html`** (additive; v1 untouched). This doc is now the *vision/strategy* record (how we got here + what's roadmap); the concrete build = `V2_BUILD_SPEC.md`. Current-state detail: `docs/SYSTEM_MAP.md`.
>
> Working name for the expanded system: **the Opportunity Graph** (placeholder).
>
> **Terminology guardrail:** "graph" here names the *target architecture*, **not a built substrate** — there is no schema, no edge/confidence model, no traversal engine yet. Until those exist, describe the deliverable as an **opportunity map / ranked opportunity table / provenance-ready claim layer / graph-shaped roadmap**. Do not imply a live knowledge graph in any demo or pitch until it is real. (External-audit recommendation, 2026-06-05.)

---

## 1. Why we're doing this (the honest trigger)

What we shipped — *"predict drug response from methylation"* — is a clean, honest **proof that our method recovers known biology**. But as a *thesis* it is **not enough**, for two blunt reasons:

1. **Pharma already does this.** GDSC was *built* for pharmacogenomic biomarker discovery; companion-diagnostic programs exist. Recovering known biology is, by definition, not a taxonomy-breaking discovery.
2. **It lives inside one box.** One drug, one tissue, one modality, within a single indication. The Pfizer challenge is explicitly about the opportunities **outside** the existing boxes.

This plan keeps the working engine as an *atomic operation* and builds the thing the challenge actually asks for around it.

## 2. The challenge, decoded into requirements

> *"…the patient subgroup hidden inside a heterogeneous indication, the mechanism that connects two conditions thought to be unrelated… the opportunities that would most change patient care are, by definition, the ones the existing taxonomy is least equipped to find… an AI-native approach [that] reveals opportunities the expert-gated model is structurally incapable of finding."*

| # | Requirement (their words → our reading) | Current demo |
|---|---|---|
| R1 | Evaluate **many** opportunities, cheap & fast (vs a handful over months) | ✗ one drug at a time |
| R2 | The **subgroup hidden in a heterogeneous indication** | ◑ within one indication only |
| R3 | The **mechanism connecting two *unrelated* conditions** (cross-indication) | ✗ never leaves one tissue |
| R4 | **Beat the taxonomy** — find what experts "structurally can't" | ✗ works inside one box |

**Load-bearing phrase:** *"the existing taxonomy is least equipped to find."* Expert review is organized by **indication** and **drug class**; the value is in what **cuts across** them.

## 3. The reframe (the whole idea in one move)

**Stop organizing by disease/drug taxonomy. Organize by molecular state.** Once you do, all three opportunity types fall out of one structure:

- **Hidden subgroup** = a molecular state that splits an indication *(our current atom)*.
- **Cross-indication connection** = the *same* state appearing in a *different* disease → a drug for disease A becomes a candidate for the marker-defined subgroup of disease B *(R3, taxonomy-breaking by construction)*.
- **Portfolio scan** = do this across all drugs × states × indications → a ranked **opportunity map** *(R1)*.

The methylation→response engine is **not the product** — it is the atom. **The product is the search over a molecular-state graph that crosses the taxonomy.**

## 4. Architecture — the composed system

Four layers (the four directions, unified):

| Layer | What it is | Direction |
|---|---|---|
| **1 · Substrate** | every sample = a **multi-omic molecular state** (methylation + expression + mutation + CNV), *not* a disease label | multi-omic state |
| **2 · Atom** | the existing marker→outcome engine relates state to response | *(reused)* |
| **3 · Graph** | nodes: states · markers · pathways · drug-targets · drugs · indications. edges: associations + mechanistic links. **Each edge a provenance-stamped claim.** | mechanism / knowledge-graph |
| **4 · Search** | an **exploratory agent** traverses the graph for high-value *non-obvious paths*, ranked into a **portfolio opportunity map** | cross-indication + exploratory agent |

**This is the parked Polymer claims-universe, now with a concrete commercial job:** a universe of pharmaco-molecular claims whose cross-taxonomy *paths* are the opportunities.

**Coverage check:** R1 (graph search ranks many) · R2 (state splits indication) · R3 (shared state + drug-target path across indications) · R4 (organized by molecular state, agent explores rather than recites). All four.

## 5. The rigorous spine (why this isn't vapor)

The within-indication thesis earned trust from a **positive control** (CDKN2A→palbociclib recovered cold). The cross-indication thesis has a *stronger* one, because **the FDA already validated "molecular state beats tissue taxonomy":**

- **MSI-high → pembrolizumab** (2017 — the first tissue-agnostic approval)
- **NTRK fusion → larotrectinib** · **BRAF → melanoma/colorectal/thyroid** · **HRD/BRCA → PARP inhibitors** (breast/ovarian/prostate/pancreatic)
- **MTAP/9p21 deletion → PRMT5/MAT2A inhibitors** (pan-cancer)

**The move:** the engine recovers one of these tissue-agnostic relationships **cold across GDSC's ~30 tissue types** (positive control), then proposes a **novel** cross-indication opportunity. Exact analog of CDKN2A, lifted to the cross-indication level.

**The thread we already hold:** **MTAP — our current #1 marker — is itself one of these pan-cancer targets.** The 9p21/MTAP/CDKN2A state crosses tumor types — and we have now **shown it cut across indications** (§5a).

## 5a. Empirical validation — R1 / R1b (2026-06-05) [CRUX RESOLVED]

The spine is no longer `[OPEN]` — it was tested. (`scripts/v2_tissue_confound_probe.py`, `scripts/v2_tissue_agnostic_scan.py`; full result in `notes/V2_TISSUE_CONFOUND.md`.)

- ✅ **The positive control is genuinely tissue-agnostic.** MTAP/CDKN2A → palbociclib survives **within-tissue** (significant in 5 unrelated tissues — lung, breast, colorectal, head & neck, ovary) **and** tissue adjustment (60–74% of the effect retained, q≈1e−6; still top-1% in a tissue-adjusted scan). The §8a safeguards pass. The cross-indication thesis is **real and confound-controlled** — evidence level **L3** for known biology.
- ✅ **The safeguards have teeth.** DAPK1 → luminespib (a v1 "novel hit") is mostly between-tissue (only 34% retained) → correctly demoted to PARTIAL/lineage.
- ⚠️ **Naive discovery is REJECTED.** Ranking all genes by tissue-adjusted association is **non-specific** (~8k of 14.6k co-significant; obscure top markers) — methylation's pervasive inter-gene co-correlation + a dominant within-tissue axis. The discovery engine **must be mechanism-anchored** (test a drug's known target/pathway genes) or global-axis-controlled — not brute-force ranking. Updates §7; converges with the external audit.

## 6. Staging — ambition *and* a Saturday demo

| Tier | Scope | Timeline |
|---|---|---|
| **Platform (the pitch / vision)** | full molecular-state graph + multi-omic + exploratory agent over patient + cell-line data | months |
| **Saturday-demoable slice** | one **worked cross-indication example**: a molecular state recovered as a *known* tissue-agnostic relationship (positive control) → a *novel* cross-indication opportunity → with a mechanism path (state→pathway→target→drug) → on a ranked **portfolio map** → narrated by the agent *exploring* | this week |
| **Roadmap (post-Saturday)** | full multi-omic state, full graph, patient cohorts, the claims-universe layer | later |

## 7. The Saturday slice, concretely

1. **Cross-tissue scan** — extend the engine: for a drug, test its marker state **across GDSC tissue types**, not within one. Surface states whose effect is tissue-agnostic. **Apply the §8a tissue-confound safeguards — pooled cross-tissue association ≠ tissue-agnostic biomarker.**
2. **Positive control** — recover a *known* tissue-agnostic relationship cold (MTAP/9p21, HRD, or MSI — whichever is cleanest, §9).
3. **Novel cross-indication hit — MECHANISM-ANCHORED (per R1b, not brute-force).** Do *not* take the top tissue-adjusted marker — it's non-specific. Instead test whether a drug's **known target / pathway genes** (from GDSC's own `PUTATIVE_TARGET` / `PATHWAY_NAME`, or ChEMBL / Reactome — pending R2) carry a *tissue-agnostic* methylation signal → a mechanistically-connected repurposing hypothesis. Label it by evidence level (§8a).
4. **Mechanism path** — annotate state→pathway→drug-target→drug so the connection has a *why*, not just a correlation. **[OPEN: which annotation sources — §9.]**
5. **Portfolio map** — the parked explorer, reframed as a **ranked opportunity map** across drugs × indications (this reuses `feat/strata-explorer`).
6. **Exploratory agent** — the agent *follows a lead* (find marker → chase across indications → connect mechanism → propose opportunity) instead of running the fixed list→rank→evaluate recipe.

## 8. Epistemic rails (carried forward + extended)

All current rails (`SYSTEM_MAP.md` Part 4) still hold. New/sharpened for v2:

- **Cross-indication signals from cell lines are HYPOTHESES**; the tissue-agnostic positive control licenses them. Never presented as validated.
- **Tissue-compartment discipline (`SYSTEM_MAP.md` Part 7) is now central** — a cross-indication signal must be attributed (tumor-intrinsic vs host vs PK) before it means anything clinically. *Cross-indication ≠ cross-compartment.*
- **Differentiate from prior art (resolved, §9.1).** Claim only **combination/framing novelty** — methylation-state ranked *across* the taxonomy + mechanism-anchored + provenance-stamped — and cite **Comms Biol 2023** + **CellMinerCDB**. The positive control (CDKN2A→palbociclib) is **confirmation, not discovery** (PLOS One 2019). Add **MGMT→TMZ as the honest negative control**; **MTAP→PRMT5i as the cross-indication narrative** (CNV-led caveat).
- **The graph is hypothesis-generating, not ground truth.** Edges are claims with provenance + confidence, not facts.

## 8a. Evidence discipline (folded in from the external audit, 2026-06-05)

Labeling/method discipline, **not new build** — these sharpen §8 and operationalize the plan's central risk (tissue confounding).

### Four evidentiary levels — never conflate them
A cross-indication claim is only as strong as the highest level it *honestly* reaches; label every opportunity by it:
1. **Cross-tissue marker association** — the state↔response correlation holds with tissues pooled. *(Weakest; may just be lineage.)*
2. **Tissue-agnostic biomarker** — the association survives **tissue adjustment / within-tissue** analysis. *(Real taxonomy-crossing signal.)*
3. **Cross-indication commercial opportunity** — a tissue-agnostic state maps to an approved/dev-stage drug in a *different* indication's subgroup. *(The actual product claim.)*
4. **Mechanism-supported repurposing hypothesis** — a provenance-stamped state→pathway→target→drug path supports it. *(Annotation, not proof.)*

### Evidence ladder — grade every opportunity
`L0` association in pooled GDSC · `L1` interpretable marker w/ known biology · `L2` not fully explained by tissue/lineage · `L3` recurs across multiple cancer types, consistent direction · `L4` mechanism path to target/drug (annotation) · `L5` external/orthogonal support (PRISM, CTRP, CCLE, literature) · `L6` clinical evidence.
**Where we are:** v1 reaches ~**L1**. The Saturday slice should target **L2–L3**, with **L4 as annotation, not proof.** Never present an L0/L1 result as if it were L3+.

### Tissue-confound safeguards (mandatory for any cross-indication claim)
The unsupervised null already proved blind structure = lineage. So any cross-tissue claim MUST:
- report the **tissue distribution** of marker-high vs marker-low groups + **per-tissue n**;
- compute **within-tissue effects** where n allows;
- fit a **tissue-adjusted** model (`AUC ~ marker + tissue` vs `AUC ~ marker`) and require the marker to survive;
- **penalize hits dominated by a single tissue/lineage**;
- call a state "recurrent across indications" **only when per-tissue sample counts support it**.

A hit that fails these is **lineage, not discovery** — and must be labeled as such.

## 9. Open questions — what the research pass must resolve [ALL OPEN]

1. **Prior art / novelty.** ✅ **RESOLVED (R2, `notes/V2_PRIOR_ART_ANNOTATION.md`):** closest prior art = **Comms Biol 2023** (GDSC methylation→drug, mechanism-anchored, but *cancer-type-stratified*, one-shot) and **CellMinerCDB** (methylation×drug, but *per-query*). NOT novel: methylation-as-feature, the tissue-agnostic paradigm, computational repurposing, even the CDKN2A→palbociclib link (published PLOS One 2019 → our control is **confirmation, not discovery**). **Defensible novelty = the *combination*:** methylation-state ranked *across* the taxonomy + mechanism-anchored (defeats the co-methylation confound, itself published in PMC7906305) + provenance-stamped, re-runnable. Combination/framing novelty only.
2. **Positive control + data feasibility.** ✅ **PARTLY RESOLVED (R1):** MTAP/9p21-CDKN2A → palbociclib is validated as tissue-agnostic in GDSC; tissue labels (30 tissues) are present and sufficient. *Still open:* additional controls (e.g. MLH1-methylation/MSI), and whether mutation/CNV labels are joinable for the multi-omic substrate. **New from R1b:** how to do *mechanism-anchored* discovery (target/pathway genes) or control the global within-tissue axis — depends on R2/§9.3 annotations.
3. **Annotation sources.** ✅ **RESOLVED (R2):** **GDSC2's own `PUTATIVE_TARGET` (87% of drugs) + `PATHWAY_NAME` (100%) columns suffice for Saturday — zero external dependency, no license risk.** Optional license-clean enrichment: DGIdb (MIT) for the 37 gaps, MSigDB Hallmark / Reactome (CC0/CC-BY) for pathway expansion. Avoid full DrugBank (CC BY-NC) + raw KEGG (paid) for any commercial build.
4. **Multi-omic substrate.** Do GDSC/CCLE expression + mutation + CNV join to our methylation lines, and is integration worth it for the Saturday slice or strictly roadmap?
5. **Demoable scope.** Given 1–4, what is the *minimum* worked example that still tells the full story?

## 10. Risks

| Risk | Mitigation |
|---|---|
| **Over-reach** — promise a platform, demo a toy | Strict staging (§6); the Saturday slice is a *complete* worked example, not a stub. |
| **Prior art overlap** — "this is just CMap" | Research pass (§9.1); pre-load the precise differentiation; lead with agent-native + mechanism attribution. |
| **Data feasibility** — annotations/labels not cleanly available | Research pass (§9.2–3) *before* building; fall back to the cleanest positive control. |
| **Honesty erosion** — cross-indication tempts overclaiming | Rails (§8); hypotheses + positive control; tissue-compartment attribution stays front and center. |
| **Breaking the working demo** — v2 destabilizes the shipped v1 | v2 is additive on a branch; the live site stays the bulletproof v1 until v2 is demo-proven. |

## 11. Explicitly NOT in scope for Saturday
- Full multi-omic integration, the full knowledge graph, patient-cohort data, the federated claims-universe, live arbitrary input. These are the *platform vision* — pitched as the horizon, not built.

## 12. Proposed next step
Run the **research pass (§9)** — prior-art differentiation + positive-control/data feasibility + a concrete demoable spec — *then* a design/implementation plan, *then* build. **This doc is the gate before any of that.**
