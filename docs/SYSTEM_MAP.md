# MethylGKB · Strata — System Map: Data Assets & Inference Rules

> A precise, source-derived inventory of **everything we have**: every data asset (what it is, where it came from, its exact shape), every inference rule (the methods that turn data into claims, with exact parameters), how they compose, and the epistemic rails that keep it honest. Written 2026-06-05 from the actual code + provenance notes, not memory.

---

## PART 1 — DATA ASSETS

All raw data is **external + public**. Everything else in the repo is *derived* from it (Part 5). Nothing here is patient-identifiable; GDSC is lab cell lines, TCGA-LAML is a public research cohort.

### 1A. GDSC — drug-response engine substrate (Cell Model Passports / Sanger)

> `cancerrxgene.org` was decommissioned; GDSC now lives at Cell Model Passports. See `notes/GDSC_SCHEMA.md`.

| Asset | File (gitignored, in `data/gdsc/`) | Source | Shape / contents | Loader |
|---|---|---|---|---|
| **Methylation** | `methylation_imputed.csv.gz` (59 MB gz, ~175 MB raw) | figshare file 43146652 (CMP multi-omics) | **Gene-level imputed β** — source 1,523 lines × **14,608 gene-symbol columns**; after join → **~1,055 lines × 14,608 genes**, β ∈ [0,1] | `load_gdsc_methylation()` |
| **Drug response** | `GDSC2_fitted_dose_response_27Oct23.xlsx` (21 MB) | `cmp.cog.sanger.ac.uk` | Long: **242,036 rows, 286 drugs, 969 lines**. Cols → `COSMIC_ID, drug_name, ln_ic50, auc`. **Lower AUC = more sensitive.** | `load_gdsc_drug_response()` |
| **Annotations** | `model_list_20260420.csv` (936 KB) | `cog.sanger.ac.uk/cmp` | COSMIC_ID-indexed: `cell_line_name, tissue, cancer_type, sanger_model_id` | `load_gdsc_annotations()` |

- **Join key:** the native key across all three is `SANGER_MODEL_ID` (SIDM). The drug file has **no COSMIC_ID** — loaders map SIDM→COSMIC_ID via `model_list` and drop rows lacking a COSMIC_ID (methylation 1,055/1,523 retained; all 969 drug lines map).
- **Critical caveat:** methylation is **gene-level imputed** (one β per gene, gene-body aggregated), **NOT CpG-island / probe-level**. This is *the* reason the textbook MGMT→temozolomide control is null here (Part 4).

### 1B. TCGA-LAML — real-patient survival axis (UCSC Xena, GDC hub)

> See `notes/TCGA_SCHEMA.md`.

| Asset | File (in `data/tcga_laml/`) | Source | Shape / contents | Loader |
|---|---|---|---|---|
| **Methylation** | `TCGA-LAML.methylation450.tsv.gz` (633 MB) | UCSC Xena GDC-LAML hub | Illumina **450K** β, probes×samples → transposed to **samples×probes**; after QC: **169 samples × ~417,326 probes**, β ∈ [0,1] | `load_tcga_methylation()` |
| **Survival** | `TCGA-LAML.survival.tsv.gz` (4 KB) | same hub | `os_time` (days), `os_event` (1=dead, 0=censored); ~190+ patients, **~169 join** with methylation | `load_tcga_clinical()` |

- **Join:** TCGA aliquot barcode (e.g. `TCGA-AB-2820-03A`), identical in both files; first aliquot per barcode kept.
- **Missing-value rule (methylation):** drop probes missing in **>20%** of samples, then fill remaining NaN with the **per-probe mean** (done in probes-as-rows orientation for speed). Result is dense.
- This is real-patient, blood/marrow myeloid — but the only outcome is **survival**, not drug response.

### 1C. Not available (named so we don't chase it)
- **BeatAML methylation + ex-vivo drug AUC** — the ideal substrate (blood methylation + real drug response) — is **dbGaP `phs001657`, controlled-access**. Not obtainable. This locked dataset is itself the pitch's "why these subgroups stay hidden."

### 1D. ⚠ Methylation RESOLUTION — gene-level now, probe/promoter-level as the upgrade

> **Under review 2026-06-05 — possible pre-Saturday change.** This is the most important data caveat in the system.

**What we have (and who did it):** the GDSC methylation we use is **provider-imputed, gene-level** — one β per gene (~14,608 genes), distributed that way by Cell Model Passports (filename `…imputed_methylation.csv`, columns are bare gene symbols). The underlying assay is probe-level (Illumina array, ~450k CpG sites); **CMP collapsed probes→genes and imputed missing values. We did NOT aggregate or impute it** — we consume it as-is. (The *only* imputation we perform anywhere is a light per-probe mean-fill on the probe-level TCGA matrix, §1B.)

**Why it matters:** the methylation that *regulates* a gene is concentrated at its **promoter CpG island**, and promoter vs gene-body methylation can differ — even oppose. Averaging the whole gene into one number **washes out the promoter-specific signal.** This is precisely why **MGMT→temozolomide is null here** (rail #4): MGMT's effect is a promoter event, diluted at gene resolution. CDKN2A/MTAP survived because their signal is strong enough to persist through the averaging. Gene-level = simpler (feature = gene) but **lower resolution.**

**Upgrade path → promoter/probe resolution (candidate sources, verify current files before relying):**
- **CCLE RRBS promoter methylation** (DepMap data page; Ghandi 2019) — ~843 cell lines × ~17,182 **gene-promoter regions**. This is *promoter-level* — the natural fix. Join to GDSC drug AUC (or DepMap **PRISM** drug response) by cell-line ID (name / Broad / Sanger crosswalk).
- **GEO Illumina 450K cell-line methylation series** — full probe/CpG-island resolution if a clean cell-line β-matrix is available.
- **Method change is minimal:** `rank_response_markers` is feature-agnostic — point it at a probe/promoter β-matrix instead of the gene matrix and the same Spearman+FDR scan runs at the higher resolution. Cost: more features (~17k–450k), more memory/wrangling, fewer lines (~843 for CCLE), and a cell-line ID join.
- **Payoff:** likely recovers MGMT→TMZ and sharpens every marker (promoter-specific, not gene-averaged). The honest, resolution-correct version of the demo.

### Constants (`strata/config.py`)
`N_TOP_FEATURES = 2000` · `RANDOM_STATE = 0` · `DATA_DIR`, `OUTPUT_DIR`.

---

## PART 2 — INFERENCE RULES

Each rule is a pure function in `strata/engine/` (plus the agent in `strata/agent/`). Exact parameters below.

### 2A. `select_variable_features(meth, n_top=2000)` — `engine/features.py`
Rank columns by variance, return the top `n_top` gene/probe names. **Used only to build the UMAP landscape**, not for discovery.

### 2B. `cluster_samples(meth, features, k=None, random_state=0)` — `engine/cluster.py`
KMeans (`n_init=10`); when `k is None`, pick k∈[2,6] by **best silhouette**. Returns a per-sample cluster label Series.
**Status: NOT the discovery mechanism.** Empirically it tracks **lineage** (blood vs solid), so it is used only for the visual landscape. (Evidence: under clustering, CDKN2A ranks 9858/14608; under the marker scan it ranks 120.)

### 2C. `rank_response_markers(meth, drug_auc, min_n=30)` — `engine/associate.py` — **THE CORE DISCOVERY RULE**
For each of the ~14,608 genes: `spearmanr(gene β, drug AUC)` over shared non-NaN samples (n ≈ 960). Then **Benjamini-Hochberg FDR** (`multipletests … fdr_bh`); keep **q < 0.1**; sort by **|rho| descending**. Returns `[gene, rho, pvalue, qvalue, n]`.
- **Sign convention:** lower AUC = more sensitive, so **negative rho = methylation marks sensitivity**; positive rho = marks resistance.
- This is a **hypothesis-generating screen** (it selects on the outcome — see Part 4).

### 2D. `discover_responder_subgroup(meth, drug_auc, marker)` — `engine/associate.py`
**Median-split** lines on the chosen marker's β (`1` = high methylation, `0` = low), then call `test_drug_association`. This turns a continuous marker into an operational responder/non-responder **subgroup** (the "companion-marker cut").

### 2E. `test_drug_association(clusters, drug_auc)` — `engine/associate.py`
Drop NaN AUC; if <2 groups → neutral result. **Kruskal-Wallis** across the clusters' AUC distributions. `responder = cluster with min mean AUC`. `effect_size = median(AUC) − responder mean AUC`. Grade via `score_responder`. Returns `AssociationResult(effect_size, pvalue, responder_label, grade)`.

### 2F. `score_responder(...)` — grade thresholds
`delta = baseline(=median AUC) − cluster mean AUC`: `++` if >0.20 · `+` if >0.07 · `neutral` in [−0.07, 0.07] · `−` if <−0.07 · `−−` if <−0.20.
(Note: real median-split effects are ~0.03–0.05, so the headline grade is usually **`neutral`** even when p is tiny — by design; see Part 4.)

### 2G. `differential_cpgs(meth, clusters, target_label)` + `annotate_cpgs(genes)` — `engine/annotate.py`
Per-gene **Mann-Whitney** subgroup-vs-rest, sorted by |Δβ|. `annotate_cpgs` is identity here (the feature **is** the gene symbol). Used for driver lists/heatmaps.

### 2H. `test_survival_association(clusters, time, event)` — `engine/associate.py`
**Multivariate log-rank** across clusters + per-cluster **Kaplan-Meier** curves. Returns `{logrank_p, km}`. Currently fed by *unsupervised* clusters → null (Part 4). A marker-driven survival scan is **not yet built** (and would need the leakage guard, Part 4).

### 2I. The agent — `strata/agent/loop.py` + `tools.py`
**Four tools** (`tools.py`, JSON-safe, never raise): `list_drugs`, `rank_markers` (wraps 2C, returns top-N + `n_significant` + `n`), `evaluate_subgroup` (wraps 2D, returns grade/p/effect/responder/n), `write_opportunity_brief` (echo).
**Protocol** (`SYSTEM_PROMPT`): `list_drugs → rank_markers → pick the most biologically interpretable top marker (prefer known biology, e.g. CDKN2A/MTAP for CDK4/6i) → evaluate_subgroup → write brief`.
**Two paths, identical output shape** `{drug, goal, transcript[], brief, mode}`:
- `_run_llm` — Anthropic tool-use loop, model **`claude-opus-4-8`**, `MAX_TURNS=10`, `MAX_TOKENS=2000`, cached system block. Used when `ANTHROPIC_API_KEY` set.
- `_run_scripted` — deterministic emulation (calls `dispatch_tool` directly). Used with no key.
- **Safety:** any LLM error (auth/credit/network) is caught → falls back to scripted, records `llm_error`. Never crashes.
**Positive-control rule:** `_KNOWN_CONTROL_GENES = {CDKN2A, MTAP}`; when the top marker is in this set, the "recovers known biology" framing fires.

### 2J. Scale projection — `scripts/scale_benchmark.py` → `site/demo/scale.json`
Measures real per-correlation CPU cost (the 2C kernel) and **linearly extrapolates** to biobank scale (20k genes × 1k outcomes × 1M cohort). Today: **4.18M correlations / 746 s CPU**; projected: **~43 CPU-days** for one pass (≈5,000×). Labeled a projection, not a GPU run.

---

## PART 3 — HOW IT COMPOSES (the flow)

```
DATA                      DISCOVERY (per drug)                         OUTPUT
────                      ────────────────────                         ──────
GDSC methylation  ┐
(1055×14608 genes)├─► rank_response_markers (2C)  ──► ranked markers ──┐
GDSC drug AUC ────┘     Spearman per gene + FDR        (rho, q)        │
   (286 drugs)                    │                                    │
                                  ▼ pick interpretable top marker      │
                       discover_responder_subgroup (2D)                ├─► agent brief
                          median-split → Kruskal (2E/2F)               │   (LLM or scripted, 2I)
                                  │                                     │
                                  ▼                                     │
                       figures (UMAP/AUC/heatmap) + driver genes (2G) ──┘
                                                                        └─► site/demo/{drug}.json

TCGA methylation + survival ─► cluster_samples (2B, unsupervised) ─► test_survival_association (2H)
   (169 patients)                                                      └─► KM + log-rank  (currently NULL)
```

**One-line essence:** for a given outcome column, correlate every gene's methylation against it, FDR-rank the markers, cut the responder subgroup at the top marker, and have an agent narrate + write it up — with a positive control recovering known biology as the trust anchor.

---

## PART 4 — EPISTEMIC RULES (the honesty rails — most important section)

These are *load-bearing*. They are encoded in the code/copy and must survive a sharp reviewer.

1. **Trust comes from recovering known biology, not from the subgroup p.** The screen recovers **CDKN2A/MTAP → Palbociclib** cold (positive control). That licenses the novel hits.
2. **The screen selects on the outcome → the median-split subgroup p is CIRCULAR.** It is *not* independent validation. **Always lead with the continuous Spearman rho / q-value**; present the split as the operational cut only. (Encoded in `SYSTEM_PROMPT` framing rules.)
3. **Effect sizes are small (~0.03–0.05 ΔAUC); grades read `neutral` despite tiny p** (p inflated by n≈960). Never headline the split effect.
4. **Methylation is gene-level imputed (provider-side), not CpG-island.** Consequence: **MGMT→temozolomide is NULL** here (rho=−0.018, p=0.58) — promoter-island silencing is averaged out. This is why the positive control is CDKN2A/MTAP, not MGMT. (`notes/MGMT_GATE.md`.) **→ resolution upgrade path documented in §1D (under review, possibly pre-Saturday).**
5. **Unsupervised clustering is a lineage detector, not the discovery axis** (blood vs solid). Discovery is the marker→response scan, not clustering.
6. **The TCGA survival split is an honest NULL** (log-rank p=0.84). Kept as a real-patient *landscape*, never claimed as a survival result.
7. **GDSC = lab cell lines, public data.** This is a **proof-of-method**, not a clinical product. Novel hits are **leads to test**.
8. **A future survival/any-outcome marker scan inherits the same leakage** (select-on-outcome). It must be framed exploratory, with the rho/q leading and (ideally) held-out or external validation before any claim.
9. **Don't conflate tissue compartments.** Tumor (PD), liver (PK), and blood (host/immune) methylomes are largely **tissue-specific and not interchangeable**; blood is a poor proxy for liver/tumor at regulatory loci. GDSC = **tumor-only** (clean, single compartment). In a clinical cohort a methylation subgroup is a **correlation, not a mechanism**, without metadata — see Part 7.

### Established results (from `notes/DEMO_HITS.md`, `NOVEL_HITS.md`)
| Role | Marker → Drug | rho | q | Note |
|---|---|---|---|---|
| **Positive control** | CDKN2A → Palbociclib | −0.336 | 7.8e−25 | rank 120/12,333; **top marker MTAP** rho=−0.440, q=6.3e−43 (9p21 co-deletion) |
| Novel hit A | DAPK1 → Luminespib | −0.346 | 3.2e−26 | HSP90i; top marker COBL (−0.438). Hypothesis. |
| Novel hit B | FKBP1A / RASSF1 → Vorinostat | +0.607 / +0.419 | 1.3e−93 / 1.4e−40 | HDACi; resistance-direction. Hypothesis. |

---

## PART 5 — DERIVED ARTIFACTS (built by us, from Part 1)

- **`site/demo/`** (committed, served on Vercel): `Palbociclib.json` (**live agent**, `mode:llm`), `Luminespib.json` (scripted), 3 figures per drug, `scale.json`, `index.json`. The public page replays these — **no API calls / no key in the browser**.
- **`outputs/`** (gitignored): per-drug briefs + figures from CLI runs.
- **`notes/`**: `GDSC_SCHEMA.md`, `TCGA_SCHEMA.md`, `MGMT_GATE.md`, `DEMO_HITS.md`, `NOVEL_HITS.md`, `NVIDIA_PITCH.md` — provenance + results + pitch.
- **`docs/`**: the design spec + plan, `strata-explainer.html` (→ PDF), this map.
- **Live site:** `https://methylkb.vercel.app` (static).

---

## PART 6 — WHAT GENERALIZES (and the open build)

- **The method is outcome-agnostic by construction:** `rank_response_markers` only needs `(methylation matrix, outcome vector keyed by sample)`. Swap drug AUC for survival, toxicity, progression → same loop. *That* is the "general engine → clinical decision support" thesis.
- **Currently wired:** drug response (286 drugs, fully). Survival exists as data + an unsupervised null; a **marker-driven survival scan is not yet built**.
- **In-flight (parked on branch `feat/strata-explorer`):** a static explorer over all 286 drugs + survival axis + an extensible outcome manifest, plus a Fly-deployable FastAPI backend for live on-demand scans (front-end gets a static↔live switch). Nothing committed yet.

---

## PART 7 — Tissue compartments & mechanistic attribution (the clinical-translation constraint)

DNA methylation is **strongly tissue-specific** — most variable CpGs differ between tissues, so one compartment's methylome is a **poor proxy** for another's. Blood↔liver concordance is **locus-dependent**: a minority of (often genetically-driven, shared-mQTL) CpGs track across tissues, while regulatory, environmentally-shaped CpGs do not. Pan-tissue epigenetic clocks work *because* they are built from a curated cross-tissue-stable CpG set — the exception that proves the rule. For drug metabolism specifically, **CYP450/UGT/transporter regulation is hepatic** and is **not** readable from peripheral blood.

**Three compartments, three distinct roles — do not conflate:**

| Compartment | Methylation governs | Role in the drug-response chain | In our data? |
|---|---|---|---|
| **Tumor** | cell-intrinsic sensitivity (MGMT, CDKN2A…) | **Pharmacodynamics** — does the target respond | ✅ GDSC cell lines model this |
| **Liver** | drug-metabolizing enzymes (CYP450, UGT), transporters | **Pharmacokinetics** — how much drug reaches the target | ❌ hepatic-only; not readable from blood |
| **Blood / immune** | hematopoiesis, immune state, CHIP, inflammation | **Host context** — immune-mediated efficacy + toxicity | ❌ (the PB thesis; a different assay) |

- **GDSC is the tumor compartment only** — a single clean compartment, no PK/host confound. That is *why* it is a clean proof-of-method, and why the honest scope is "tumor-intrinsic sensitivity in cell lines."
- **A real clinical cohort mixes all three** (a blood draw = immune cells; cfDNA = tumor + normal; a biopsy = tumor + stroma). A methylation-defined responder subgroup could be driven by any compartment — or by a confound (age, cell-type composition, genotype/mQTL, batch).

> ### Principle — Clinical Cohort Response Subgroups Require Metadata for Mechanistic Attribution
> A methylation-defined response subgroup is, on its own, a **black-box correlation**. Moving from *which* subgroup responds to *why* requires metadata that attributes the signal to a compartment and a mechanism: **tissue source · cell-type deconvolution · genotype/mQTL · age/clock · mutation/CNV · treatment · batch.** Without it you have a predictor, not a mechanism — and nothing safe to act on clinically. **The engine finds the subgroup; the metadata layer earns the attribution.**

This turns the honest limitation into the roadmap: it specifies *what to collect* in a real assay (paired compartments + covariates), and it aligns with the paired-design Mol-CBC architecture (paired sampling = the attribution layer).
