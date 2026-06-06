# Polymer — Model Architecture, Inputs, and Training Results

Polymer predicts **drug response from DNA methylation** and explains *why*: given a tumor/cell-line
methylome it ranks every drug by predicted sensitivity, and surfaces the hyper-/hypo-methylated genes
driving each prediction. This document describes the model, its inputs, how it is trained, and how it
performs.

> **Scope/honesty.** Trained on **GDSC cancer cell lines** (Cell Model Passports), not patients.
> Methylation is **gene-level** beta (promoter-CpG detail is averaged out, except where the promoter
> matrix is used explicitly). This is a methodology demonstration on open data, not a clinical device.

---

## 1. What the model does

For a methylome \(x\) (per-gene beta values) and a drug \(d\), Polymer predicts the dose-response
**AUC** (lower AUC = more sensitive):

\[
\widehat{\text{AUC}}(x, d) \;=\; b_d \;+\; \big\langle\, f_\theta(x),\; e_d \,\big\rangle
\]

- \(f_\theta(x)\) — a shared **methylation encoder** (MLP) → 32-dim cell embedding.
- \(e_d\) — a learned **drug embedding** (32-dim), \(b_d\) — a learned **per-drug bias**.
- \(b_d\) absorbs each drug's *average* potency (the same for everyone); the dot product
  \(\langle f_\theta(x), e_d\rangle\) is the **personalized** interaction — the only term that can say a
  drug is differentially good *for this methylome*.

The encoder is **pretrained across all 286 drugs at once** (multi-task), so a single representation
serves every drug and the model can rank the full panel for a new methylome.

---

## 2. Inputs

### Predictor input (features)
| Item | Detail |
|---|---|
| Methylation matrix | Gene-level beta, COSMIC_ID-keyed (`methylation_genelevel_bycosmic` / `methylation_imputed`), ~**14,608** gene-symbol columns, beta ≈ [0,1] |
| Feature selection | Top-**3,000** most-variable genes, selected on **train lines only** |
| Standardization | z-scored per gene using **train** mean/SD (no leakage) |
| Drug identity | 286 GDSC2 compounds, integer-encoded → learned embedding |

### Label
| Item | Detail |
|---|---|
| Target | GDSC2 fitted dose-response **AUC** (one value per cell-line × drug), standardized on train pairs |

### Auxiliary inputs (not used by the predictor)
| Item | Use |
|---|---|
| Promoter methylation (`methylation_promoter_bycosmic`, ~17,181 promoters, has NaNs) | Focal **hyper/hypo annotation** of genes of interest; per-gene association (recovers MGMT→temozolomide) |
| `model_list` annotations (tissue, cancer type) | Baselines, per-cancer analysis |
| Drug target / pathway (from `gdsc2_merged.csv`) | Display labels and interpretation |

### Data scale (default `imputed` source)
- **1,055** cell lines × **14,608** genes; **286** drugs; **242,036** (line, drug) AUC pairs.

---

## 3. Architecture

```
methylome x (3,000 z-scored genes)
        │
        ▼
  ┌──────────────────────── encoder f_θ ────────────────────────┐
  │ Linear(3000 → 256) → ReLU → BatchNorm1d(256) → Dropout(0.3)  │
  │ Linear(256 → 32)                                            │  → cell embedding (32-d)
  └─────────────────────────────────────────────────────────────┘
        │
        │   drug d ──► Embedding(286 → 32) = e_d ,  Embedding(286 → 1) = b_d
        ▼
  pred_AUC = b_d + dot(cell_embedding, e_d)
```

| Component | Shape / params |
|---|---|
| Encoder L1 | Linear 3000→256 (+BatchNorm, Dropout 0.3) |
| Encoder L2 | Linear 256→32 |
| Drug embedding | Embedding(286, 32) |
| Drug bias | Embedding(286, 1) |
| Cell / drug embedding dim | 32 |

This is **content-based matrix factorization**: a recommender where the "user" (cell line) is encoded
from methylation content and the "item" (drug) is a free embedding.

---

## 4. Training

| Setting | Value |
|---|---|
| Loss | MSE on standardized AUC |
| Optimizer | Adam, lr = 2e-3, weight_decay = 1e-5 |
| Batch size | 4,096 pairs |
| Epochs | 40 (early stopping on validation RMSE) |
| Split | **Cell-line-blind**: 70% train / 15% val / 15% test lines (seed 1337) — every test line is unseen |
| Leakage controls | Gene selection, z-scoring, AUC standardization, and per-drug means all fit on **train only** |
| Promoter NaNs | Imputed with **train-line** column means (leakage-safe) when promoter is used |

Reproduce:
```bash
python train_drug_recommender.py                 # gene-level (default)
python train_drug_recommender.py --source promoter
python polymer_benchmark.py leaderboard           # scored leaderboard vs baselines
jupyter notebook polymer_demo.ipynb               # training curves + interpretability
```

---

## 5. Results (held-out cell lines)

### 5.1 Headline benchmark (gene-level, 158 test lines, 34,762 pairs)

| Model | RMSE | R² | Pearson r | **personalized R²** | within-line ρ (de-meaned) | top-5 recovery |
|---|---|---|---|---|---|---|
| **Polymer** | **0.074** | **0.737** | **0.859** | **0.382** | **0.488** | 3.98 / 5 |
| PerDrug + Tissue | 0.083 | 0.671 | 0.819 | 0.228 | 0.357 | 3.89 |
| PerDrug mean | 0.094 | 0.574 | 0.757 | ≈0 (by design) | — | 3.77 |
| Global mean | 0.144 | 0.000 | — | −1.35 | 0.13 | 0.11 |

- **`personalized R²`** is the honest metric: R² computed on each drug's *mean-removed residuals*, so a
  per-drug-mean model scores 0. Polymer reaches **0.38** there — real personalization beyond "which
  drugs are potent on average."
- The raw R² (0.74) is dominated by global drug potency (a few drugs, e.g. Sepantronium, are sensitive
  on nearly every line); personalized R² strips that out.

### 5.2 Methylation predictability across all 286 drugs
- **280 / 286** drugs are more methylation-predictable than the tissue baseline.
- Mean Pearson r rises **0.32 → 0.48**.
- Most-predictable drug class: **PI3K/AKT/mTOR inhibitors** (Pictilisib, Taselisib, Buparlisib,
  Dactolisib, Alpelisib, …) plus EGFR (Afatinib) — a coherent biological signal.

### 5.3 Association validation (positive + recovered controls)
| Association | Resolution | Result |
|---|---|---|
| CDKN2A → Palbociclib | gene-level | ρ = −0.336, p = 6×10⁻²⁷ (textbook control, recovered) |
| MGMT → Temozolomide | gene-level | ρ = −0.018, p = 0.58 (null — averaged out) |
| MGMT → Temozolomide | **promoter** | ρ = −0.109, p = 0.011 (**recovered** at promoter resolution) |

### 5.4 Methylation source comparison (identical 557-line cohort + split)
| Feature set | features | personalized R² |
|---|---|---|
| **gene-level only** | 3,000 | **0.346** |
| concat (gene + promoter) | 6,000 | 0.271 |
| two-tower gated fusion | 3,000+3,000 | 0.220 |
| promoter only | 3,000 | 0.196 |

**Finding:** gene-level alone wins. Naive fusion and gated two-tower fusion both *hurt* — promoter, as a
bulk feature set on this smaller cohort, adds noise. Promoter's value is **focal** (specific silencing
markers like MGMT), not distributed; it is best used as a biomarker lens, not bulk model input. Polymer
therefore uses **gene-level** as its predictor.

---

## 6. Interpretability — genes of interest

For a given methylome and a chosen drug, Polymer attributes the prediction to genes via
**gradient × methylation deviation**:

\[
\text{contribution}_i \;=\; \frac{\partial\,\widehat{\text{AUC}}}{\partial z_i}\;\cdot\; z_i
\qquad (z_i = \text{the gene's standardized beta} = \text{its hyper/hypo score})
\]

- \(z_i > 0\) → **hyper-methylated** vs cohort; \(z_i < 0\) → **hypo-methylated**.
- **Negative contribution** → that gene's (hyper/hypo)methylation **drives sensitivity** (lowers AUC);
  positive → drives resistance.

This yields an auditable, gene-level rationale for each recommendation (see `polymer_demo.ipynb`,
Stages 4–5), combined with the genes' promoter-beta annotation.

---

## 7. Files

| File | Role |
|---|---|
| `train_drug_recommender.py` | Train the recommender (`--source imputed|genelevel|promoter`); saves `outputs/polymer_model*.pt` |
| `polymer_infer.py` | Load artifact, parse a pasted methylome, return ranked drugs |
| `polymer_app.py` | Streamlit paste-in demo UI |
| `polymer_benchmark.py` | Model-agnostic benchmark + leaderboard (`--source ...`) |
| `polymer_demo.ipynb` | End-to-end notebook: training curves, prediction, gene attribution (matplotlib) |
| `methylation_auc_association.py` | Per-gene methylation↔AUC association (`--source ...`) |
| `train_combined_methylation.py` | Gene vs promoter vs concat vs two-tower comparison |
| `strata/data/gdsc.py` | Data loaders (gene-level, promoter, drug response, annotations) |

---

## 8. Limitations

- **Cell lines, not patients** — clinical translation requires patient methylomes (e.g. controlled-access
  cohorts). The framing is right; the clinical claim is not yet earned.
- **Gene-level resolution** averages promoter-CpG silencing (mitigated by the promoter matrix for focal
  markers).
- **Not cold-start for drugs** — all 286 drugs are seen in training; predicting an entirely new compound
  would need a chemistry-based drug tower.
- **Attribution is first-order/local** — gradient×deviation is a local explanation, not a causal claim.
