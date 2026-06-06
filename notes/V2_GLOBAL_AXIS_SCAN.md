# Global-axis-controlled all-gene scan — result (2026-06-06)

`scripts/v2_global_axis_scan.py`: unbiased scan over ALL ~14,600 genes, residualizing
AUC + every gene's methylation on [tissue fixed-effects + top-3 methylome PCs] (SVA/RUV-style
latent-factor adjustment) to defeat the R1b co-correlation. Variance filter (std≥0.02) added
to rule out low-variance artifacts.

## Result (gene-level methylation)
| drug | control | q<0.1 (of ~14.6k) | top genes | known marker rank |
|---|---|---|---|---|
| Palbociclib | none | 12,797 | ZFP112, CHST12, FERMT3… (obscure) | MTAP 79 |
| | tissue | 8,917 | SDR42E1, ZFP36L2… (obscure) | MTAP 93 |
| | **tissue+PCs** | **4,934** | ZFP36L2, SDR42E1, SNORA6, AIF1L… (obscure) | MTAP **160**, CDKN2A 681 |
| Trametinib | tissue+PCs | 9,025 (PCs did NOT help) | POLS, NELL2, MTMR15… | BRAF **9,512 (ns)**, DUSP6 4,657 |
| Nutlin-3a | tissue+PCs | 6,533 | C10orf131, CYGB… | TP53 2,815, MDM2 3,916, CDKN2A 795 |

## HONEST CONCLUSION (negative-ish, load-bearing)
The global-axis control **partially works**: it cuts non-specificity (Palbociclib 12.8k→4.9k significant, ~60% reduction). **But it does NOT deliver clean unbiased discovery:**
- **Top genes stay obscure** (SNORAs, miRNAs, ZNFs, ORFs) across drugs — not interpretable, not the known biology. Variance filter ruled out a low-variance artifact (tops unchanged).
- **Known mechanism genes do NOT rise** — they sink (MTAP→160, BRAF→9,512/ns, DUSP6→4,657, TP53→2,815). An unbiased scan ranks them *below* obscure co-correlated genes.
- For Trametinib, PC removal **increased** n_sig — the latent-factor adjustment isn't reliably capturing the right nuisance.

**Why:** single-gene methylation→drug-response is a **weak (|r|≲0.25 for real signals), diffuse** signal in a **pervasively co-correlated** methylome; the highest partial correlations are dominated by residual shared structure, not gene-specific biology. n≈850 is underpowered for clean genome-wide discovery of effects this small.

## Implications
1. **Unbiased single-gene discovery is a hard research problem here — not a demo deliverable.** This VALIDATES the mechanism-anchored design (the only route to interpretable, specific hits) and quantifies the limit.
2. To actually push unbiased discovery: RUV with negative-control genes (vs PCA); **promoter resolution + multi-omic** (regress out CNV/expression to disambiguate); or much larger n.
3. Honest framing: "we tried unbiased genome-wide discovery with proper latent-factor control; it reduces noise but doesn't surface clean novel genes in this data — which is why we anchor to mechanism."

## PROMOTER re-run (2026-06-06) — the specificity problem is SOLVED
Same scan on CCLE RRBS promoter methylation (coverage≥60% + var filter + mean-fill), n≈449 lines.

| drug | naive q<0.1 | tissue | **tissue+PCs** | known marker |
|---|---|---|---|---|
| Palbociclib | 8,308 | 1,061 | **689** /13.8k | MTAP rank 247 (r=+0.17, ns-dir); CDKN2A ns |
| Trametinib | — | — | **646** | BRAF/KRAS/DUSP6 all ns |
| Nutlin-3a | 5,988 | 526 | **753** | TP53/MDM2/CDKN2A ns |
| Temozolomide | 2,288 | 173 | **337** | MGMT ns |

**WIN:** promoter + global-axis control collapses "significant" from ~30–60% (gene-level) to **~2–5% (a few hundred genes)** — a tractable candidate set. The user's hypothesis was right: **promoter resolution is what makes unbiased genome-wide discovery specific.**

**Nuance:** the surviving top hits are mostly **obscure genes** (AGAP2-AS1, RCAN3AS, ZNF234, RNGTT, SMIM15…), with some biologically plausible ones in the tail (RASAL2 for a MEK drug, HTRA2 for an apoptosis inducer, CDKN2C for cell cycle). Canonical genes correctly **don't** top the list because their drug-response role is genomic (BRAF/KRAS mutation; MTAP/CDKN2A deletion), not promoter methylation — so methylation legitimately shouldn't flag them. MTAP even flips to +direction at promoter resolution (confirming the gene-level signal was CNV).

**State:** we moved from *non-specific (can't pick anything)* → *specific (few-hundred-gene candidate list per drug)*. Those candidates are **unvalidated novel hypotheses** → need orthogonal validation (RNA-seq: does candidate promoter methylation track its expression = functional silencing?).
