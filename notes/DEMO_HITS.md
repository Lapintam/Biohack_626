# Demo hits — marker→response discovery (Rung-A, marker-driven)

Date: 2026-06-04. Mechanism: `rank_response_markers` Spearman-correlates each gene's
methylation against drug AUC (BH-FDR, q<0.1), then `discover_responder_subgroup`
median-splits lines on the chosen marker and grades via the existing Kruskal machinery.
**Lower AUC = more drug-sensitive**, so negative rho = methylation marks SENSITIVITY,
positive rho = methylation marks RESISTANCE.

**Honest caveat on effect size:** the strength of every hit lives in the *continuous*
Spearman rho (and its tiny q), not in the dichotomized median-split. After median-split
the AUC effect sizes are modest (~0.04–0.05) and the grade lands at `neutral`, even
though the split p-values are astronomically significant (driven by large n≈960). The
demo should lead with rho/q, and present the median-split as the operational
companion-marker cut, not as the headline effect.

---

## (1) Positive control — CDKN2A → Palbociclib  [RECOVERED]
- **CDKN2A**: rank **120 / 12333** significant markers, rho = **-0.336**, q = **7.8e-25**
  (matches the established gate value rho=-0.336). Median-split subgroup p=7.3e-18,
  effect 0.033, grade neutral.
- **Top marker for Palbociclib = MTAP** (rank 1, rho=-0.440, q=6.3e-43). MTAP sits
  immediately adjacent to CDKN2A on 9p21 and is co-deleted / co-methylated with it —
  so the #1 marker IS the CDKN2A locus, just the neighboring gene.
- **Verdict:** the marker mechanism recovers the positive control cleanly (CDKN2A in
  the top ~1% with the exact expected rho), unlike blind clustering where it ranked
  9858/14608. CDKN2A is NOT top-20 only because ~119 genes (led by its own 9p21
  neighbor MTAP) correlate marginally more strongly.
- Mechanism one-liner: 9p21 (CDKN2A/MTAP) methylation/silencing tracks CDK4/6-inhibitor
  (Palbociclib) sensitivity.

## (2) Novel hit A — DAPK1 → Luminespib  [recommended novel demo]
- **Luminespib top marker = COBL** (rho=-0.438, q=2.1e-42); split p=7.8e-31, effect
  0.044, grade neutral.
- **DAPK1** (the gate's candidate) ranks **87 / 11385**, rho=**-0.346**, q=3.2e-26 —
  negative, i.e. DAPK1 methylation marks Luminespib SENSITIVITY.
- Mechanism one-liner: DAPK1 promoter hypermethylation (a canonical tumor-suppressor
  silencing event) marks sensitivity to the HSP90 inhibitor Luminespib — a methylation-
  defined responder stratum that cuts across histology.
- Why DAPK1 over COBL for the demo: DAPK1 promoter methylation is a well-established,
  interpretable epigenetic silencing event; COBL is the strongest correlate but less
  mechanistically legible.

## (3) Novel hit B — FKBP1A → Vorinostat  [strongest correlation overall]
- **Vorinostat top marker = FKBP1A**, rho = **+0.607**, q = **1.3e-93** (the single
  strongest marker→drug correlation in the whole scan). Positive rho → high FKBP1A
  methylation marks RESISTANCE to Vorinostat. Split p=4.3e-61, effect 0.049, neutral.
- The gate's candidate **RASSF1** ranks 406/11904, rho=+0.419, q=1.4e-40 (also a
  resistance-direction marker; RASSF1 silencing is a classic methylation event).
- Mechanism one-liner: FKBP1A (and RASSF1) methylation state stratifies response to the
  HDAC inhibitor Vorinostat — an epigenetic-marker × epigenetic-drug axis.

---

## Recommended demo trio
1. **CDKN2A/MTAP → Palbociclib** (positive control, canonical CDK4/6 biology, 9p21).
2. **DAPK1 → Luminespib** (novel, most interpretable: tumor-suppressor silencing → HSP90i sensitivity).
3. (optional) **FKBP1A/RASSF1 → Vorinostat** (novel, strongest correlation; epigenetic-marker × HDACi).

Note: the gate (blind-clustering) and this marker method disagree on the *single top*
gene per drug (MTAP vs CDKN2A; COBL vs DAPK1; FKBP1A vs RASSF1) — but the marker method
places every gate candidate well inside the significant top-tier, and surfaces the same
biological loci. The marker mechanism is the right choice for the demo.
