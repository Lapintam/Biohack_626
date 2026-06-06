# Promoter-level validation (#1 payoff) — 2026-06-06

`scripts/v2_promoter_validation.py`: re-ran key marker→drug pairs on GDSC gene-level (imputed) vs CCLE RRBS promoter methylation.

| gene → drug | gene-level (pooled / r_adj / lvl) | promoter (pooled / r_adj / lvl) |
|---|---|---|
| MGMT → Temozolomide | −0.014 / +0.078 / null | **−0.117 / −0.069 / weak** (direction recovers, stays weak) |
| CDKN2A → Palbociclib | −0.303 / −0.164 / **L3** | −0.270 / −0.059 / L1 (pooled holds, tissue-agnostic collapses) |
| MTAP → Palbociclib | **−0.461** / −0.196 / **L3** | **−0.065** / +0.159 / L1 (signal vanishes) |
| INPP4B → Rapamycin | −0.166 / −0.213 / L3 | −0.127 / −0.042 / L1 (weakens) |
| SFRP1 → SB216763 | −0.183 / −0.169 / L3 | −0.098 / −0.106 / L2 (weakens, holds partly) |

## KEY FINDING (load-bearing, honest)
The gene-level **imputed** methylation matrix **inflated the methylation→response signal**, especially at **9p21 (MTAP/CDKN2A)** — a **deletion** locus where palbociclib sensitivity is driven by **copy-number loss, not promoter methylation**. The imputed "methylation" value at a deleted gene tracks the deletion, so the gene-level positive control was **substantially a CNV proxy**. At true promoter (RRBS) resolution the 9p21→palbociclib signal largely vanishes (MTAP −0.461→−0.065). R2 had flagged the MTAP/CNV caveat.

- MGMT→TMZ **partially recovers** (direction correct, ~8× stronger pooled) but stays weak — matches the literature (messy in cell lines).
- Novel leads (INPP4B, SFRP1) **weaken** at promoter resolution (SFRP1 holds partly at L2).
- Caveat on the comparison: promoter n is smaller (~277–421 vs ~840) due to RRBS coverage/NaN, which reduces tissue-adjusted power — but the **pooled** collapses too (MTAP), which power alone can't explain → CNV-proxy interpretation stands.

## Implications
1. **Promoter-level is the honest substrate; the truthful signal is weaker.** Gene-level results overstate methylation's role.
2. **Need copy-number** (multi-omic, roadmap #3) to disambiguate methylation-silencing from deletion. **RNA-seq** (also available, same lines) would confirm silencing (expression↓) vs CNV.
3. **External validation (#2)** now even more important for the surviving leads.
4. Decision (demo-day): keep gene-level demo WITH this caveat, or reframe around the honest promoter picture + the "we caught our own CNV artifact" rigor narrative.
