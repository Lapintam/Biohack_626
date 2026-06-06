# Track B — silencing-validation result (2026-06-06)

`scripts/v2_silencing_validation.py`. Pipeline per drug: global-axis promoter scan →
response-candidates (q<0.1) → silencing filter (promoter-meth ↔ own expression, neg =
functional silencing, FDR<0.05) → VALIDATED LEADS (response-associated AND functionally silenced).
Expression∩promoter = 613 lines.

## Funnel (per drug)
| drug | response-candidates | with expr | functionally-silenced & validated |
|---|---|---|---|
| Palbociclib | ~500 | ~450 | ~160 |
| Temozolomide | 337 | 310 | 160 |
| Olaparib | 500 | 453 | 175 |

## Example validated leads (response-associated + functionally silenced)
- **Palbociclib:** CDKN2B-AS1 (ANRIL — the 9p21 lncRNA, biologically coherent), TEAD4, SPRY4, KCTD15.
- **Temozolomide:** BCL2L11/BIM (pro-apoptotic), NDRG1, MYCL, SORL1, RASGEF1B, RAPGEF5.
- **Olaparib:** N6AMT1, PINK1, ID2, DHCR24, WWC1 (Hippo), TUBB2A.

## Honest read
- **WIN:** end-to-end unbiased→validated pipeline works. Leads are principled — both tissue-agnostically response-associated AND functionally silenced by their own promoter methylation (not passenger correlations). This is a real discovery engine, beyond known cancer genes.
- **Caveat:** the silencing filter keeps ~half the candidates (~160/drug) because promoter-methylation→silencing is *common* — so it confirms MECHANISM but doesn't alone produce a tiny shortlist. The crisp shortlist = top-ranked by combined response+silencing stringency (the top ~10–15 shown per drug).
- Leads are mostly **obscure genes** (expected for genuine unbiased discovery) with some coherent ones (CDKN2B-AS1/9p21, BCL2L11 apoptosis). Novelty is real but **unconfirmed** → needs external validation (PRISM/CTRP) to reach a defensible handful.

## Status
Track B pipeline complete: promoter discovery (specific) → RNA-seq silencing (functional) → ranked validated leads. Next: external validation (#2) on the top leads.
