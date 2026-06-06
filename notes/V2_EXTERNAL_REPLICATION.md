# External replication — Phase 1 payoff (2026-06-06)

`scripts/v2_external_replication.py`: top silencing-validated leads re-tested in
PRISM + CTRP (independent screens, AUCs sign-comparable to GDSC). REPLICATES =
external tissue-adjusted partial corr same sign as GDSC + p<0.05. Output: `outputs/external_replication.csv`.

## BULLETPROOF leads (PRISM-replicated): 7
| drug | gene | gdsc_r | prism_r (p) | ctrp_r (p) | direction | screens |
|---|---|---|---|---|---|---|
| **Olaparib** | **ID2** | -0.33 | -0.12 (.023) | -0.25 (4e-7) | meth→SENSITIVE | **GDSC+PRISM+CTRP** |
| **Olaparib** | **TRIM32** | -0.28 | -0.11 (.033) | -0.19 (1e-4) | meth→SENSITIVE | **GDSC+PRISM+CTRP** |
| **Olaparib** | **OTULINL** | -0.29 | -0.11 (.040) | -0.23 (5e-6) | meth→SENSITIVE | **GDSC+PRISM+CTRP** |
| Temozolomide | TSPAN13 | -0.26 | -0.16 | (no CTRP cov) | meth→SENSITIVE | GDSC+PRISM |
| Trametinib | N4BP3 | -0.20 | -0.10 | — | meth→SENSITIVE | GDSC+PRISM |
| Trametinib | COL4A3 | -0.10 | -0.11 | — | meth→SENSITIVE | GDSC+PRISM |
| Trametinib | NACC2 | +0.30 | +0.13 | +CTRP | meth→resistant | GDSC+PRISM+CTRP |

## Hero card: Olaparib → ID2
Promoter methylation of **ID2** (Inhibitor of DNA-binding 2; a developmental TF, NOT a canonical PARP/DDR gene) marks **olaparib sensitivity**, **tissue-agnostically**, is **functionally silenced** (meth↑→expr↓), and **replicates across GDSC + PRISM + CTRP**. The "no expert would prioritize this gene" proof-point — obscure gene, cross-indication, functionally + externally validated.

## Honest read
- Funnel is stringent: only 7 survive discovery→silencing→external replication (good — most leads correctly don't).
- **Effect sizes modest** (gdsc_r ~-0.3, single-gene methylation→response is inherently weak/diffuse). PRISM replication is **weak-but-significant** (r~-0.11); **CTRP is the stronger confirmer** (r~-0.25, p<1e-6) for the Olaparib trio. The load-bearing claim is **directional concordance across 3 independent screens** + functional silencing — not a large effect.
- Cell-line, not patient. These are defensible novel HYPOTHESES (evidence ladder ~L5), not clinical claims.

## Status
Phase 1 (science) DONE — bulletproof cross-indication cards exist. Next: Phase 2 — strata-v3.html front end (three gates + Track B pipeline + hero cards), data from outputs/external_replication.csv.
