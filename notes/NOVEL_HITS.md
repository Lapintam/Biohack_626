# Novel responder-subgroup scan (Rung-A unsupervised pipeline)

Date: 2026-06-04. Pipeline: top-2000 variable gene-level methylation features →
auto-k KMeans (silhouette) → Kruskal association of cluster vs drug AUC →
differential drivers (Mann-Whitney on Δβ).

## Critical caveat (read first)

Auto-k clustering of the GDSC methylome yields **only 2 clusters: 164 vs 891 lines**.
The minority cluster (cluster 0, n≈150–164) is **the same subgroup for every drug** and
its top differential drivers are an identical pan-leukocyte / hematopoietic signature:
`TJP1, INPP5D, FERMT3, PPIC, KLF4, DLC1, SPN, PTPRCAP, CD38`. INPP5D (SHIP1), FERMT3
(kindlin-3), SPN (CD43), PTPRCAP, CD38 are all immune/blood-lineage markers.

**Conclusion: the unsupervised pipeline is rediscovering ONE lineage axis (blood /
hematologic cell lines), not drug-specific methylation programs.** Every "hit" below is
the blood-cluster being sensitive to a given drug. The drivers are NOT drug-specific —
they are the same lineage barcode each time. This is honest signal (blood lines really
are more/less sensitive to these agents) but it is NOT a methylation-defined cross-
indication stratum invisible to current taxonomy. It is lineage, which taxonomy already sees.

## Candidate hits (significant p<0.05, grade + responder), ranked by AUC effect

1. **Venetoclax** — grade `+`, p=1.5e-3, AUC effect=0.143, responder n=150 (cluster 0).
   Drivers: TJP1, INPP5D, FERMT3, PPIC, KLF4, DLC1.
   *Biologically expected:* Venetoclax is a BCL2 inhibitor used clinically in blood
   cancers; cluster 0 = blood lines. This is a correct-but-known lineage effect, largest effect size in the scan.

2. **Nutlin-3a (-)** — grade `+`, p=8.2e-6, AUC effect=0.091, responder n=152 (cluster 0).
   Drivers: TJP1, INPP5D, FERMT3, PPIC, KLF4, DLC1.
   MDM2 inhibitor; blood-line cluster again. Same lineage barcode.

3. **Trametinib** — grade `+`, p=2.3e-11, AUC effect=0.074, responder n=814 (cluster 1,
   the MAJORITY / non-blood cluster). Drivers: same lineage set.
   MEK inhibitor; here the non-blood cluster is the responder. Still the same 2-way split.

## Other drugs tested (not clean + hits)
- Sorafenib: p=4.2e-20 but grade=neutral (effect 0.056) — significant split, small effect.
- Erlotinib, Olaparib, Alpelisib: p<0.05 but grade=neutral, tiny effects (0.02–0.03).
- Dasatinib: ns (p=0.49). Lapatinib: ns (p=0.93).

## Implication for the demo
A genuinely "novel, taxonomy-invisible" methylation stratum requires finer-grained
clustering than a 2-way lineage split — e.g. higher fixed k, sub-clustering within a
lineage, or a supervised/biomarker discovery step. The current unsupervised Rung-A is a
lineage detector. See diagnostic verdict (scripts/diag_positive_control.py): even the
established CDKN2A→Palbociclib positive control is NOT recovered unsupervised (CDKN2A
ranks 9858/14608 as a driver) — only the targeted median-split recovers it (p=3.7e-18).
