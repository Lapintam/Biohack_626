"""Promoter-vs-gene-level validation (#1 payoff check).

Re-run key marker->drug pairs on BOTH methylation sources:
  - GDSC gene-level imputed  (current)
  - CCLE RRBS promoter/TSS    (the upgrade)
Headline question: does MGMT -> Temozolomide RECOVER at promoter resolution?
(It was null at gene-level: rho=-0.018.) And do the demo markers sharpen?

Run: PYTHONPATH=. uv run python scripts/v2_promoter_validation.py
"""

from __future__ import annotations

import pandas as pd
from strata.data import gdsc
from strata.data import ccle_methylation as ccle
from scripts.v2_mechanism_anchored import eval_gene

PAIRS = [
    ("MGMT", "Temozolomide", "negative control -> RECOVERY test"),
    ("CDKN2A", "Palbociclib", "positive control"),
    ("MTAP", "Palbociclib", "positive control (9p21 neighbor)"),
    ("INPP4B", "Rapamycin", "novel lead A"),
    ("SFRP1", "SB216763", "novel lead B"),
]


def _auc(drug, name):
    s = (drug[drug["drug_name"].str.lower() == name.lower()]
         .drop_duplicates("COSMIC_ID").set_index("COSMIC_ID")["auc"])
    return s[~s.index.duplicated(keep="first")].astype(float)


def main():
    drug = gdsc.load_gdsc_drug_response()
    tissue = gdsc.load_gdsc_annotations()["tissue"]
    m_gene = gdsc.load_gdsc_methylation()
    m_prom = ccle.load_ccle_promoter_methylation()

    print(f"{'gene -> drug':<26}{'source':<10}{'pooled':>8}{'r_adj':>8}{'wsig':>6}{'lvl':>5}{'n':>6}")
    print("-" * 74)
    for gene, dn, label in PAIRS:
        auc = _auc(drug, dn)
        for src, meth in (("gene-lvl", m_gene), ("promoter", m_prom)):
            if gene not in meth.columns:
                print(f"{gene+' -> '+dn:<26}{src:<10}  (gene absent)")
                continue
            e = eval_gene(meth, auc, tissue, gene)
            if e is None:
                print(f"{gene+' -> '+dn:<26}{src:<10}  (insufficient n)")
                continue
            print(f"{gene+' -> '+dn:<26}{src:<10}{e['rho_pooled']:>8.3f}{e['r_adj']:>8.3f}"
                  f"{e['within_sig']:>6}{e['level']:>5}{e['n']:>6}")
        print(f"   ^ {label}")
    print("\n(neg rho/r_adj = methylation marks SENSITIVITY; recovery = MGMT goes clearly negative)")


if __name__ == "__main__":
    main()
