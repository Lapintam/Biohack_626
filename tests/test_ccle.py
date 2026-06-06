import numpy as np

from strata.data import ccle_methylation as cm


def test_promoter_methylation_loads():
    m = cm.load_ccle_promoter_methylation()
    # Real shape is ~614 lines x ~17,181 genes (see notes/CCLE_RRBS_SCHEMA.md).
    assert m.shape[0] > 300 and m.shape[1] > 3000   # lines x genes
    assert m.index.dtype == object                  # COSMIC_ID str


def test_demo_genes_present():
    m = cm.load_ccle_promoter_methylation()
    for g in ["MGMT", "CDKN2A"]:
        assert g in m.columns


def test_values_are_betas():
    m = cm.load_ccle_promoter_methylation()
    vals = m.values
    # Betas in [0, 1] where present; NaN allowed for uncovered windows.
    present = vals[~np.isnan(vals)]
    assert ((present >= 0) & (present <= 1)).mean() > 0.99


def test_drop_in_orientation_matches_gdsc():
    # Same orientation/keys as the GDSC gene-level loader: rows=COSMIC_ID, cols=gene.
    from strata.data import gdsc
    m = cm.load_ccle_promoter_methylation()
    g = gdsc.load_gdsc_methylation()
    assert m.index.name == g.index.name == "COSMIC_ID"
    assert len(set(m.index) & set(g.index)) > 100
