from strata.data import tcga_laml


def test_clinical_has_survival():
    c = tcga_laml.load_tcga_clinical()
    assert {"os_time", "os_event"} <= set(c.columns)
    assert c["os_event"].isin([0, 1]).all()
    assert (c["os_time"] >= 0).all()


def test_methylation_and_clinical_overlap():
    m = tcga_laml.load_tcga_methylation()
    c = tcga_laml.load_tcga_clinical()
    common = m.index.intersection(c.index)
    assert len(common) >= 100  # LAML has ~190+ patients with methylation
    assert ((m.values >= 0) & (m.values <= 1)).mean() > 0.99
