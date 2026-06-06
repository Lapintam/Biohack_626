from strata.data import gdsc


def test_methylation_shape():
    m = gdsc.load_gdsc_methylation()
    assert m.shape[0] > 500 and m.shape[1] > 1000
    assert m.index.dtype == object               # COSMIC_ID as str
    assert ((m.values >= 0) & (m.values <= 1)).mean() > 0.99  # betas


def test_drug_response_columns():
    d = gdsc.load_gdsc_drug_response()
    assert {"COSMIC_ID", "drug_name", "ln_ic50", "auc"} <= set(d.columns)
    assert (d["drug_name"].str.lower() == "temozolomide").any()


def test_annotations():
    a = gdsc.load_gdsc_annotations()
    assert {"cell_line_name", "tissue"} <= set(a.columns)
    assert a.index.dtype == object
    assert len(a) > 500


def test_mgmt_methylation_feature_present():
    m = gdsc.load_gdsc_methylation()
    assert "MGMT" in m.columns
    mgmt = m["MGMT"]
    assert ((mgmt >= 0) & (mgmt <= 1)).mean() > 0.99
