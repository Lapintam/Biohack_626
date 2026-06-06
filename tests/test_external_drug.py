from strata.data import external_drug

_DEMO_DRUGS = [
    "palbociclib",
    "trametinib",
    "temozolomide",
    "olaparib",
    "dabrafenib",
    "dasatinib",
]


def test_prism_loads_long_schema():
    d = external_drug.load_prism_drug_response()
    assert list(d.columns) == ["COSMIC_ID", "drug_name", "auc"]
    # collapsed to one row per (line, drug); GDSC-scale signal expected.
    assert len(d) > 100_000
    assert d["COSMIC_ID"].dtype == object
    assert d["COSMIC_ID"].map(type).eq(str).all()
    assert d["drug_name"].str.islower().all()


def test_prism_has_demo_drug():
    d = external_drug.load_prism_drug_response()
    names = d["drug_name"]
    assert any(names.str.contains(drug, case=False, regex=False).any() for drug in _DEMO_DRUGS)


def test_ctrp_loads_long_schema():
    d = external_drug.load_ctrp_drug_response()
    assert list(d.columns) == ["COSMIC_ID", "drug_name", "auc"]
    assert len(d) > 300
    assert d["COSMIC_ID"].map(type).eq(str).all()
    names = d["drug_name"]
    assert any(names.str.contains(drug, case=False, regex=False).any() for drug in _DEMO_DRUGS)
