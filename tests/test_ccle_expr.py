from strata.data import ccle_expression as ce


def test_expression_loads():
    e = ce.load_ccle_expression()
    # Real shape: 733 lines x 57138 unique gene symbols (CCLE 19Q4 full).
    assert e.shape[0] > 600 and e.shape[1] > 10000
    assert e.index.dtype == object


def test_demo_genes_present():
    e = ce.load_ccle_expression()
    for g in ["MGMT", "CDKN2A", "MTAP", "INPP4B", "SFRP1"]:
        assert g in e.columns
