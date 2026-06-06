import numpy as np
import pandas as pd
from strata.engine.associate import test_drug_association as _test_drug_association
from strata.engine.associate import AssociationResult


def test_finds_responder_cluster():
    clusters = pd.Series([0, 0, 0, 1, 1, 1], index=list("abcdef"))
    auc = pd.Series([0.1, 0.15, 0.12, 0.8, 0.85, 0.82], index=list("abcdef"))
    res = _test_drug_association(clusters, auc)
    assert isinstance(res, AssociationResult)
    assert res.responder_label == 0
    assert res.pvalue < 0.05
    assert res.grade in {"++", "+"}


def test_drug_association_drops_nan_auc():
    clusters = pd.Series([0, 0, 1, 1], index=list("abcd"))
    auc = pd.Series([0.1, 0.15, 0.8, np.nan], index=list("abcd"))
    res = _test_drug_association(clusters, auc)
    assert res.pvalue == res.pvalue  # not NaN
    assert res.responder_label == 0
