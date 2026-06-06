import numpy as np, pandas as pd
from strata.engine.associate import rank_response_markers


def test_top_marker_is_the_correlated_gene():
    rng = np.random.default_rng(0)
    n = 80
    auc = pd.Series(rng.uniform(0, 1, n), index=[f"S{i}" for i in range(n)])
    meth = pd.DataFrame(rng.uniform(0, 1, (n, 50)),
                        index=auc.index, columns=[f"G{j}" for j in range(50)])
    # make G7 strongly anti-correlated with auc (high methylation -> low auc)
    meth["G7"] = np.clip(1 - auc.values + rng.normal(0, 0.05, n), 0, 1)
    res = rank_response_markers(meth, auc)
    assert res.iloc[0]["gene"] == "G7"
    assert res.iloc[0]["rho"] < 0
