import pandas as pd
from strata.engine.annotate import differential_cpgs


def test_finds_driver_cpgs(synthetic_meth):
    meth, truth = synthetic_meth
    diff = differential_cpgs(meth, truth, target_label=1)
    top10 = set(diff.head(10)["cpg"])
    assert top10 == {f"cg{i:03d}" for i in range(10)}
