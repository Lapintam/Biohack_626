import numpy as np
import pandas as pd
import pytest

@pytest.fixture
def synthetic_meth():
    """50 samples x 100 CpGs; two latent groups differing at 10 CpGs."""
    rng = np.random.default_rng(0)
    n, p = 50, 100
    base = rng.uniform(0.1, 0.9, size=(n, p))
    group = np.array([0] * 25 + [1] * 25)
    base[group == 1, :10] += 0.4  # group-1 hypermethylated at first 10 CpGs
    base = np.clip(base, 0, 1)
    samples = [f"S{i}" for i in range(n)]
    cpgs = [f"cg{i:03d}" for i in range(p)]
    return pd.DataFrame(base, index=samples, columns=cpgs), pd.Series(group, index=samples)
