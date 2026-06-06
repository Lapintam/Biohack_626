import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def differential_cpgs(
    meth: pd.DataFrame, clusters: pd.Series, target_label: int
) -> pd.DataFrame:
    common = meth.index.intersection(clusters.index)
    meth, clusters = meth.loc[common], clusters.loc[common]
    in_grp = clusters == target_label
    rows = []
    for cpg in meth.columns:
        a, b = meth.loc[in_grp, cpg], meth.loc[~in_grp, cpg]
        try:
            _, p = mannwhitneyu(a, b)
        except ValueError:
            p = 1.0
        rows.append((cpg, float(a.mean() - b.mean()), float(p)))
    return (
        pd.DataFrame(rows, columns=["cpg", "delta_beta", "pvalue"])
        .sort_values("delta_beta", key=lambda s: s.abs(), ascending=False)
        .reset_index(drop=True)
    )


def annotate_cpgs(cpgs: list[str]) -> pd.DataFrame:
    """Feature columns in the GDSC gene-level matrix ARE gene symbols, so the gene
    is the feature id itself."""
    return pd.DataFrame({"cpg": cpgs, "gene": cpgs, "region": ["gene-level"] * len(cpgs)})
