import pandas as pd


def select_variable_features(meth: pd.DataFrame, n_top: int = 2000) -> list[str]:
    variances = meth.var(axis=0).sort_values(ascending=False)
    return variances.head(n_top).index.tolist()
