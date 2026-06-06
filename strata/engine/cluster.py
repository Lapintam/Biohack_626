import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


def cluster_samples(
    meth: pd.DataFrame,
    features: list[str],
    k: int | None = None,
    random_state: int = 0,
) -> pd.Series:
    X = meth[features].values
    if k is None:
        best_k, best_s = 2, -1.0
        for cand in range(2, 7):
            lab = KMeans(n_clusters=cand, random_state=random_state, n_init=10).fit_predict(X)
            s = silhouette_score(X, lab)
            if s > best_s:
                best_k, best_s = cand, s
        k = best_k
    labels = KMeans(n_clusters=k, random_state=random_state, n_init=10).fit_predict(X)
    return pd.Series(labels, index=meth.index, name="cluster")
