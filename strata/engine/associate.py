from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy.stats import kruskal, spearmanr
from statsmodels.stats.multitest import multipletests
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test


@dataclass
class AssociationResult:
    effect_size: float
    pvalue: float
    responder_label: int
    grade: str


def score_responder(auc_by_cluster: dict[int, float], baseline: float) -> dict[int, str]:
    grades = {}
    for label, auc in auc_by_cluster.items():
        delta = baseline - auc
        if delta > 0.20:
            grades[label] = "++"
        elif delta > 0.07:
            grades[label] = "+"
        elif delta < -0.20:
            grades[label] = "--"
        elif delta < -0.07:
            grades[label] = "-"
        else:
            grades[label] = "neutral"
    return grades


def test_drug_association(clusters: pd.Series, drug_auc: pd.Series) -> AssociationResult:
    common = clusters.index.intersection(drug_auc.index)
    clusters, drug_auc = clusters.loc[common], drug_auc.loc[common]
    keep = drug_auc.notna()
    clusters, drug_auc = clusters[keep], drug_auc[keep]
    if clusters.nunique() < 2:
        return AssociationResult(effect_size=0.0, pvalue=1.0,
                                 responder_label=int(clusters.iloc[0]) if len(clusters) else 0,
                                 grade="neutral")
    groups = [drug_auc[clusters == c].values for c in sorted(clusters.unique())]
    stat, p = kruskal(*groups)
    means = {int(c): float(drug_auc[clusters == c].mean()) for c in clusters.unique()}
    baseline = float(drug_auc.median())
    grades = score_responder(means, baseline)
    responder = min(means, key=means.get)
    effect = baseline - means[responder]
    return AssociationResult(
        effect_size=effect,
        pvalue=float(p),
        responder_label=responder,
        grade=grades[responder],
    )


def test_survival_association(
    clusters: pd.Series, time: pd.Series, event: pd.Series
) -> dict:
    common = clusters.index.intersection(time.index).intersection(event.index)
    clusters, time, event = clusters.loc[common], time.loc[common], event.loc[common]
    lr = multivariate_logrank_test(time, clusters, event)
    km = {}
    for c in sorted(clusters.unique()):
        m = clusters == c
        kmf = KaplanMeierFitter().fit(time[m], event[m])
        km[int(c)] = (kmf.timeline.tolist(), kmf.survival_function_.iloc[:, 0].tolist())
    return {"logrank_p": float(lr.p_value), "km": km}


def rank_response_markers(meth: pd.DataFrame, drug_auc: pd.Series, min_n: int = 30) -> pd.DataFrame:
    """For each gene's methylation, Spearman-correlate against drug AUC.
    Returns DataFrame[gene, rho, pvalue, qvalue, n] sorted by |rho| desc among q<0.1.
    Lower AUC = more sensitive, so a negative rho = methylation marks sensitivity."""
    common = meth.index.intersection(drug_auc.index)
    m, a = meth.loc[common], drug_auc.loc[common]
    a = a.astype(float)
    keep = a.notna()
    m, a = m.loc[keep], a.loc[keep]
    n = len(a)
    rows = []
    vals = a.values
    for gene in m.columns:
        rho, p = spearmanr(m[gene].values, vals)
        rows.append((gene, float(rho), float(p), n))
    df = pd.DataFrame(rows, columns=["gene", "rho", "pvalue", "n"]).dropna(subset=["pvalue"])
    df["qvalue"] = multipletests(df["pvalue"], method="fdr_bh")[1]
    sig = df[df["qvalue"] < 0.1].copy()
    sig = sig.reindex(sig["rho"].abs().sort_values(ascending=False).index).reset_index(drop=True)
    return sig


def discover_responder_subgroup(meth: pd.DataFrame, drug_auc: pd.Series, marker: str):
    """Median-split lines on `marker` methylation; return (clusters, AssociationResult).
    clusters: 1 = high methylation, 0 = low. Reuses test_drug_association."""
    common = meth.index.intersection(drug_auc.index)
    g = meth.loc[common, marker]
    clusters = (g > g.median()).astype(int)
    clusters.name = "cluster"
    assoc = test_drug_association(clusters, drug_auc.loc[common])
    return clusters, assoc
