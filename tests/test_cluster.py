from strata.engine.features import select_variable_features
from strata.engine.cluster import cluster_samples
from sklearn.metrics import adjusted_rand_score


def test_recovers_latent_groups(synthetic_meth):
    meth, truth = synthetic_meth
    feats = select_variable_features(meth, n_top=10)
    labels = cluster_samples(meth, feats, k=2, random_state=0)
    assert adjusted_rand_score(truth, labels.loc[truth.index]) > 0.8
