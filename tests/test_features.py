from strata.engine.features import select_variable_features


def test_selects_most_variable(synthetic_meth):
    meth, _ = synthetic_meth
    # Signal CpGs are cg000-cg009 (group-1 hypermethylated).
    # cg009 gets clipped to 1.0 in many group-1 samples (lower observed variance).
    # Use n_top=20 and verify >=9 of the 10 signal CpGs are recovered.
    feats = select_variable_features(meth, n_top=20)
    assert len(feats) == 20
    signal_cpgs = {f"cg{i:03d}" for i in range(10)}
    recovered = signal_cpgs & set(feats)
    assert len(recovered) >= 9, f"Expected >=9 signal CpGs in top-20, got {recovered}"
