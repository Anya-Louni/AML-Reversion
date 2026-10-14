"""Phase 1 unit tests.

All arrays in this file are MOCK DATA generated for testing only.
They are never used by the real pipeline.
"""
import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

from pipeline.phase1_landscape import landscape as L
from pipeline.phase1_landscape.mapping import LandscapeMapper, lognorm
from pipeline.phase1_landscape.signature import define_signature, spearman_with

rng = np.random.default_rng(0)


def test_lognorm_matches_definition():
    X = np.array([[1, 3, 0], [0, 0, 0], [2, 2, 4]], dtype=float)  # MOCK counts
    out = lognorm(X).toarray()
    exp = np.log1p(X / np.maximum(X.sum(1, keepdims=True), 1) * 1e4)
    np.testing.assert_allclose(out, exp, rtol=1e-5)


@pytest.fixture
def mapper():
    # MOCK: 300 cells on a 1D trajectory embedded in 20 genes
    t = rng.uniform(0, 1, 300)
    W = rng.normal(size=(1, 20))
    X = t[:, None] @ W + rng.normal(scale=0.05, size=(300, 20))
    mean = X.mean(0)
    pca = PCA(5, random_state=0).fit(X - mean)
    Y = np.c_[t, t * 2, -t]
    genes = np.array([f"g{i}" for i in range(20)])
    return LandscapeMapper.fit(X, genes, pca, mean, Y, k=10), X, t, genes


def test_mapper_recovers_score(mapper):
    m, X, t, genes = mapper
    P = m.predict(X)
    assert spearmanr(P.score, t)[0] > 0.95


def test_align_fills_missing_with_mean(mapper):
    m, X, t, genes = mapper
    sub = genes[::-1][:15]  # reversed order, 5 genes missing
    Xq = pd.DataFrame(X, columns=genes)[sub].values
    out, frac = m.align(Xq, sub)
    assert frac == pytest.approx(0.75)
    missing = ~np.isin(m.genes, sub)
    np.testing.assert_allclose(out[:, missing], np.tile(m.mean[missing], (len(X), 1)))
    np.testing.assert_allclose(out[:, ~missing], X[:, ~missing])


def test_to_score_clips():
    pt = np.r_[np.linspace(0, 1, 98), 50.0, -50.0]  # MOCK with outliers
    s, (lo, hi) = L.to_score(pt, (0.01, 0.99))
    assert s.min() == 0 and s.max() == 1
    assert 0.3 < s[49] < 0.7  # middle value stays in the middle despite outliers


def test_order_agreement():
    ct = pd.Series(["HSC", "Prog", "GMP", "Mono", "cDC"])
    assert L.order_agreement(np.array([0, 1, 2, 3, 99.0]), ct, ["HSC", "Prog", "GMP", "Mono"]) == pytest.approx(1.0)


def test_potential_grid_slopes_downhill_with_score():
    # MOCK: score increases along x
    xy = rng.uniform(0, 10, size=(4000, 2))
    score = xy[:, 0] / 10
    g = L.potential_grid(xy, score, {"grid": 60, "margin": 0.05, "bandwidth": 0.05, "weight_pseudotime": 0.9})
    assert g["height"].shape == (60, 60)
    h = L.height_at(g, np.c_[[1, 9], [5, 5]])
    assert h[0] > h[1]  # stem side (low score) is higher
    assert 0 <= g["height"].min() and g["height"].max() <= 1


def test_spearman_with_matches_scipy():
    X = rng.poisson(2, size=(200, 5)).astype(float)  # MOCK
    y = X[:, 0] + rng.normal(size=200)
    out = spearman_with(X, y, [f"g{i}" for i in range(5)], min_frac=0).set_index("gene")
    for i in range(5):
        assert out.loc[f"g{i}", "spearman"] == pytest.approx(spearmanr(X[:, i], y)[0], abs=1e-6)


def test_define_signature_sizes():
    corr = pd.DataFrame({"gene": list("abcdef"), "spearman": [0.9, 0.5, 0.1, -0.1, -0.5, -0.9]})
    s = define_signature(corr, 2)
    assert list(s[s.direction == "up_with_maturation"].gene) == ["a", "b"]
    assert list(s[s.direction == "down_with_maturation"].gene) == ["f", "e"]
