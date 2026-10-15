"""Experiment 3 unit tests: dataset parsing, simulated change vector, and the F1 metrics.

All arrays and tables in this file are MOCK DATA generated for testing only.
They are never used by the real pipeline.
"""
import numpy as np
import pandas as pd
import pytest

from pipeline.exp3 import datasets as DS
from pipeline.exp3 import f1 as F

rng = np.random.default_rng(0)
GENES = [f"g{i}" for i in range(60)]


def _sim(n_drugs=40, norm=10.0):
    D = pd.DataFrame(rng.normal(size=(n_drugs, len(GENES))), index=[f"d{i}" for i in range(n_drugs)], columns=GENES)

    class M:
        genes = np.array(GENES)
    return F.Simulator(M, D, np.ones(len(GENES)), set(GENES), norm, 1.0), D


def test_log2cpm_matches_definition():
    X = pd.DataFrame({"a": [10.0, 0.0, 90.0], "b": [50.0, 50.0, 0.0]}, index=list("xyz"))
    out = DS.log2cpm(X)
    exp = np.log2(X / X.sum(axis=0) * 1e6 + 1)
    np.testing.assert_allclose(out.values, exp.values)
    assert out.loc["y", "a"] == 0.0            # zero counts -> log2(1) = 0


def test_real_change_is_difference_of_means():
    X = pd.DataFrame(rng.integers(1, 500, size=(20, 4)).astype(float),
                     index=[f"g{i}" for i in range(20)], columns=["t1", "t2", "c1", "c2"])
    meta = pd.DataFrame({"sample": ["t1", "t2", "c1", "c2"], "arm": ["ATRA", "ATRA", "control", "control"]})
    r = DS.real_change(X, meta, "ATRA")
    L = DS.log2cpm(X)
    np.testing.assert_allclose(r.values, (L[["t1", "t2"]].mean(axis=1) - L[["c1", "c2"]].mean(axis=1)).values)


def test_real_change_drops_all_zero_genes():
    X = pd.DataFrame({"t1": [5.0, 0.0], "c1": [5.0, 0.0]}, index=["keep", "drop"])
    meta = pd.DataFrame({"sample": ["t1", "c1"], "arm": ["ATRA", "control"]})
    assert list(DS.real_change(X, meta, "ATRA").index) == ["keep"]


def test_simulated_change_has_fixed_length_and_direction():
    sim, D = _sim(norm=7.0)
    c = sim.change("d0")
    assert np.linalg.norm(c.values) == pytest.approx(7.0)
    # direction follows the signature when sd weighting is uniform
    assert F.cosine(c.values, D.loc["d0"].values) == pytest.approx(1.0, abs=1e-6)


def test_change_matrix_matches_change():
    sim, _ = _sim()
    V = sim.change_matrix(["d1", "d2"])
    np.testing.assert_allclose(V.loc["d1"].values, sim.change("d1").values)


def test_direction_test_ranks_true_drug_first_on_its_own_signal():
    sim, D = _sim(n_drugs=50)
    real = sim.change("d7") * 3.0          # real change proportional to drug d7's simulated change
    out = F.direction_test(sim, real, "d7", list(D.index), n_perm=200)
    assert out["cosine"] == pytest.approx(1.0, abs=1e-6)
    assert out["rank_in_null"] == 1 and out["p_rank"] <= 0.05
    assert out["perm_p"] < 0.05


def test_direction_test_null_drug_is_not_significant():
    sim, D = _sim(n_drugs=200)
    real = pd.Series(rng.normal(size=len(GENES)), index=GENES)
    ranks = [F.direction_test(sim, real, d, list(D.index), n_perm=50)["p_rank"] for d in D.index[:12]]
    assert np.mean([r <= 0.05 for r in ranks]) < 0.3       # unrelated drugs rarely pass


def test_signature_shift_sign():
    sig = pd.DataFrame({"gene": ["g0", "g1", "g2", "g3"],
                        "direction": ["up_with_maturation"] * 2 + ["down_with_maturation"] * 2})
    ch = pd.Series({"g0": 2.0, "g1": 2.0, "g2": -1.0, "g3": -1.0})
    assert F.signature_shift(ch, sig) == pytest.approx(3.0)
    assert F.signature_shift(-ch, sig) == pytest.approx(-3.0)


def test_maturity_test_detects_agreement():
    sim, _ = _sim()
    sig = pd.DataFrame({"gene": GENES[:10] + GENES[10:20],
                        "direction": ["up_with_maturation"] * 10 + ["down_with_maturation"] * 10})
    real = sim.change("d0") * 2.0
    out = F.maturity_test(sim, real, "d0", sig)
    assert out["same_direction"]
    out_opp = F.maturity_test(sim, -real, "d0", sig)
    assert not out_opp["same_direction"]


def test_trajectory_on_exact_line_has_no_off_axis_component():
    sim, _ = _sim()
    u = sim.change("d3")
    reals = {"t1": u * 0.5, "t2": u * 1.0, "t3": u * 2.0}
    tr = F.trajectory_test(sim, reals, "d3")
    assert tr.off_axis_fraction.max() < 1e-6
    assert tr.parallel.is_monotonic_increasing


def test_trajectory_detects_off_axis_movement():
    sim, _ = _sim()
    u = sim.change("d3")
    perp = pd.Series(rng.normal(size=len(GENES)), index=GENES)
    perp = perp - u * (perp @ u) / (u @ u)          # orthogonalize
    tr = F.trajectory_test(sim, {"t": u + perp * (np.linalg.norm(u) / np.linalg.norm(perp))}, "d3")
    assert tr.off_axis_fraction.iloc[0] == pytest.approx(np.sqrt(0.5), abs=0.02)


# ---------------------------------------------------------------- F2
from pipeline.exp3 import f2 as F2  # noqa: E402


def test_reduced_rank_map_recovers_a_linear_relation():
    n, p, k = 60, 25, 2
    Z = rng.normal(size=(n, k))
    A, B = rng.normal(size=(k, p)), rng.normal(size=(k, p))
    Xs, Ys = Z @ A, Z @ B                       # MOCK: both sides share a rank-2 latent cause
    m = F2.ReducedRankMap(rank=k, alpha=0.01).fit(Xs, Ys)
    pred = m.predict(Xs)
    cos = [F.cosine(pred[i], Ys[i]) for i in range(n)]
    assert np.mean(cos) > 0.9


def test_reduced_rank_map_predict_shape_for_rank_one():
    Xs, Ys = rng.normal(size=(30, 12)), rng.normal(size=(30, 12))
    out = F2.ReducedRankMap(rank=1, alpha=1.0).fit(Xs, Ys).predict(Xs[:3])
    assert out.shape == (3, 12)


def test_select_map_never_sees_excluded_drug():
    idx = [f"d{i}" for i in range(20)]
    cols = [f"g{i}" for i in range(8)]
    Xs = pd.DataFrame(rng.normal(size=(20, 8)), index=idx, columns=cols)
    Ys = pd.DataFrame(rng.normal(size=(20, 8)), index=idx, columns=cols)
    Ys.loc["d0"] = 1e6                          # a wild value that would distort any fit that saw it
    a, _ = F2.select_map(Xs, Ys, [2], [1.0], exclude=["d0"])
    b, _ = F2.select_map(Xs.drop(index="d0"), Ys.drop(index="d0"), [2], [1.0], exclude=[])
    assert a["cv_cosine"] == pytest.approx(b["cv_cosine"])


def test_time_interpolate_is_linear_and_bracketed():
    a = pd.Series({"g1": 0.0, "g2": 10.0})
    b = pd.Series({"g1": 4.0, "g2": 0.0})
    out = F2.time_interpolate({4.0: a, 96.0: b}, 50.0)
    w = (50.0 - 4.0) / (96.0 - 4.0)
    assert out["g1"] == pytest.approx(4.0 * w) and out["g2"] == pytest.approx(10.0 * (1 - w))
    with pytest.raises(ValueError):
        F2.time_interpolate({4.0: a, 96.0: b}, 200.0)


def test_consensus_excluding_leaves_out_that_line():
    obs = pd.DataFrame({"drug": ["x", "x", "y"], "cell_iname": ["HL60", "THP1", "HL60"]})
    X = np.array([[10.0, 10.0], [2.0, 2.0], [5.0, 5.0]])
    g = ["g1", "g2"]
    ex = F2.consensus_excluding(obs, X, g, "HL60")
    assert list(ex.index) == ["x"] and ex.loc["x", "g1"] == pytest.approx(2.0)
    hl = F2.per_cell_line_signatures(obs, X, g, "HL60")
    assert hl.loc["x", "g1"] == pytest.approx(10.0)
