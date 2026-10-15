"""Unit tests for Phases 2-5 (scoring functions and evaluation metrics).

All arrays in this file are MOCK DATA generated for testing only.
They are never used by the real pipeline.
"""
import numpy as np
import pandas as pd
import pytest

from pipeline.phase2_patients import deconvolution as D
from pipeline.phase2_patients import placement as P
from pipeline.phase2_patients import stuck as S
from pipeline.phase3_drugs import scoring as SC
from pipeline.phase4_validation import models as M
from pipeline.phase5_simulation_export import simulate as SIM

rng = np.random.default_rng(0)


# ---------------------------------------------------------------- phase 2
def test_maturity_threshold_separates():
    score = np.r_[rng.uniform(0, 0.3, 200), rng.uniform(0.6, 1, 200)]
    ct = pd.Series(["HSC"] * 200 + ["Mono"] * 200)
    t = P.maturity_threshold(score, ct, ["HSC"], ["Mono"])
    assert 0.3 <= t["threshold"] <= 0.6 and t["balanced_accuracy"] == pytest.approx(1.0)


def test_nnls_recovers_mixture():
    ref = pd.DataFrame(rng.gamma(2, 50, size=(300, 4)), columns=list("abcd"), index=[f"g{i}" for i in range(300)])
    true = np.array([0.5, 0.3, 0.2, 0.0])
    bulk = pd.DataFrame({"s": ref.values @ true}, index=ref.index)
    f = D.nnls_fractions(bulk, ref).iloc[0].values
    np.testing.assert_allclose(f, true, atol=1e-6)


def test_malignant_shares_renormalize():
    fr = pd.DataFrame({"A-like": [0.2], "B-like": [0.6], "T_NK": [0.2]})
    sh = D.malignant_shares(fr, ["A-like", "B-like"])
    assert sh.iloc[0].tolist() == pytest.approx([0.25, 0.75])


def test_bulk_mature_expectation_linear():
    m = pd.Series(np.linspace(0, 0.8, 50), index=[f"s{i}" for i in range(50)])
    expr = pd.DataFrame({"up": 2 + 3 * m.values, "down": 5 - 4 * m.values}, index=m.index)
    fit = S.bulk_mature_expectation(expr, m)
    assert fit.loc["up", "expected_mature"] == pytest.approx(5.0)
    assert fit.loc["down", "expected_mature"] == pytest.approx(1.0)


# ---------------------------------------------------------------- phase 3
def test_cosine_and_connectivity_signs():
    q = rng.normal(size=200)
    Dm = np.vstack([q + rng.normal(scale=0.3, size=200), -q, rng.normal(size=200)])
    cos = SC.cosine(Dm, q[None, :])[0]
    assert cos[0] > 0.8 and cos[1] == pytest.approx(-1)
    up, dn = SC.query_sets(q, 30)
    con = SC.connectivity(Dm, up, dn)
    assert con[0] > 0.5 and con[1] < -0.5 and abs(con[2]) < 0.4


def test_auroc_basic():
    s = pd.Series([0.9, 0.8, 0.1, 0.2, 0.3], index=list("abcde"))
    assert SC.auroc(s, {"a", "b"}) == pytest.approx(1.0)
    assert SC.auroc(s, {"c"}) == pytest.approx(0.0)
    assert SC.auroc_perm_p(s, {"a", "b"}, 2000, 0) < 0.2


# ---------------------------------------------------------------- phase 4
def test_patient_folds_keep_patients_together():
    subj = pd.Series(["p1", "p1", "p2", "p3", "p3", "p4", "p5"])
    f = M.patient_folds(subj, 3, 0)
    for p in subj.unique():
        assert len(set(f[(subj == p).values])) == 1


def test_perm_test_detects_signal_and_null():
    n_s, n_d = 40, 30
    base = rng.normal(size=(n_s, n_d))
    Y = pd.DataFrame(base + rng.normal(scale=0.5, size=(n_s, n_d)))
    S_good = pd.DataFrame(base)
    S_null = pd.DataFrame(rng.normal(size=(n_s, n_d)))
    _, p_good, _ = M.perm_test(S_good, Y, 20, 500, 0)
    _, p_null, _ = M.perm_test(S_null, Y, 20, 500, 0)
    assert p_good < 0.01 and p_null > 0.01


def test_global_baseline_is_constant_per_drug_and_centred_metric_zero():
    Y = pd.DataFrame(rng.normal(size=(30, 25)), index=[f"s{i}" for i in range(30)])
    folds = np.arange(30) % 5
    P_, Ptr = M.fit_predict("global", Y, folds)
    met = M.metrics(P_, Ptr, Y, folds, 20)
    assert np.allclose(met["per_sample_centred"].dropna(), 0)


# ---------------------------------------------------------------- phase 5
class _Mapper:
    """MOCK mapper: score = mean expression, position = first two genes."""
    genes = np.array(["a", "b", "c"])

    class knn:
        @staticmethod
        def predict(Z):
            return np.c_[Z.mean(1), Z[:, 0], Z[:, 1]]

    @staticmethod
    def embed(X):
        return X


def test_simulate_zero_shift_is_static_and_outcomes():
    X0 = np.array([[0.1, 0.1, 0.1], [0.9, 0.9, 0.9]])
    traj = SIM.simulate(_Mapper, X0, np.zeros(3), np.linspace(0, 1, 5))
    assert np.allclose(traj[0], traj[-1])
    traj2 = SIM.simulate(_Mapper, X0, np.ones(3), np.linspace(0, 1, 5))
    oc = SIM.outcome(traj2, 0.5)
    assert oc["already_mature"] == 1 and oc["matured"] == 1 and oc["resistant"] == 0


def test_mutation_tags_parse():
    obs = pd.DataFrame({"patient": ["P1", "P1", "P2"], "MutTranscripts": ["NPM1.W288fs/21 FLT3.D835Y/3", "NPM1.W288fs/9", np.nan]})
    tags = SIM.mutation_tags(obs)
    assert tags == {"P1": ["FLT3 D835Y", "NPM1 W288fs"]}
