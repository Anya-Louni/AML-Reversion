"""Experiment 2 unit tests: curve fits, ceiling decomposition, split logic, test lock, models and metrics.

All arrays and tables in this file are MOCK DATA generated for testing only.
They are never used by the real pipeline.
"""
import numpy as np
import pandas as pd
import pytest
from scipy.special import ndtr

from pipeline.exp2 import ceiling as C
from pipeline.exp2 import curves
from pipeline.exp2 import models as M
from pipeline.exp2 import split as SP
from pipeline.exp2.export import hierarchy_predictions

rng = np.random.default_rng(0)
CONC = np.array([0.0137, 0.0412, 0.123, 0.37, 1.11, 3.33, 10.0])


# ---------------------------------------------------------------- curves
def test_auc_of_flat_full_viability_is_max():
    assert curves.auc(20.0, 0.0, CONC.min(), CONC.max()) == pytest.approx(100 * np.log10(10 / 0.0137), rel=1e-3)


def test_probit_fit_recovers_parameters():
    a, b = 0.3, -1.2
    via = 100 * ndtr(a + b * np.log10(CONC))
    fa, fb = curves.fit_probit(CONC, via)
    assert fa == pytest.approx(a, abs=0.05) and fb == pytest.approx(b, abs=0.05)


# ---------------------------------------------------------------- ceiling
def test_additive_fit_recovers_effects_with_missing():
    rows, cols = rng.normal(0, 10, 40), rng.normal(0, 20, 15)
    Y = pd.DataFrame(100 + rows[:, None] + cols[None, :])
    Y = Y.mask(rng.random(Y.shape) < 0.2)
    f = C.additive_fit(Y)
    assert np.nanmax(np.abs(f["resid"].values)) < 1e-6
    sh = C.variance_shares(Y, f)
    assert sh["residual"] == pytest.approx(0, abs=1e-6)


def test_additive_fit_handles_empty_row():
    Y = pd.DataFrame([[1.0, 2.0], [np.nan, np.nan], [3.0, 5.0]])
    f = C.additive_fit(Y)
    assert np.isfinite(f["row"]).all() and np.isfinite(f["col"]).all()


# ---------------------------------------------------------------- split and lock
def test_split_keeps_patients_together_and_uses_wave_rule():
    samples = pd.Index(["s1", "s2", "s3", "s4", "s5"])
    subj = pd.Series(["p1", "p1", "p2", "p3", "p3"], index=samples)
    cohort = pd.Series(["Waves3+4", "Waves3+4", "Waves1+2", "Waves3+4", "Both"], index=samples)
    df = SP.make_split(samples, subj, cohort).set_index("sample")
    assert df.loc[["s1", "s2"], "split"].eq("test").all()          # all samples in waves 3+4 -> test
    assert df.loc["s3", "split"] == "development"
    assert df.loc[["s4", "s5"], "split"].eq("development").all()   # mixed labels -> development
    SP.check_split(df.reset_index())


def test_check_split_rejects_patient_on_both_sides():
    bad = pd.DataFrame({"subject": ["p1", "p1"], "split": ["test", "development"]})
    with pytest.raises(ValueError):
        SP.check_split(bad)


def test_test_lock(tmp_path):
    SP.assert_test_not_evaluated(tmp_path)
    SP.mark_test_evaluated(tmp_path)
    with pytest.raises(RuntimeError):
        SP.assert_test_not_evaluated(tmp_path)


# ---------------------------------------------------------------- models
def _toy():
    n, k, d = 120, 4, 12
    F = pd.DataFrame(rng.dirichlet(np.ones(k), n), columns=[f"t{i}" for i in range(k)], index=[f"s{i}" for i in range(n)])
    S_true = rng.normal(0, 30, (d, k))
    Y = pd.DataFrame(F.values @ S_true.T + rng.normal(0, 3, (n, d)), index=F.index, columns=[f"d{j}" for j in range(d)])
    return F, Y


def test_pure_type_sensitivity_reproduces_predictions():
    F, Y = _toy()
    means = Y.mean()
    R = Y - means
    m = M.PerDrugRidge(1.0).fit(R, F, F.index)
    S = m.pure_type_sensitivity(F, means)
    pred = m.predict(F, F.index[:10]) + means
    for s in F.index[:10]:
        direct = hierarchy_predictions(F.loc[s], S, means)["pred_response"]
        np.testing.assert_allclose(direct.values, pred.loc[s, S.index].values, atol=1e-6)


def test_hierarchy_model_generalizes_on_signal():
    F, Y = _toy()
    means = Y.iloc[:80].mean()
    R = Y - means
    m = M.PerDrugRidge(1.0).fit(R, F, F.index[:80])
    P = m.predict(F, F.index[80:])
    assert M.pooled_pearson(P, R.loc[F.index[80:]]) > 0.5


def test_metrics_basic():
    R = pd.DataFrame(rng.normal(size=(6, 25)), index=list("abcdef"))
    subj = pd.Series(["p1", "p1", "p2", "p3", "p4", "p5"], index=R.index)
    ps = M.per_sample_spearman(R, R)
    assert np.allclose(ps, 1)
    pp = M.per_patient(ps, subj)
    assert len(pp) == 5
    assert M.top_k_hits(R, R).eq(1.0).all()
    mean, lo, hi = M.bootstrap_ci(pd.Series(rng.normal(0.3, 0.1, 200)), 2000)
    assert lo < mean < hi and lo > 0.2
    assert M.sign_flip_test(pd.Series(rng.normal(0.3, 0.1, 50)), 2000) < 0.01
    assert M.sign_flip_test(pd.Series(rng.normal(0.0, 0.1, 50)), 2000, alternative="two-sided") > 0.01


def test_patient_folds_deterministic_and_grouped():
    subj = pd.Series(["p1", "p1", "p2", "p3", "p4", "p4", "p5"], index=list("abcdefg"))
    f1, f2 = M.patient_folds(subj, 3, 0), M.patient_folds(subj, 3, 0)
    assert f1.equals(f2)
    assert f1["a"] == f1["b"] and f1["e"] == f1["f"]
