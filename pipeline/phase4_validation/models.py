"""Models and metrics for Phase 4 (see reports/phase4_plan.md)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.linear_model import LinearRegression, RidgeCV

ALPHAS = np.logspace(-2, 4, 13)


def patient_folds(subjects: pd.Series, k: int, seed: int) -> np.ndarray:
    u = np.array(sorted(subjects.unique()))
    rng = np.random.default_rng(seed)
    rng.shuffle(u)
    fold_of = {s: i % k for i, s in enumerate(u)}
    return subjects.map(fold_of).values


def _std(Xtr, Xte):
    mu, sd = Xtr.mean(0), Xtr.std(0)
    sd[sd == 0] = 1
    return (Xtr - mu) / sd, (Xte - mu) / sd


def fit_predict(kind: str, Y: pd.DataFrame, folds: np.ndarray, R: pd.DataFrame | None = None,
                F: pd.DataFrame | None = None, seed: int = 0) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (test predictions, in-sample training predictions per fold) as samples x drugs.

    Training predictions are needed to centre scores with training-fold statistics (metric 2).
    `train_pred` is stored as a dict fold -> DataFrame (train samples x drugs).
    """
    P = pd.DataFrame(np.nan, index=Y.index, columns=Y.columns)
    train = {}
    for f in np.unique(folds):
        te, tr = folds == f, folds != f
        Ptr = pd.DataFrame(np.nan, index=Y.index[tr], columns=Y.columns)
        for d in Y.columns:
            y = Y[d].values
            ok = tr & ~np.isnan(y)
            if kind == "global":
                m = np.nanmean(y[ok]); P.loc[te, d] = m; Ptr[d] = m
            elif kind == "raw":
                P.loc[te, d] = R[d].values[te]; Ptr[d] = R[d].values[tr]
            elif kind == "random":
                rng = np.random.default_rng(seed + f)
                P.loc[te, d] = rng.normal(size=te.sum()); Ptr[d] = rng.normal(size=tr.sum())
            else:
                if kind == "calibrated":
                    X = R[[d]].values
                    model = LinearRegression()
                else:
                    X = F.values if kind in ("mutation", "fractions") else np.c_[R[[d]].values, F.values]
                    model = RidgeCV(alphas=ALPHAS)
                if ok.sum() < 10:
                    continue
                Xtr, Xte = _std(X[ok], X[te])
                _, Xall_tr = _std(X[ok], X[tr])
                model.fit(Xtr, y[ok])
                P.loc[te, d] = model.predict(Xte); Ptr[d] = model.predict(Xall_tr)
        train[f] = Ptr
    return P, train


def _spearman_rows(A: np.ndarray, B: np.ndarray, min_n: int) -> np.ndarray:
    out = np.full(A.shape[0], np.nan)
    for i in range(A.shape[0]):
        ok = ~np.isnan(A[i]) & ~np.isnan(B[i])
        if ok.sum() < min_n:
            continue
        a, b = A[i, ok], B[i, ok]
        if a.std() == 0 or b.std() == 0:
            out[i] = 0.0
            continue
        out[i] = np.corrcoef(rankdata(a), rankdata(b))[0, 1]
    return out


def centre(M: pd.DataFrame, Y: pd.DataFrame, folds: np.ndarray, train_ref: dict | None) -> pd.DataFrame:
    """Centre each drug column of M on test samples using the training-fold mean/SD of `train_ref` (or of Y)."""
    out = pd.DataFrame(np.nan, index=M.index, columns=M.columns)
    for f in np.unique(folds):
        te = folds == f
        ref = train_ref[f] if train_ref is not None else Y[folds != f]
        mu, sd = ref.mean(), ref.std()
        z = (M[te] - mu) / sd.replace(0, np.nan)
        out[te] = z.fillna(0.0).where(M[te].notna())
    return out


def metrics(P: pd.DataFrame, Ptrain: dict, Y: pd.DataFrame, folds: np.ndarray, min_drugs: int) -> dict:
    m1 = _spearman_rows(P.values, Y.values, min_drugs)
    Pc, Yc = centre(P, Y, folds, Ptrain), centre(Y, Y, folds, None)
    m2 = _spearman_rows(Pc.values, Yc.values, min_drugs)
    m3 = _spearman_rows(P.values.T, Y.values.T, 20)
    return {"per_sample_raw": pd.Series(m1, index=Y.index), "per_sample_centred": pd.Series(m2, index=Y.index),
            "per_drug": pd.Series(m3, index=Y.columns), "Pc": Pc, "Yc": Yc}


def perm_test(S: pd.DataFrame, Y: pd.DataFrame, min_drugs: int, n_perm: int, seed: int) -> tuple[float, float, np.ndarray]:
    """Within-sample permutation of drug labels. Statistic: mean per-sample Spearman. One-sided p (greater)."""
    rng = np.random.default_rng(seed)
    obs, null_sum, n = [], np.zeros(n_perm), 0
    for i in range(len(Y)):
        ok = ~np.isnan(S.values[i]) & ~np.isnan(Y.values[i])
        if ok.sum() < min_drugs:
            continue
        a, b = rankdata(S.values[i, ok]), rankdata(Y.values[i, ok])
        a = (a - a.mean()) / a.std(); b = (b - b.mean()) / b.std()
        obs.append((a * b).mean())
        perms = np.argsort(rng.random((n_perm, len(a))), axis=1)
        null_sum += (a[perms] * b).mean(1)
        n += 1
    stat = float(np.mean(obs)); null = null_sum / n
    return stat, float((np.sum(null >= stat) + 1) / (n_perm + 1)), null


def boot_ci(values: pd.Series, groups: pd.Series, n_boot: int, seed: int) -> tuple[float, float, float]:
    v = values.dropna()
    g = groups.loc[v.index]
    by = v.groupby(g).apply(list)
    keys = by.index.values
    rng = np.random.default_rng(seed)
    stats = []
    for _ in range(n_boot):
        pick = rng.choice(len(keys), len(keys), replace=True)
        stats.append(np.mean([x for k in pick for x in by.iloc[k]]))
    return float(v.mean()), float(np.quantile(stats, 0.025)), float(np.quantile(stats, 0.975))


def boot_ci_simple(values: pd.Series, n_boot: int, seed: int) -> tuple[float, float, float]:
    v = values.dropna().values
    rng = np.random.default_rng(seed)
    s = [rng.choice(v, len(v), replace=True).mean() for _ in range(n_boot)]
    return float(v.mean()), float(np.quantile(s, 0.025)), float(np.quantile(s, 0.975))
