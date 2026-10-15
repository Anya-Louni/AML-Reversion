"""Experiment 2, Step 1: how much patient-specific signal does Beat AML ex vivo response contain?

Within one patient, ranking drugs is unaffected by that patient's overall sensitivity, so the only
signal a patient-specific method can use is the drug x patient interaction: the residual of
    AUC = grand mean + drug effect + patient effect + residual.
The ceiling compares the variance of that residual with the variance of measurement noise.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from pipeline.exp2.curves import curve_auc


def additive_fit(Y: pd.DataFrame, n_iter: int = 100, tol: float = 1e-8) -> dict:
    """Alternating least squares for Y = mu + row + col + e on observed entries only."""
    M = Y.values.astype(float)
    obs = ~np.isnan(M)
    mu = np.nanmean(M)
    r = np.zeros(M.shape[0]); c = np.zeros(M.shape[1])
    def _mean(X, axis):  # mean over observed entries; rows/columns with no data get effect 0
        cnt = obs.sum(axis)
        return np.where(cnt > 0, np.nansum(np.where(obs, X, 0.0), axis) / np.maximum(cnt, 1), 0.0)

    for _ in range(n_iter):
        c_new = _mean(M - mu - r[:, None], 0)
        r_new = _mean(M - mu - c_new[None, :], 1)
        r_new -= r_new.mean()
        done = max(np.abs(c_new - c).max(), np.abs(r_new - r).max()) < tol
        r, c = r_new, c_new
        if done:
            break
    fitted = mu + r[:, None] + c[None, :]
    resid = np.where(obs, M - fitted, np.nan)
    return {"mu": mu, "row": pd.Series(r, index=Y.index), "col": pd.Series(c, index=Y.columns),
            "resid": pd.DataFrame(resid, index=Y.index, columns=Y.columns)}


def variance_shares(Y: pd.DataFrame, fit: dict) -> dict:
    M = Y.values; obs = ~np.isnan(M)
    rows = np.broadcast_to(fit["row"].values[:, None], M.shape)[obs]
    cols = np.broadcast_to(fit["col"].values[None, :], M.shape)[obs]
    res = fit["resid"].values[obs]
    tot = np.var(M[obs])
    return {"total_var": float(tot), "drug": float(np.var(cols) / tot), "patient": float(np.var(rows) / tot),
            "residual": float(np.var(res) / tot), "cross_terms": float(1 - (np.var(cols) + np.var(rows) + np.var(res)) / tot)}


def replicate_aucs(raw: pd.DataFrame) -> pd.DataFrame:
    """AUC of every separate technical replicate curve, for sample-drug pairs measured more than once.

    A replicate curve = one (replicate, plate, run) group with at least 5 concentrations.
    """
    raw = raw[raw.passed_qc].copy()
    raw["dbgap_rnaseq_sample"] = raw.dbgap_rnaseq_sample.fillna("NA")
    key = ["dbgap_subject_id", "dbgap_rnaseq_sample", "inhibitor"]
    rep = ["replicate", "plate_num", "run_index"]
    n_curves = raw.groupby(key)[rep].apply(lambda d: len(d.drop_duplicates()))
    multi = n_curves[n_curves >= 2].index
    sub = raw.set_index(key).loc[multi].reset_index()
    out = []
    for k, d in sub.groupby(key + rep):
        if d.well_concentration.nunique() >= 5:
            out.append((*k, curve_auc(d.well_concentration.values, d.normalized_viability.values)))
    return pd.DataFrame(out, columns=key + rep + ["auc"])


def n_curves_per_pair(raw: pd.DataFrame) -> pd.Series:
    raw = raw[raw.passed_qc].copy()
    raw["dbgap_rnaseq_sample"] = raw.dbgap_rnaseq_sample.fillna("NA")
    key = ["dbgap_rnaseq_sample", "inhibitor"]
    return raw.groupby(key)[["replicate", "plate_num", "run_index"]].apply(lambda d: len(d.drop_duplicates()))


def per_row_spearman(A: pd.DataFrame, B: pd.DataFrame, min_n: int) -> pd.Series:
    out = {}
    for i in A.index.intersection(B.index):
        a, b = A.loc[i], B.loc[i]
        ok = a.notna() & b.notna()
        if ok.sum() >= min_n and a[ok].std() > 0 and b[ok].std() > 0:
            out[i] = np.corrcoef(rankdata(a[ok]), rankdata(b[ok]))[0, 1]
    return pd.Series(out, dtype=float)
