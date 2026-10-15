"""Experiment 3, F4: combination prediction (plan: reports/exp3_f4_plan.md).

The model carried forward from F2 is linear, so applying two drugs in sequence is identical to
summing their shifts. Only composition rules that differ for a linear model are testable here.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def compose(a: pd.Series, b: pd.Series, rule: str, norm: float) -> pd.Series:
    """Combine two single-drug shift vectors. All rules are renormalized to the same length."""
    idx = a.index.intersection(b.index)
    a, b = a.loc[idx], b.loc[idx]
    if rule == "additive":
        v = a + b
    elif rule == "elementwise max":
        v = pd.Series(np.where(np.abs(a.values) >= np.abs(b.values), a.values, b.values), index=idx)
    elif rule == "mean":
        v = (a + b) / 2
    else:
        raise ValueError(f"unknown rule {rule}")
    return v / max(np.linalg.norm(v.values), 1e-9) * norm


def signature_weights(sig: pd.DataFrame, genes: pd.Index) -> pd.Series:
    """Linear functional giving the maturation shift: +1/n_up on up genes, -1/n_down on down genes."""
    up = [g for g in sig.loc[sig.direction.str.startswith("up"), "gene"] if g in genes]
    dn = [g for g in sig.loc[sig.direction.str.startswith("down"), "gene"] if g in genes]
    w = pd.Series(0.0, index=genes)
    if up:
        w[up] = 1.0 / len(up)
    if dn:
        w[dn] = -1.0 / len(dn)
    return w


def all_pair_maturation_scores(V: np.ndarray, w: np.ndarray, norm: float) -> np.ndarray:
    """Maturation shift of every additive pair, computed in closed form.

    For unit-length-scaled rows, score(i, j) = (w.v_i + w.v_j) / ||v_i + v_j||, and
    ||v_i + v_j||^2 = ||v_i||^2 + ||v_j||^2 + 2 v_i.v_j, so the whole matrix follows from
    the per-drug scores and the Gram matrix.
    """
    a = V @ w
    G = V @ V.T
    d = np.diag(G)
    denom = np.sqrt(np.maximum(d[:, None] + d[None, :] + 2 * G, 1e-18))
    return (a[:, None] + a[None, :]) / denom
