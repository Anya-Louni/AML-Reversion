"""Reversion scores: how strongly a drug signature points along a patient's push-to-maturity vector.

Two methods:
  cosine       : cosine similarity between the drug signature and the push vector (= -stuck), on shared genes.
  connectivity : rank-based KS connectivity score (Lamb et al. 2006, CMap). The query is the top-N genes
                 the patient needs raised (up set) and lowered (down set); the score is positive when the
                 drug raises the up set and lowers the down set.
Higher = pushes toward maturity, for both.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def cosine(D: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """D: drugs x genes, Q: queries x genes. Returns queries x drugs."""
    Dn = D / np.maximum(np.linalg.norm(D, axis=1, keepdims=True), 1e-12)
    Qn = Q / np.maximum(np.linalg.norm(Q, axis=1, keepdims=True), 1e-12)
    return Qn @ Dn.T


def _ks(rank_pos: np.ndarray, n: int) -> np.ndarray:
    """rank_pos: drugs x t, 1-based positions of the query genes in each drug's descending ranking."""
    V = np.sort(rank_pos, axis=1)
    t = V.shape[1]
    j = np.arange(1, t + 1)
    a = (j / t - V / n).max(1)
    b = (V / n - (j - 1) / t).max(1)
    return np.where(a > b, a, -b)


def connectivity(D: np.ndarray, up_idx: np.ndarray, down_idx: np.ndarray) -> np.ndarray:
    """D: drugs x genes. up_idx / down_idx: gene column indices of the query sets. Returns one score per drug."""
    n = D.shape[1]
    order = np.argsort(-D, axis=1)
    rank = np.empty_like(order)
    rank[np.arange(D.shape[0])[:, None], order] = np.arange(1, n + 1)
    ks_up, ks_dn = _ks(rank[:, up_idx], n), _ks(rank[:, down_idx], n)
    return np.where(np.sign(ks_up) == np.sign(ks_dn), 0.0, (ks_up - ks_dn) / 2)


def query_sets(push: np.ndarray, size: int) -> tuple[np.ndarray, np.ndarray]:
    o = np.argsort(push)
    return o[::-1][:size], o[:size]


def score_patients(drug_sig: pd.DataFrame, push: pd.DataFrame, size: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """push: patients x genes. Returns (cosine, connectivity) as patients x drugs, on shared genes."""
    genes = drug_sig.columns.intersection(push.columns)
    D, Q = drug_sig[genes].values, push[genes].values
    cos = pd.DataFrame(cosine(D, Q), index=push.index, columns=drug_sig.index)
    con = {}
    for i, p in enumerate(push.index):
        up, dn = query_sets(Q[i], min(size, len(genes) // 2))
        con[p] = connectivity(D, up, dn)
    return cos, pd.DataFrame(con, index=drug_sig.index).T


def auroc(scores: pd.Series, positives: set) -> float:
    y = scores.index.isin(list(positives))
    if y.sum() == 0 or (~y).sum() == 0:
        return np.nan
    r = scores.rank().values
    return float((r[y].sum() - y.sum() * (y.sum() + 1) / 2) / (y.sum() * (~y).sum()))


def auroc_perm_p(scores: pd.Series, positives: set, n_perm: int, seed: int) -> float:
    obs = auroc(scores, positives)
    k = int(scores.index.isin(list(positives)).sum())
    rng = np.random.default_rng(seed)
    vals = scores.values
    r = pd.Series(vals).rank().values
    n = len(vals)
    null = np.empty(n_perm)
    for i in range(n_perm):
        pos = rng.choice(n, k, replace=False)
        null[i] = (r[pos].sum() - k * (k + 1) / 2) / (k * (n - k))
    return float((np.sum(null >= obs) + 1) / (n_perm + 1))
