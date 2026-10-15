"""Experiment 3, F2: models that try to beat the linear shift (plan: reports/exp3_f2_plan.md)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge

MYELOID = ["HL60", "THP1", "U937", "SKM1", "NOMO1"]


def per_cell_line_signatures(sig_obs: pd.DataFrame, X: np.ndarray, genes, line: str) -> pd.DataFrame:
    """Mean signature per drug within one cell line."""
    m = (sig_obs.cell_iname == line).values
    if m.sum() == 0:
        return pd.DataFrame(columns=genes)
    df = pd.DataFrame(X[m], columns=genes)
    df["drug"] = sig_obs.drug.values[m]
    return df.groupby("drug").mean()


def consensus_excluding(sig_obs: pd.DataFrame, X: np.ndarray, genes, exclude_line: str,
                        lines=MYELOID) -> pd.DataFrame:
    """Mean signature per drug over myeloid lines other than `exclude_line`."""
    keep = sig_obs.cell_iname.isin([c for c in lines if c != exclude_line]).values
    df = pd.DataFrame(X[keep], columns=genes)
    df["drug"] = sig_obs.drug.values[keep]
    return df.groupby("drug").mean()


class ReducedRankMap:
    """Low-rank ridge map from one signature space to another, for small n.

    Both sides are projected onto `rank` principal components fitted on the training drugs; a ridge
    model maps the source components to the target components, and the prediction is reconstructed.
    """

    def __init__(self, rank: int = 3, alpha: float = 1.0):
        self.rank, self.alpha = rank, alpha

    def fit(self, Xs: np.ndarray, Ys: np.ndarray) -> "ReducedRankMap":
        k = min(self.rank, Xs.shape[0] - 1, Xs.shape[1])
        self.px = PCA(k, random_state=0).fit(Xs)
        self.py = PCA(k, random_state=0).fit(Ys)
        self.ridge = Ridge(alpha=self.alpha).fit(self.px.transform(Xs), self.py.transform(Ys))
        return self

    def predict(self, Xs: np.ndarray) -> np.ndarray:
        z = self.ridge.predict(self.px.transform(Xs))
        return self.py.inverse_transform(np.asarray(z).reshape(len(Xs), -1))


def select_map(Xs: pd.DataFrame, Ys: pd.DataFrame, ranks, alphas, exclude: list[str]) -> tuple[dict, pd.DataFrame]:
    """Leave-one-drug-out over the training drugs to pick rank and penalty.

    `exclude` drugs are removed before any fitting or selection, so they are never seen.
    """
    drugs = [d for d in Xs.index if d not in exclude]
    Xs, Ys = Xs.loc[drugs], Ys.loc[drugs]
    rows = []
    for r in ranks:
        for a in alphas:
            cos = []
            for i, d in enumerate(drugs):
                tr = [x for x in drugs if x != d]
                m = ReducedRankMap(r, a).fit(Xs.loc[tr].values, Ys.loc[tr].values)
                p = m.predict(Xs.loc[[d]].values)[0]
                y = Ys.loc[d].values
                cos.append(float(p @ y / max(np.linalg.norm(p) * np.linalg.norm(y), 1e-12)))
            rows.append({"rank": r, "alpha": a, "mean_cosine": float(np.mean(cos))})
    grid = pd.DataFrame(rows)
    best = grid.loc[grid.mean_cosine.idxmax()]
    return {"rank": int(best["rank"]), "alpha": float(best["alpha"]), "cv_cosine": float(best.mean_cosine)}, grid


def baseline_consensus_cosine(Xs: pd.DataFrame, Ys: pd.DataFrame, exclude: list[str]) -> float:
    """How well the raw consensus already matches the cell-line signature, for reference."""
    drugs = [d for d in Xs.index if d not in exclude]
    c = []
    for d in drugs:
        x, y = Xs.loc[d].values, Ys.loc[d].values
        c.append(float(x @ y / max(np.linalg.norm(x) * np.linalg.norm(y), 1e-12)))
    return float(np.mean(c))


def time_interpolate(changes: dict[float, pd.Series], target_t: float) -> pd.Series:
    """Linear interpolation in time between the two nearest fitted timepoints."""
    ts = sorted(changes)
    lo = max([t for t in ts if t < target_t], default=None)
    hi = min([t for t in ts if t > target_t], default=None)
    if lo is None or hi is None:
        raise ValueError("target timepoint is not bracketed by the fitted timepoints")
    w = (target_t - lo) / (hi - lo)
    a, b = changes[lo], changes[hi]
    idx = a.index.intersection(b.index)
    return a.loc[idx] * (1 - w) + b.loc[idx] * w
