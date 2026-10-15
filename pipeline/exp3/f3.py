"""Experiment 3, F3: which cell states are over-represented after treatment (plan: reports/exp3_f3_plan.md).

The target is a substitute for the brief's definition: single-cell data never sees the same cell twice,
so we ask which untreated states are enriched after treatment rather than which cells resist.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

MARKERS = ["CD34", "PROM1", "CD14", "LYZ", "CSF1R", "MPO", "ELANE"]
FEATURE_SETS = {
    "S1 maturity score only": ["score"],
    "S2 cell type only": ["__celltype__"],
    "F full state": ["score", "x", "y", "knn_dist", "__celltype__"] + MARKERS,
}


def load_10x(raw: Path, prefix: str) -> tuple[sp.csr_matrix, pd.Index, pd.Index]:
    """Read a 10x mtx triple. Returns cells x genes counts, cell barcodes, gene symbols."""
    with gzip.open(raw / f"{prefix}_matrix.mtx.gz", "rb") as f:
        M = scipy.io.mmread(f).tocsr()                      # genes x cells
    with gzip.open(raw / f"{prefix}_barcodes.tsv.gz", "rt") as f:
        cells = pd.Index([l.strip() for l in f])
    with gzip.open(raw / f"{prefix}_features.tsv.gz", "rt") as f:
        genes = pd.Index([l.split("\t")[1].strip() for l in f])
    return M.T.tocsr(), cells, genes


def collapse_duplicate_genes(X: sp.csr_matrix, genes: pd.Index) -> tuple[sp.csr_matrix, pd.Index]:
    if not genes.duplicated().any():
        return X, genes
    codes, uniq = pd.factorize(genes)
    agg = sp.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)),
                        shape=(len(codes), len(uniq)), dtype=np.float32)
    return (X @ agg).tocsr(), pd.Index(uniq)


def place(mapper, X: sp.csr_matrix, genes: pd.Index, ref_Z: np.ndarray) -> pd.DataFrame:
    """Map cells onto the Phase 1 landscape, with the Phase 2 off-manifold distance."""
    from pipeline.phase1_landscape.mapping import lognorm
    A, _ = mapper.align(lognorm(X), genes)
    Z = mapper.embed(A)
    out = pd.DataFrame(mapper.knn.predict(Z), columns=list(mapper.targets))
    d, _ = mapper.knn.kneighbors(Z)
    out["knn_dist"] = d.mean(1)
    return out


def marker_frame(X: sp.csr_matrix, genes: pd.Index) -> pd.DataFrame:
    from pipeline.phase1_landscape.mapping import lognorm
    L = lognorm(X)
    cols = {}
    for g in MARKERS:
        i = genes.get_indexer([g])[0]
        cols[g] = np.asarray(L[:, i].todense()).ravel() if i >= 0 else np.zeros(X.shape[0])
    return pd.DataFrame(cols)


def design(df: pd.DataFrame, feats: list[str]) -> np.ndarray:
    """Build a model matrix; '__celltype__' expands to one-hot columns."""
    parts = []
    for f in feats:
        if f == "__celltype__":
            parts.append(pd.get_dummies(df["cell_type"].astype(str), prefix="ct").astype(float))
        else:
            parts.append(df[[f]].astype(float))
    return pd.concat(parts, axis=1).values


def fit_eval(df: pd.DataFrame, feats: list[str], folds: np.ndarray, seed: int = 0) -> dict:
    """AUROC for separating treated from untreated cells, on held-out folds."""
    X, y = design(df, feats), df["treated"].values.astype(int)
    pred = np.full(len(y), np.nan)
    for f in np.unique(folds):
        te, tr = folds == f, folds != f
        if len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2:
            continue
        m = HistGradientBoostingClassifier(max_iter=200, random_state=seed).fit(X[tr], y[tr])
        pred[te] = m.predict_proba(X[te])[:, 1]
    ok = ~np.isnan(pred)
    per_fold = []
    for f in np.unique(folds):
        s = ok & (folds == f)
        if s.sum() > 10 and len(np.unique(y[s])) == 2:
            per_fold.append(roc_auc_score(y[s], pred[s]))
    return {"auroc": float(roc_auc_score(y[ok], pred[ok])) if ok.sum() and len(np.unique(y[ok])) == 2 else float("nan"),
            "auroc_per_fold_mean": float(np.mean(per_fold)) if per_fold else float("nan"),
            "n_folds_scored": len(per_fold), "n_cells": int(ok.sum()), "pred": pred}


def maturity_direction(df: pd.DataFrame, threshold: float) -> dict:
    """Direction of the maturity shift: the readout least vulnerable to the batch confound."""
    a, b = df[~df.treated.astype(bool)], df[df.treated.astype(bool)]
    return {"n_untreated": len(a), "n_treated": len(b),
            "mean_score_untreated": float(a.score.mean()), "mean_score_treated": float(b.score.mean()),
            "delta_mean_score": float(b.score.mean() - a.score.mean()),
            "stuck_frac_untreated": float((a.score < threshold).mean()),
            "stuck_frac_treated": float((b.score < threshold).mean()),
            "delta_stuck_frac": float((b.score < threshold).mean() - (a.score < threshold).mean())}


def permutation_importance_auroc(df: pd.DataFrame, feats: list[str], folds: np.ndarray,
                                 n_repeats: int = 5, seed: int = 0) -> pd.Series:
    """Drop in held-out AUROC when each feature is shuffled."""
    base = fit_eval(df, feats, folds, seed)["auroc"]
    rng = np.random.default_rng(seed)
    out = {}
    for f in feats:
        drops = []
        for _ in range(n_repeats):
            d = df.copy()
            col = "cell_type" if f == "__celltype__" else f
            d[col] = d[col].values[rng.permutation(len(d))]
            drops.append(base - fit_eval(d, feats, folds, seed)["auroc"])
        out["cell type" if f == "__celltype__" else f] = float(np.mean(drops))
    return pd.Series(out).sort_values(ascending=False)
