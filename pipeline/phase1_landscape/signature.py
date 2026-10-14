"""Gene-level maturation signature: correlation of each gene with the differentiation score."""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import rankdata


def spearman_with(X, y: np.ndarray, genes, min_frac: float = 0.01, chunk: int = 2000) -> pd.DataFrame:
    """Spearman correlation of every gene (columns of X) with y, over genes detected in >= min_frac of cells."""
    X = sp.csc_matrix(X)
    det = np.asarray((X > 0).mean(0)).ravel()
    keep = np.where(det >= min_frac)[0]
    ry = rankdata(y); ry = (ry - ry.mean()) / ry.std()
    rho = np.empty(len(keep))
    for s in range(0, len(keep), chunk):
        cols = keep[s:s + chunk]
        R = np.apply_along_axis(rankdata, 0, X[:, cols].toarray())
        R = (R - R.mean(0)) / np.where(R.std(0) > 0, R.std(0), 1)
        rho[s:s + chunk] = (R * ry[:, None]).mean(0)
    return pd.DataFrame({"gene": np.asarray(genes)[keep], "spearman": rho, "frac_detected": det[keep]}) \
        .sort_values("spearman", ascending=False).reset_index(drop=True)


def define_signature(corr: pd.DataFrame, top_n: int) -> pd.DataFrame:
    up = corr.nlargest(top_n, "spearman").assign(direction="up_with_maturation")
    down = corr.nsmallest(top_n, "spearman").assign(direction="down_with_maturation")
    return pd.concat([up, down], ignore_index=True)
