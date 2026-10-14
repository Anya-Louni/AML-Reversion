"""Map any cell (or shifted expression vector) onto the Phase 1 landscape.

The map is a fixed, non-transductive function:
  log-normalized expression (HVGs) -> centre with reference means -> PCA projection
  -> distance-weighted kNN regression onto (differentiation score, landscape x, y).
Being fixed matters: Phase 5 maps drug-shifted expression frame by frame, so the
projection must not be refitted or batch-corrected per query.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsRegressor


def lognorm(X, target_sum: float = 1e4):
    """Counts -> log1p(counts / total * target_sum), same as scanpy normalize_total + log1p."""
    X = sp.csr_matrix(X, dtype=np.float32)
    tot = np.asarray(X.sum(1)).ravel()
    tot[tot == 0] = 1
    X = sp.diags(target_sum / tot) @ X
    X.data = np.log1p(X.data)
    return X.tocsr()


@dataclass
class LandscapeMapper:
    genes: np.ndarray          # HVGs, in model order
    mean: np.ndarray           # reference mean of log-normalized HVG expression
    pca: PCA
    knn: KNeighborsRegressor
    targets: tuple = ("score", "x", "y")

    @classmethod
    def fit(cls, Xlog_hvg: np.ndarray, genes, pca: PCA, mean, Y: np.ndarray, k: int) -> "LandscapeMapper":
        Z = pca.transform(Xlog_hvg - mean)
        knn = KNeighborsRegressor(n_neighbors=k, weights="distance").fit(Z, Y)
        return cls(np.asarray(genes), np.asarray(mean), pca, knn)

    def align(self, X, genes) -> tuple[np.ndarray, float]:
        """Reorder a (cells x genes) log-normalized matrix to model genes.

        Missing genes are filled with the reference mean, so they contribute nothing
        after centring. Returns the dense matrix and the fraction of model genes present.
        """
        idx = pd.Index(genes).get_indexer(self.genes)
        present = idx >= 0
        X = X.toarray() if sp.issparse(X) else np.asarray(X)
        out = np.tile(self.mean, (X.shape[0], 1)).astype(np.float32)
        out[:, present] = X[:, idx[present]]
        return out, float(present.mean())

    def embed(self, Xlog_hvg: np.ndarray) -> np.ndarray:
        return self.pca.transform(Xlog_hvg - self.mean)

    def predict(self, X, genes=None) -> pd.DataFrame:
        if genes is not None:
            X, _ = self.align(X, genes)
        P = self.knn.predict(self.embed(np.asarray(X)))
        return pd.DataFrame(P, columns=list(self.targets))


def evaluate(Xlog_hvg, Y, groups, pca, mean, genes, k, test_fraction, seed, strata) -> dict:
    """Held-out mapping error: (1) random split stratified by cell type, (2) leave-one-sample-out."""
    from scipy.stats import spearmanr
    from sklearn.model_selection import train_test_split

    rng_idx = np.arange(len(Y))
    tr, te = train_test_split(rng_idx, test_size=test_fraction, random_state=seed, stratify=strata)
    m = LandscapeMapper.fit(Xlog_hvg[tr], genes, pca, mean, Y[tr], k)
    P = m.knn.predict(m.embed(Xlog_hvg[te]))
    span = np.ptp(Y[:, 1:], axis=0).mean()
    random_split = {
        "n_test": len(te),
        "score_mae": float(np.mean(np.abs(P[:, 0] - Y[te, 0]))),
        "score_spearman": float(spearmanr(P[:, 0], Y[te, 0])[0]),
        "xy_rmse_frac_of_range": float(np.sqrt(np.mean(np.sum((P[:, 1:] - Y[te, 1:]) ** 2, 1))) / span),
    }
    rows = []
    for g in np.unique(groups):
        te = groups == g
        if te.sum() < 30:
            continue
        m = LandscapeMapper.fit(Xlog_hvg[~te], genes, pca, mean, Y[~te], k)
        P = m.knn.predict(m.embed(Xlog_hvg[te]))
        rows.append({"held_out_sample": g, "n": int(te.sum()),
                     "score_mae": float(np.mean(np.abs(P[:, 0] - Y[te, 0]))),
                     "score_spearman": float(spearmanr(P[:, 0], Y[te, 0])[0]),
                     "score_sd_true": float(Y[te, 0].std()),
                     "xy_rmse_frac_of_range": float(np.sqrt(np.mean(np.sum((P[:, 1:] - Y[te, 1:]) ** 2, 1))) / span)})
    loso = pd.DataFrame(rows)
    return {"random_split": random_split, "loso": loso,
            "loso_weighted_mae": float(np.average(loso.score_mae, weights=loso.n)),
            "loso_weighted_xy": float(np.average(loso.xy_rmse_frac_of_range, weights=loso.n))}
