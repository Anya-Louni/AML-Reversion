"""Place van Galen malignant cells on the Phase 1 landscape."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from pipeline.phase1_landscape.mapping import LandscapeMapper, lognorm


def maturity_threshold(score: np.ndarray, cell_type: pd.Series, immature, mature) -> dict:
    """Score cut that best separates immature from mature normal cells (max balanced accuracy)."""
    imm = cell_type.isin(immature).values
    mat = cell_type.isin(mature).values
    grid = np.linspace(0, 1, 1001)
    ba = np.array([(np.mean(score[imm] < t) + np.mean(score[mat] >= t)) / 2 for t in grid])
    t = float(grid[np.argmax(ba)])
    return {"threshold": t, "balanced_accuracy": float(ba.max()),
            "sensitivity_immature": float(np.mean(score[imm] < t)), "specificity_mature": float(np.mean(score[mat] >= t))}


def map_cells(mapper: LandscapeMapper, counts, genes, ref_Z: np.ndarray) -> tuple[pd.DataFrame, dict]:
    """Map raw-count cells; also return each cell's mean distance to its k reference neighbours.

    Off-manifold reference: the same distance for every reference cell computed leave-one-out.
    """
    X, frac = mapper.align(lognorm(counts), genes)
    Z = mapper.embed(X)
    pred = pd.DataFrame(mapper.knn.predict(Z), columns=list(mapper.targets))
    d, _ = mapper.knn.kneighbors(Z)
    pred["knn_dist"] = d.mean(1)
    k = mapper.knn.n_neighbors
    dref, _ = mapper.knn.kneighbors(ref_Z, n_neighbors=k + 1)
    ref_dist = dref[:, 1:].mean(1)
    return pred, {"gene_fraction_present": frac, "ref_dist": ref_dist}


def like_order_agreement(score, like_type: pd.Series, order) -> float:
    r = like_type.map({c: i for i, c in enumerate(order)})
    ok = r.notna().values
    return float(spearmanr(score[ok], r.values[ok])[0])


def patient_summary(obs: pd.DataFrame, threshold: float, min_cells: int) -> pd.DataFrame:
    g = obs.groupby("patient", observed=True)
    out = g.agg(n_malignant=("score", "size"), median_score=("score", "median"),
                mean_score=("score", "mean"), q25=("score", lambda s: s.quantile(0.25)),
                q75=("score", lambda s: s.quantile(0.75)),
                stuck_fraction=("score", lambda s: float(np.mean(s < threshold))),
                offmanifold_fraction=("off_manifold", "mean"))
    comp = obs.groupby(["patient", "CellType"], observed=True).size().unstack(fill_value=0)
    comp = comp.div(comp.sum(1), axis=0).add_prefix("frac_")
    out = out.join(comp)
    out["low_cell_count"] = out.n_malignant < min_cells
    return out.sort_values("stuck_fraction", ascending=False)
