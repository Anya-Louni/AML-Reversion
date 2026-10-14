"""Normal myeloid differentiation landscape from van Galen normal cells."""
from __future__ import annotations

import anndata as ad
import harmonypy
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.ndimage import gaussian_filter
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

from pipeline.common import get_logger
from pipeline.phase1_landscape.mapping import lognorm

log = get_logger("phase1")
sc.settings.verbosity = 0


def select_reference(vg: ad.AnnData, p1: dict) -> tuple[ad.AnnData, dict]:
    o = vg.obs
    keep = (o.PredictionRefined == "normal") & o.CellType.isin(p1["cell_types"]) & (o.sample_type != "cell_line")
    if not p1["include_aml_normal_cells"]:
        keep &= o.sample_type == "healthy_BM"
    a = vg[keep.values].copy()
    a.obs["source"] = np.where(a.obs.sample_type == "healthy_BM", "healthy", "AML_normal")
    mt = a.var_names.str.startswith("MT-")
    a.obs["pct_mt"] = np.asarray(a[:, mt].X.sum(1)).ravel() / a.obs.n_counts.values * 100
    q = p1["qc"]
    n0 = a.n_obs
    a = a[(a.obs.n_genes_detected >= q["min_genes"]) & (a.obs.pct_mt < q["max_pct_mt"])].copy()
    g0 = a.n_vars
    sc.pp.filter_genes(a, min_cells=q["min_cells_per_gene"])
    stats = {"n_cells_selected": int(n0), "n_cells_after_qc": a.n_obs, "n_genes_before": int(g0), "n_genes_after": a.n_vars,
             "umis_min": float(a.obs.n_counts.min()), "pct_mt_max": float(a.obs.pct_mt.max()), "n_mt_genes": int(mt.sum())}
    return a, stats


def preprocess(a: ad.AnnData, p1: dict, seed: int) -> tuple[PCA, np.ndarray, np.ndarray]:
    """Log-normalize, choose HVGs on counts, fit PCA, run Harmony. Returns (pca, hvg mean, hvg genes)."""
    # Samples with too few cells cannot be modelled as their own batch (loess / Harmony need data); pool them.
    n = a.obs[p1["harmony_batch"]].astype(str).map(a.obs[p1["harmony_batch"]].astype(str).value_counts())
    a.obs["batch"] = np.where(n >= p1["min_batch_cells"], a.obs[p1["harmony_batch"]].astype(str), "pooled_small_samples")
    a.layers["counts"] = a.X.copy()
    a.X = lognorm(a.X)
    # 'seurat' flavor (binned dispersion on log data) is used per batch; the loess fit of 'seurat_v3'
    # is numerically unstable within several of these small, sorted batches.
    sc.pp.highly_variable_genes(a, n_top_genes=p1["n_hvg"], flavor="seurat", batch_key="batch")
    genes = a.var_names[a.var.highly_variable].values
    Xh = a[:, genes].X.toarray()
    mean = Xh.mean(0)
    pca = PCA(n_components=p1["n_pcs"], random_state=seed, svd_solver="full").fit(Xh - mean)
    a.obsm["X_pca"] = pca.transform(Xh - mean)
    ho = harmonypy.run_harmony(a.obsm["X_pca"], a.obs, "batch", random_state=seed, verbose=False)
    a.obsm["X_harmony"] = np.asarray(ho.Z_corr).T
    return pca, mean, genes


def root_candidates(a: ad.AnnData, n: int) -> np.ndarray:
    """Healthy-donor HSCs ranked by the authors' independent HSC classifier score."""
    ok = np.where((a.obs.CellType == "HSC").values & (a.obs.source == "healthy").values)[0]
    return ok[np.argsort(-a.obs.Score_HSC.values[ok])][:n]


def run_dpt(a: ad.AnnData, rep: str, root: int, p1: dict, seed: int) -> np.ndarray:
    sc.pp.neighbors(a, n_neighbors=p1["n_neighbors"], use_rep=rep, random_state=seed, key_added="ref")
    sc.tl.diffmap(a, n_comps=15, neighbors_key="ref")
    a.uns["iroot"] = int(root)
    sc.tl.dpt(a, neighbors_key="ref")
    return a.obs["dpt_pseudotime"].values.astype(float)


def terminal_cells(a: ad.AnnData) -> pd.Series:
    """Terminal states from the authors' classifier scores (independent of our pseudotime)."""
    out = {}
    for ct, score in [("Mono", "Score_Mono"), ("cDC", "Score_cDC")]:
        m = (a.obs.CellType == ct).values
        out[a.obs_names[m][np.argmax(a.obs[score].values[m])]] = ct
    return pd.Series(out)


def run_palantir(a: ad.AnnData, rep: str, root: int, p1: dict, seed: int,
                 terminal: pd.Series | None = None) -> tuple[np.ndarray, pd.DataFrame]:
    import palantir
    b = ad.AnnData(obs=pd.DataFrame(index=a.obs_names))
    b.obsm["X_pca"] = a.obsm[rep]
    palantir.utils.run_diffusion_maps(b, n_components=10, knn=p1["n_neighbors"], seed=seed)
    palantir.utils.determine_multiscale_space(b)
    palantir.core.run_palantir(b, a.obs_names[root], terminal_states=terminal, knn=p1["n_neighbors"],
                               num_waypoints=p1["palantir_waypoints"], n_jobs=1, seed=seed)
    fate = b.obsm["palantir_fate_probabilities"]
    fate = fate if isinstance(fate, pd.DataFrame) else pd.DataFrame(fate, index=a.obs_names)
    return b.obs["palantir_pseudotime"].values.astype(float), fate


def order_agreement(pt: np.ndarray, ct: pd.Series, order: list[str]) -> float:
    rank = ct.map({c: i for i, c in enumerate(order)})
    ok = rank.notna().values
    return float(spearmanr(pt[ok], rank.values[ok])[0])


def to_score(pt: np.ndarray, q: tuple[float, float]) -> tuple[np.ndarray, tuple[float, float]]:
    lo, hi = np.quantile(pt, q)
    return np.clip((pt - lo) / (hi - lo), 0, 1), (float(lo), float(hi))


def embed(a: ad.AnnData, rep: str, seed: int) -> tuple[np.ndarray, str]:
    try:
        import phate
        Y = phate.PHATE(n_components=2, knn=15, random_state=seed, n_jobs=1, verbose=0).fit_transform(a.obsm[rep])
        return np.asarray(Y), "PHATE"
    except Exception as e:  # fall back to a force-directed layout of the same graph
        log.warning("PHATE failed (%s); using force-directed layout", e)
        sc.tl.draw_graph(a, layout="fr", neighbors_key="ref", random_state=seed)
        return a.obsm["X_draw_graph_fr"], "force-directed (FR)"


def potential_grid(xy: np.ndarray, score: np.ndarray, cfg: dict) -> dict:
    """Height surface on a regular grid.

    height = w * (1 - smoothed score) + (1 - w) * normalized(-log density), inside the data;
    it rises to a rim outside the data so the terrain is a closed basin.
    """
    n, m, bw, w = cfg["grid"], cfg["margin"], cfg["bandwidth"], cfg["weight_pseudotime"]
    lo, hi = xy.min(0), xy.max(0)
    span = hi - lo
    lo, hi = lo - m * span, hi + m * span
    gx, gy = np.linspace(lo[0], hi[0], n), np.linspace(lo[1], hi[1], n)
    # Bin cells, then Gaussian-smooth counts and score-weighted counts (a Nadaraya-Watson smoother).
    ix = np.clip(((xy[:, 0] - lo[0]) / (hi[0] - lo[0]) * (n - 1)).round().astype(int), 0, n - 1)
    iy = np.clip(((xy[:, 1] - lo[1]) / (hi[1] - lo[1]) * (n - 1)).round().astype(int), 0, n - 1)
    cnt = np.zeros((n, n)); s_sum = np.zeros((n, n))
    np.add.at(cnt, (iy, ix), 1); np.add.at(s_sum, (iy, ix), score)
    sigma = bw * n / (1 + 2 * m)
    dens = gaussian_filter(cnt, sigma); ssm = gaussian_filter(s_sum, sigma)
    score_g = np.where(dens > 1e-9, ssm / np.maximum(dens, 1e-12), np.nan)
    # Fill score far from data with the nearest smoothed value (wide kernel) so the surface is defined everywhere.
    wide = gaussian_filter(cnt, sigma * 4); wide_s = gaussian_filter(s_sum, sigma * 4)
    score_far = wide_s / np.maximum(wide, 1e-12)
    score_g = np.where(np.isnan(score_g) | (dens < 1e-6), score_far, score_g)
    density = dens / dens.sum()
    neglog = -np.log(density + density[density > 0].min())
    ref = neglog[dens > np.quantile(dens[dens > 0], 0.05)]
    neglog_n = np.clip((neglog - ref.min()) / (ref.max() - ref.min()), 0, 1)
    h_data = w * (1 - score_g) + (1 - w) * neglog_n
    support = dens / (dens + np.quantile(dens[dens > 0], 0.25))
    height = support * h_data + (1 - support) * 1.15
    height = gaussian_filter(height, sigma * 0.5)
    height = (height - height.min()) / (height.max() - height.min())
    return {"x": gx, "y": gy, "height": height, "score": np.clip(score_g, 0, 1), "density": density, "support": support}


def height_at(grid: dict, xy: np.ndarray) -> np.ndarray:
    from scipy.interpolate import RegularGridInterpolator
    f = RegularGridInterpolator((grid["y"], grid["x"]), grid["height"], bounds_error=False, fill_value=None)
    return f(np.c_[xy[:, 1], xy[:, 0]])
