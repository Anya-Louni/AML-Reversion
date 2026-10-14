"""Phase 1: normal myeloid differentiation landscape."""
from __future__ import annotations

import json

import anndata as ad
import joblib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from pipeline.common import get_logger, path
from pipeline.phase1_landscape import landscape as L
from pipeline.phase1_landscape.mapping import LandscapeMapper, evaluate
from pipeline.phase1_landscape.report import write_report
from pipeline.phase1_landscape.signature import define_signature, spearman_with

log = get_logger("phase1")


def run(cfg: dict) -> None:
    p1, seed = cfg["phase1"], cfg["seed"]
    proc = path(cfg, "processed", "x").parent
    vg = ad.read_h5ad(proc / "vangalen_raw.h5ad")
    a, qc = L.select_reference(vg, p1)
    del vg
    log.info("reference: %d cells x %d genes", a.n_obs, a.n_vars)
    pca, mean, hvg = L.preprocess(a, p1, seed)
    ct = a.obs.CellType.astype(str)
    order = p1["lineage_order"]
    st: dict = {"qc": qc, "n_hvg": len(hvg), "pca_var_explained": float(pca.explained_variance_ratio_.sum())}

    # ---- batch structure before / after Harmony (kNN same-label fractions)
    from sklearn.neighbors import NearestNeighbors
    def mixing(rep):
        _, idx = NearestNeighbors(n_neighbors=31).fit(a.obsm[rep]).kneighbors(a.obsm[rep])
        out = {}
        for col in ["sample", "source"]:
            lab = np.asarray(a.obs[col].astype(str), dtype=object)
            out[col] = float((lab[idx[:, 1:]] == lab[:, None]).mean())
        lab = np.asarray(ct, dtype=object)
        out["CellType"] = float((lab[idx[:, 1:]] == lab[:, None]).mean())
        return out
    exp = {c: float((a.obs[c].astype(str).value_counts(normalize=True) ** 2).sum()) for c in ["sample", "source", "CellType"]}
    st["mixing"] = {"uncorrected_pca": mixing("X_pca"), "harmony": mixing("X_harmony"), "expected_if_random": exp}

    # ---- pseudotime: primary = DPT on Harmony PCs; second method = Palantir
    roots = L.root_candidates(a, p1["n_root_candidates"])
    dpt = L.run_dpt(a, "X_harmony", roots[0], p1, seed)
    if not np.isfinite(dpt).all():
        raise RuntimeError(f"{(~np.isfinite(dpt)).sum()} cells disconnected from the root in the kNN graph")
    root_rho = [float(spearmanr(dpt, L.run_dpt(a, "X_harmony", r, p1, seed))[0]) for r in roots[1:]]
    dpt_unc = L.run_dpt(a, "X_pca", roots[0], p1, seed)
    dpt = L.run_dpt(a, "X_harmony", roots[0], p1, seed)  # recompute primary last so the graph in `a` matches it
    log.info("running Palantir")
    pal_auto, fate_auto = L.run_palantir(a, "X_harmony", roots[0], p1, seed)
    term = L.terminal_cells(a)
    pal, fate = L.run_palantir(a, "X_harmony", roots[0], p1, seed, terminal=term)
    a.obs["dpt"], a.obs["dpt_uncorrected"], a.obs["palantir"] = dpt, dpt_unc, pal
    for c in fate.columns:
        a.obs[f"fate_{term.get(c, c)}"] = fate[c].values
    score, clip = L.to_score(dpt, tuple(p1["score_clip_quantiles"]))
    a.obs["score"] = score

    lin = ct.isin(order).values
    st["pseudotime"] = {
        "root_cell": a.obs_names[roots[0]], "root_sample": str(a.obs["sample"].iloc[roots[0]]),
        "order_rho": {"dpt": L.order_agreement(dpt, ct, order), "palantir": L.order_agreement(pal, ct, order),
                      "dpt_uncorrected": L.order_agreement(dpt_unc, ct, order)},
        "dpt_vs_palantir_spearman": float(spearmanr(dpt, pal)[0]),
        "dpt_vs_palantir_spearman_lineage_only": float(spearmanr(dpt[lin], pal[lin])[0]),
        "dpt_vs_uncorrected_spearman": float(spearmanr(dpt, dpt_unc)[0]),
        "root_sensitivity_min_rho": min(root_rho), "root_sensitivity_median_rho": float(np.median(root_rho)),
        "score_clip": clip,
        "palantir_auto_terminal_states": [str(ct.loc[c]) if c in ct.index else str(c) for c in fate_auto.columns],
        "palantir_auto_order_rho": L.order_agreement(pal_auto, ct, order),
        "palantir_auto_vs_fixed_spearman": float(spearmanr(pal_auto, pal)[0]),
        "palantir_terminal_states": {k: v for k, v in term.items()},
    }
    fate_cols = [c for c in a.obs.columns if c.startswith("fate_")]
    st["fate_by_type"] = a.obs.groupby("CellType", observed=True)[fate_cols].median().reindex(p1["cell_types"])
    med = a.obs.groupby(["CellType", "source"], observed=True)[["score", "palantir"]].median().unstack()
    st["median_by_type_source"] = med

    # ---- embedding and potential surface
    xy, method = L.embed(a, "X_harmony", seed)
    a.obsm["X_landscape"] = xy
    grid = L.potential_grid(xy, score, p1["landscape"])
    a.obs["height"] = L.height_at(grid, xy)
    st["embedding_method"] = method
    st["height_vs_score_spearman"] = float(spearmanr(a.obs.height, score)[0])
    centroids = pd.DataFrame(xy, index=a.obs_names).groupby(ct.values).median()
    st["centroids"] = centroids

    # ---- mapping expression -> (score, x, y)
    Xh = a[:, hvg].X.toarray()
    Y = np.c_[score, xy]
    ev = evaluate(Xh, Y, a.obs["sample"].astype(str).values, pca, mean, hvg, p1["mapping"]["k"],
                  p1["mapping"]["test_fraction"], seed, ct.values)
    mapper = LandscapeMapper.fit(Xh, hvg, pca, mean, Y, p1["mapping"]["k"])
    st["mapping"] = ev

    # ---- maturation signature (main lineage only, cDC excluded)
    corr = spearman_with(a.X[lin], score[lin], a.var_names)
    excl = corr.gene.str.match(p1["signature_exclude_regex"])
    corr["excluded_from_primary"] = excl.values
    sig = define_signature(corr[~excl], p1["signature_top_n"])
    sig_full = define_signature(corr, p1["signature_top_n"])
    st["signature"] = {"n_excluded_genes": int(excl.sum()),
                       "full_down_frac_excluded": float(sig_full[sig_full.direction.str.startswith("down")].gene.str.match(p1["signature_exclude_regex"]).mean()),
                       "full_up_frac_excluded": float(sig_full[sig_full.direction.str.startswith("up")].gene.str.match(p1["signature_exclude_regex"]).mean()),
                       "overlap_primary_full": int(len(set(sig.gene) & set(sig_full.gene)))}
    mk = p1["markers"]
    st["markers"] = {g: (float(corr.set_index("gene").spearman.get(g, np.nan)), "stem" if g in mk["stem"] else "mature")
                     for g in mk["stem"] + mk["mature"]}
    st["markers_pass"] = all((r < 0) if kind == "stem" else (r > 0) for r, kind in st["markers"].values())

    # ---- save
    a.write_h5ad(proc / "landscape_reference.h5ad", compression="gzip")
    np.savez_compressed(proc / "landscape_grid.npz", **grid)
    joblib.dump(mapper, proc / "landscape_mapper.joblib")
    corr.to_csv(proc / "maturation_correlations.csv", index=False)
    sig.to_csv(proc / "maturation_signature.csv", index=False)
    sig_full.to_csv(proc / "maturation_signature_full.csv", index=False)
    centroids.to_csv(proc / "landscape_centroids.csv")
    ser = {k: v for k, v in st.items() if k not in ("median_by_type_source", "centroids", "fate_by_type")}
    ser["mapping"] = {k: v for k, v in ev.items() if k != "loso"}
    (proc / "phase1_stats.json").write_text(json.dumps(ser, indent=2, default=str))
    write_report(cfg, a, grid, corr, sig, st, fate)
    log.info("phase 1 done")
