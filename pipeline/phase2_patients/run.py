"""Phase 2: place patients on the landscape, deconvolve bulk samples, build stuck signatures."""
from __future__ import annotations

import json

import anndata as ad
import joblib
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from pipeline.common import get_logger, path
from pipeline.phase1_landscape.mapping import lognorm
from pipeline.phase2_patients import deconvolution as D
from pipeline.phase2_patients import placement as P
from pipeline.phase2_patients import stuck as S
from pipeline.phase2_patients.report import write_report

log = get_logger("phase2")


def run(cfg: dict) -> None:
    p2 = cfg["phase2"]
    proc = path(cfg, "processed", "x").parent
    st: dict = {}
    vg = ad.read_h5ad(proc / "vangalen_raw.h5ad")
    ref = ad.read_h5ad(proc / "landscape_reference.h5ad")
    mapper = joblib.load(proc / "landscape_mapper.joblib")
    sig = pd.read_csv(proc / "maturation_signature.csv")

    # ---------- maturity threshold, fixed on normal reference cells
    thr = P.maturity_threshold(ref.obs.score.values, ref.obs.CellType.astype(str), p2["immature_types"], p2["mature_types"])
    st["threshold"] = thr
    log.info("maturity threshold %.3f (balanced accuracy %.3f)", thr["threshold"], thr["balanced_accuracy"])

    # ---------- map all malignant AML cells (all days) + cell lines
    mal = vg.obs.PredictionRefined.eq("malignant").values & vg.obs.CellType.isin(p2["malignant_like_types"]).values
    m = vg[mal].copy()
    ref_Z = mapper.embed(mapper.align(ref.X, ref.var_names)[0])
    pred, info = P.map_cells(mapper, m.X, m.var_names, ref_Z)
    cut = float(np.quantile(info["ref_dist"], p2["offmanifold_quantile"]))
    for c in pred.columns:
        m.obs[c] = pred[c].values
    m.obs["off_manifold"] = m.obs.knn_dist > cut
    m.obs["stuck"] = m.obs.score < thr["threshold"]
    st["mapping"] = {"gene_fraction_present": info["gene_fraction_present"], "offmanifold_cut": cut,
                     "ref_dist_median": float(np.median(info["ref_dist"]))}

    aml = m.obs.sample_type.eq("AML").values
    d0 = aml & (m.obs.day == 0).values
    lo = p2["like_order"]
    st["like_order_rho_all"] = P.like_order_agreement(m.obs.score.values[aml], m.obs.CellType.astype(str)[aml], lo)
    st["like_order_rho_d0"] = P.like_order_agreement(m.obs.score.values[d0], m.obs.CellType.astype(str)[d0], lo)
    st["score_by_like_type"] = m.obs[aml].groupby("CellType", observed=True).score.describe()[["count", "25%", "50%", "75%"]] \
        .reindex(p2["malignant_like_types"])
    st["offmanifold_by_like_type"] = m.obs[aml].groupby("CellType", observed=True).off_manifold.mean().reindex(p2["malignant_like_types"])
    st["patients_d0"] = P.patient_summary(m.obs[d0], thr["threshold"], p2["min_cells_per_patient"])
    st["samples_all"] = m.obs[aml].groupby("sample", observed=True).agg(
        patient=("patient", "first"), day=("day", "first"), n=("score", "size"),
        median_score=("score", "median"), stuck_fraction=("stuck", "mean"))
    st["cell_lines"] = m.obs[~aml].groupby("sample", observed=True).agg(n=("score", "size"), median_score=("score", "median"),
                                                                      stuck_fraction=("stuck", "mean"), off_manifold=("off_manifold", "mean"))
    m.obs[["patient", "sample", "day", "sample_type", "CellType", "score", "x", "y", "knn_dist", "off_manifold", "stuck"]] \
        .to_parquet(proc / "malignant_cells_mapped.parquet")

    # ---------- single-cell stuck signatures (day-0 patients)
    sig_genes = [g for g in sig.gene if g in vg.var_names]
    mono = ref[ref.obs.CellType == "Mono"]
    mono_mean = np.asarray(mono[:, sig_genes].X.mean(0)).ravel()   # reference X is already log-normalized
    Xl = lognorm(m.X)
    gi = m.var_names.get_indexer(sig_genes)
    rows = {}
    for p, idx in m.obs[d0].groupby("patient", observed=True).indices.items():
        cells = np.where(d0)[0][idx]
        rows[p] = S.sc_stuck(Xl[cells][:, gi].toarray(), mono_mean, sig_genes)
    sc_stuck = pd.DataFrame(rows).T
    sc_stuck.to_parquet(proc / "stuck_signatures_vangalen.parquet")
    direction = sig.set_index("gene").loc[sig_genes].spearman.pipe(np.sign)
    st["sc_stuck_alignment"] = {p: float(spearmanr(-sc_stuck.loc[p], direction)[0]) for p in sc_stuck.index}

    # ---------- deconvolution reference, validated on pseudo-bulk with patient held out
    dc = p2["deconvolution"]
    labels = D.label_cells(vg.obs, p2["malignant_like_types"], dc["normal_groups"])
    X_cpm = D.cpm(vg.X)
    full_ref = D.build_reference(X_cpm, vg.var_names, labels, dc["n_markers_per_type"], dc["min_mean_cpm"])
    st["deconv_ref_shape"] = full_ref.shape
    val = []
    aml_samples = vg.obs[vg.obs.sample_type == "AML"].groupby("sample", observed=True).size()
    norm_to_like = {"HSC": "HSC-like", "Prog": "Prog-like", "GMP": "GMP-like", "ProMono": "ProMono-like", "Mono": "Mono-like", "cDC": "cDC-like"}
    for s in aml_samples[aml_samples >= 200].index:
        o = vg.obs[vg.obs["sample"] == s]
        pt = o.patient.iloc[0]
        keep = (vg.obs.patient != pt).values
        r = D.build_reference(X_cpm[keep], vg.var_names, labels[keep], dc["n_markers_per_type"], dc["min_mean_cpm"])
        pb = pd.DataFrame({s: np.asarray(X_cpm[(vg.obs["sample"] == s).values].sum(0)).ravel()}, index=vg.var_names)
        pb = pb / pb.sum() * 1e6
        est = D.malignant_shares(D.nnls_fractions(pb, r), p2["malignant_like_types"]).iloc[0]
        truth_lab = o.CellType.astype(str).where(o.PredictionRefined == "malignant", o.CellType.astype(str).map(norm_to_like))
        truth = truth_lab[truth_lab.isin(p2["malignant_like_types"])].value_counts(normalize=True).reindex(p2["malignant_like_types"]).fillna(0)
        for t in p2["malignant_like_types"]:
            val.append({"sample": s, "type": t, "true": truth[t], "est": est.get(t, np.nan)})
    val = pd.DataFrame(val)
    st["pseudobulk"] = val
    st["pseudobulk_overall_r"] = float(pearsonr(val.true, val.est)[0])
    coarse = val.assign(mature=val.type.isin(p2["mature_like_types"])).query("mature").groupby("sample")[["true", "est"]].sum()
    st["pseudobulk_mature_share"] = {"pearson": float(pearsonr(coarse.true, coarse.est)[0]),
                                     "mae": float(np.mean(np.abs(coarse.true - coarse.est))), "n_samples": len(coarse)}
    st["pseudobulk_by_type"] = val.groupby("type").apply(lambda g: pd.Series({
        "pearson": pearsonr(g.true, g.est)[0] if g.true.std() > 0 else np.nan,
        "mae": np.mean(np.abs(g.true - g.est))}), include_groups=False).reindex(p2["malignant_like_types"])

    # ---------- Beat AML deconvolution
    be = ad.read_h5ad(proc / "beataml_expr.h5ad")
    bulk = pd.DataFrame((2.0 ** be.X).T, index=be.var_names, columns=be.obs_names)
    frac = D.nnls_fractions(bulk, full_ref)
    shares = D.malignant_shares(frac, p2["malignant_like_types"])
    st["deconv_genes_used"] = int(len(full_ref.index.intersection(bulk.index)))
    out = frac.add_prefix("frac_").join(shares.add_prefix("share_"))
    out["mature_share"] = shares[[c for c in p2["mature_like_types"]]].sum(1)
    like_med = st["score_by_like_type"]["50%"]
    out["bulk_maturity_index"] = (shares * like_med.reindex(shares.columns).values).sum(1)
    out["rna_control"] = be.obs.rna_control.reindex(out.index).values
    out.to_parquet(proc / "beataml_deconvolution.parquet")

    vgs = pd.read_parquet(proc / "beataml_vg_celltype_scores.parquet").set_index("dbgap_rnaseq_sample")
    name = {"HSC-like": "HSC-like", "Prog-like": "Progenitor-like", "GMP-like": "GMP-like",
            "ProMono-like": "Promono-like", "Mono-like": "Monocyte-like", "cDC-like": "cDC-like"}
    common = out.index.intersection(vgs.index)
    st["vs_bottomly"] = pd.DataFrame([{"type": t, "n": len(common), "spearman": spearmanr(out.loc[common, f"frac_{t}"], vgs.loc[common, name[t]])[0]}
                                      for t in p2["malignant_like_types"]])

    # ---------- bulk stuck signatures
    aml_idx = out.index[(out.rna_control.isin(["No", "", "nan"])) | out.rna_control.isna()]
    E = pd.DataFrame(be[aml_idx].X, index=aml_idx, columns=be.var_names)
    bsig = [g for g in sig.gene if g in E.columns]
    fit = S.bulk_mature_expectation(E[bsig], out.loc[aml_idx, "mature_share"])
    bstuck = S.bulk_stuck(E[bsig], fit.expected_mature)
    bstuck.to_parquet(proc / "stuck_signatures_beataml.parquet")
    fit.to_csv(proc / "beataml_mature_expectation.csv")
    sgn = sig.set_index("gene").loc[bsig].spearman
    st["bulk_signature_transfer"] = {"n_genes": len(bsig),
                                     "sign_agreement": float(np.mean(np.sign(fit.slope) == np.sign(sgn))),
                                     "spearman_slope_vs_sc_rho": float(spearmanr(fit.slope, sgn)[0])}
    st["n_bulk_aml"] = len(aml_idx)
    (proc / "phase2_stats.json").write_text(json.dumps({k: v for k, v in st.items() if isinstance(v, (dict, float, int, tuple))},
                                                        indent=2, default=str))
    write_report(cfg, st, m, out, ref)
    log.info("phase 2 done")
