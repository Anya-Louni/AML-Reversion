"""Phase 5: simulate drug-driven shifts for van Galen patients and export compact data for the web app."""
from __future__ import annotations

import datetime as dt
import json
import re

import anndata as ad
import joblib
import numpy as np
import pandas as pd

from pipeline.common import ROOT, get_logger, path
from pipeline.phase1_landscape.mapping import lognorm
from pipeline.phase5_simulation_export import simulate as SIM

log = get_logger("phase5")
WEB = ROOT / "web" / "public" / "data"


def slug(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")


def q16(a: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return np.round((np.clip(a, lo, hi) - lo) / (hi - lo) * 65535).astype(np.uint16)


def r(x, n=4):
    return [round(float(v), n) for v in np.asarray(x).ravel()]


def clean_nan(o):
    """Replace NaN/inf with None recursively so the output is strict JSON."""
    if isinstance(o, dict):
        return {k: clean_nan(v) for k, v in o.items()}
    if isinstance(o, list):
        return [clean_nan(v) for v in o]
    if isinstance(o, (float, np.floating)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def validation_export(proc, seed: int) -> dict:
    s4 = json.loads((proc / "phase4_stats.json").read_text())
    s3 = json.loads((proc / "phase3_stats.json").read_text())
    models = pd.read_csv(proc / "phase4_model_table.csv")
    pdg = pd.read_csv(proc / "phase4_per_drug_m1.csv").dropna(subset=["spearman"])
    ctl = pd.read_csv(proc / "phase3_positive_controls.csv")
    idh = pd.read_csv(proc / "phase4_idh.csv") if (proc / "phase4_idh.csv").exists() else pd.DataFrame()
    R = pd.read_parquet(proc / "reversion_cosine_beataml.parquet")
    resp = pd.read_parquet(proc / "beataml_drug_response.parquet").dropna(subset=["dbgap_rnaseq_sample"])
    Y = resp.pivot_table(index="dbgap_rnaseq_sample", columns="inhibitor", values="auc", aggfunc="mean")
    pairs = R.stack().rename("s").to_frame().join((-Y).stack().rename("y"), how="inner").reset_index()
    pairs = pairs.dropna(subset=["s", "y"])
    pairs = pairs[pairs.iloc[:, 1].isin(pdg.drug)].sample(n=min(2500, len(pairs)), random_state=seed)
    G = pd.read_parquet(proc / "drug_generic_scores.parquet")
    bg = G[G.tier == "AML/monocytic"].cosine
    counts, edges = np.histogram(bg, bins=40)
    return {
        "summary": {k: s4[k] for k in ("n_samples", "n_patients", "n_drugs", "n_pairs", "P1", "P2", "claims")},
        "models": models.to_dict(orient="records"),
        "per_drug": pdg[["drug", "spearman", "q", "tier", "n"]].fillna({"tier": ""}).to_dict(orient="records"),
        "scatter": [{"s": round(a, 4), "y": round(b, 2), "d": d} for a, b, d in zip(pairs.s, pairs.y, pairs.iloc[:, 1])],
        "controls": [{"drug": c.drug, "tier": c.tier, "percentile": None if pd.isna(c.get("percentile_in_myeloid_screen")) else float(c.percentile_in_myeloid_screen),
                      "cosine": None if pd.isna(c.get("cosine")) else float(c.cosine)} for _, c in ctl.iterrows()],
        "control_tests": s3["control_tests"],
        "idh": idh.to_dict(orient="records"),
        "venetoclax": s4.get("venetoclax"),
        "generic_hist": [{"bin": float((edges[i] + edges[i + 1]) / 2), "count": int(c)} for i, c in enumerate(counts)],
    }


def display_height(grid: dict, inside: float = 0.35, falloff_cells: float = 14.0, depth: float = 0.55) -> np.ndarray:
    """Heightfield for rendering only: the landscape height where there are cells, extended outward from the
    data edge and falling off smoothly with distance, so the terrain reads as an island without cliffs.
    Cell positions and scores are unaffected (they use the same surface in the web app)."""
    from scipy.ndimage import distance_transform_edt, gaussian_filter
    h, sup = grid["height"], grid["support"]
    mask = sup >= inside
    dist, (iy, ix) = distance_transform_edt(~mask, return_indices=True)
    edge_h = h[iy, ix]                                   # nearest in-data height
    drop = depth * np.clip(dist / falloff_cells, 0, 1) ** 1.5
    out = np.where(mask, h, edge_h - drop)
    return gaussian_filter(out, 1.2)


def run(cfg: dict) -> None:
    p5, seed = cfg["phase5"], cfg["seed"]
    proc = path(cfg, "processed", "x").parent
    (WEB / "traj").mkdir(parents=True, exist_ok=True)
    (WEB / "patients").mkdir(parents=True, exist_ok=True)

    mapper = joblib.load(proc / "landscape_mapper.joblib")
    ref = ad.read_h5ad(proc / "landscape_reference.h5ad")
    grid = dict(np.load(proc / "landscape_grid.npz"))
    s2 = json.loads((proc / "phase2_stats.json").read_text())
    thr, cut = s2["threshold"]["threshold"], s2["mapping"]["offmanifold_cut"]
    D = pd.read_parquet(proc / "drug_consensus_signatures.parquet")
    meta = pd.read_parquet(proc / "drug_signature_meta.parquet")
    rank = pd.read_parquet(proc / "drug_rankings_vangalen.parquet")
    lincs_var = ad.read_h5ad(proc / "lincs_signatures_all.h5ad", backed="r").var
    reliable = set(lincs_var.index[lincs_var.feature_space.isin(["landmark", "best inferred"])])
    D = D[[g for g in D.columns if g in reliable]]
    vg = ad.read_h5ad(proc / "vangalen_raw.h5ad")
    tags = SIM.mutation_tags(vg.obs[vg.obs.sample_type == "AML"])

    # per-gene SD of reference log expression, in model gene order
    Xr, _ = mapper.align(ref.X, ref.var_names)
    gene_sd = Xr.std(0)

    # Beat AML context for badges
    resp = pd.read_parquet(proc / "beataml_drug_response.parquet").dropna(subset=["dbgap_rnaseq_sample"])
    beat_ctx = resp.groupby("inhibitor").agg(beat_median_auc=("auc", "median"), beat_n=("auc", "size"))
    perdrug = pd.read_csv(proc / "phase4_per_drug_m1.csv").set_index("drug") if (proc / "phase4_per_drug_m1.csv").exists() else pd.DataFrame()

    alphas = np.linspace(0, p5["alpha_max"], p5["n_frames"])
    xlo, xhi = float(grid["x"][0]), float(grid["x"][-1])
    ylo, yhi = float(grid["y"][0]), float(grid["y"][-1])
    rng = np.random.default_rng(seed)
    d0 = (vg.obs.sample_type == "AML") & (vg.obs.day == 0) & (vg.obs.PredictionRefined == "malignant") \
        & vg.obs.CellType.isin(cfg["phase2"]["malignant_like_types"])
    like_codes = {t: i for i, t in enumerate(cfg["phase2"]["malignant_like_types"])}
    patients, sim_rows = [], []
    groups = vg.obs[d0].groupby("patient", observed=True).indices
    # Common shift length (see config): the largest candidate norm that keeps, on average over all
    # patient x drug pairs, at most `max_offmanifold_end` of final-frame cells off the reference manifold.
    plan = {}
    for p, idx in groups.items():
        cells = np.where(d0.values)[0][idx]
        if len(cells) > p5["max_cells_per_patient"]:
            cells = np.sort(rng.choice(cells, p5["max_cells_per_patient"], replace=False))
        rk = rank[rank.patient == p].sort_values("rank")
        cand = rk[(rk.tier == "AML/monocytic") & ~rk.drug.str.startswith("BRD-")]
        drugs = list(dict.fromkeys(cand.drug.head(p5["top_n_drugs"]).tolist() + [d for d in p5["always_include"] if d in D.index]))
        plan[p] = (cells, drugs)
    scan = {}
    for norm in p5["shift_norm_candidates"]:
        fr = []
        for p, (cells, drugs) in plan.items():
            X0, _ = mapper.align(lognorm(vg.X[cells]), vg.var_names)
            for d in drugs:
                sh = SIM.drug_shift(mapper, D.loc[d], gene_sd)
                fr.append(SIM.offmanifold_fraction(mapper, np.clip(X0 + p5["alpha_max"] * sh / max(np.linalg.norm(sh), 1e-9) * norm, 0, None), cut))
        scan[norm] = float(np.mean(fr))
    ok = [n for n, v in scan.items() if v <= p5["max_offmanifold_end"]]
    target_norm = max(ok) if ok else min(scan)
    log.info("shift norm scan %s -> %s", scan, target_norm)
    for p, (cells, drugs) in plan.items():
        n_all = len(groups[p])
        rk = rank[rank.patient == p].sort_values("rank")
        X0, _ = mapper.align(lognorm(vg.X[cells]), vg.var_names)
        ct = vg.obs.CellType.iloc[cells].astype(str).map(like_codes).values.astype(int)
        pinfo = {"id": p, "n_malignant_d0": int(n_all), "n_shown": int(len(cells)), "mutations": tags.get(p, []),
                 "cell_types": ct.tolist(), "simulated_drugs": []}
        for d in drugs:
            raw_shift = SIM.drug_shift(mapper, D.loc[d], gene_sd)
            shift = raw_shift / max(np.linalg.norm(raw_shift), 1e-9) * target_norm
            traj = SIM.simulate(mapper, X0, shift, alphas)
            oc = SIM.outcome(traj, thr)
            oc["offmanifold_start"] = SIM.offmanifold_fraction(mapper, X0, cut)
            oc["offmanifold_end"] = SIM.offmanifold_fraction(mapper, np.clip(X0 + alphas[-1] * shift, 0, None), cut)
            oc["raw_shift_norm"] = float(np.linalg.norm(raw_shift))
            status = oc.pop("status")
            buf = np.stack([q16(traj[..., 1], xlo, xhi), q16(traj[..., 2], ylo, yhi), q16(traj[..., 0], 0, 1)], axis=-1)
            fn = f"{slug(p)}__{slug(d)}.bin"
            (WEB / "traj" / fn).write_bytes(buf.astype("<u2").tobytes())
            pinfo["simulated_drugs"].append({"drug": d, "file": f"traj/{fn}", "status": status.tolist(), **oc})
            sim_rows.append({"patient": p, "drug": d, **oc})
        s0 = traj[0, :, 0]
        pinfo["stuck_fraction"] = float(np.mean(s0 < thr))
        pinfo["median_score"] = float(np.median(s0))
        top = rk.head(40).copy()
        top = top.join(beat_ctx, on="drug").join(perdrug[["spearman"]].rename(columns={"spearman": "validation_rho"}), on="drug") \
            if len(perdrug) else top.join(beat_ctx, on="drug")
        pinfo["ranking"] = [{k: (None if (isinstance(v, float) and np.isnan(v)) else (round(v, 4) if isinstance(v, float) else v))
                             for k, v in row.items() if k != "patient"}
                            for row in top.to_dict(orient="records")]
        (WEB / "patients" / f"{slug(p)}.json").write_text(json.dumps(pinfo, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
        patients.append({k: pinfo[k] for k in ("id", "n_malignant_d0", "n_shown", "mutations", "stuck_fraction", "median_score")}
                        | {"file": f"patients/{slug(p)}.json", "top_drug": drugs[0]})
        log.info("%s: %d drugs simulated", p, len(drugs))

    # landscape + reference cells
    sub = rng.choice(ref.n_obs, size=min(3000, ref.n_obs), replace=False)
    rx = ref.obsm["X_landscape"][sub]
    cents = pd.read_csv(proc / "landscape_centroids.csv", index_col=0)
    landscape = {"n": int(grid["height"].shape[0]), "x": [xlo, xhi], "y": [ylo, yhi],
                 "height": r(grid["height"], 3), "support": r(grid["support"], 2), "height_vis": r(display_height(grid), 3), "threshold": thr,
                 "labels": [{"name": n, "x": float(v.iloc[0]), "y": float(v.iloc[1])} for n, v in cents.iterrows()],
                 "reference": {"x": r(rx[:, 0], 5), "y": r(rx[:, 1], 5), "score": r(ref.obs.score.values[sub], 3),
                               "type": ref.obs.CellType.astype(str).values[sub].tolist()}}
    (WEB / "landscape.json").write_text(json.dumps(landscape))
    manifest = {"generated": dt.date.today().isoformat(), "threshold": thr, "frames": p5["n_frames"], "alpha_max": p5["alpha_max"],
                "like_types": cfg["phase2"]["malignant_like_types"], "quantization": {"x": [xlo, xhi], "y": [ylo, yhi], "score": [0, 1]},
                "traj_layout": "uint16 little-endian, [frame][cell][x, y, score]",
                "status_codes": {"0": "already mature", "1": "matured", "2": "resistant", "3": "moved but not past threshold at the end"},
                "patients": sorted(patients, key=lambda d: -d["stuck_fraction"]),
                "disclaimer": "Research simulation. In-silico hypotheses from a linear expression-shift model. Not for clinical use."}
    (WEB / "manifest.json").write_text(json.dumps(manifest, indent=1))
    (WEB / "validation.json").write_text(json.dumps(clean_nan(validation_export(proc, seed)), allow_nan=False))
    sims = pd.DataFrame(sim_rows)
    sims.to_csv(proc / "phase5_simulation_summary.csv", index=False)
    size = sum(f.stat().st_size for f in WEB.rglob("*") if f.is_file())
    log.info("exported %.1f MB to %s", size / 1e6, WEB)
    from pipeline.phase5_simulation_export.report import write_report
    write_report(cfg, sims, size, thr, {"scan": scan, "target_norm": target_norm})
