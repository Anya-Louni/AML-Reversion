"""Phase 3: drug ranking by reversion score."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from pipeline.common import get_logger, path
from pipeline.phase3_drugs import scoring as S
from pipeline.phase3_drugs.report import write_report
from pipeline.phase3_drugs.signatures import consensus, fetch_wide

log = get_logger("phase3")


def run(cfg: dict) -> None:
    p3, seed = cfg["phase3"], cfg["seed"]
    proc = path(cfg, "processed", "x").parent
    raw = path(cfg, "raw", "lincs", "lincs2020", "x").parent
    sigs = fetch_wide(cfg, raw, proc)
    log.info("signatures: %d", sigs.n_obs)
    fs = sigs.var.feature_space
    space = {"landmark_best_inferred": fs.isin(["landmark", "best inferred"]), "landmark": fs.eq("landmark")}
    D_all, meta, per_cell = consensus(sigs, p3)
    cmap_to_drug = sigs.obs.groupby("cmap_name").drug.first().to_dict()
    st: dict = {"n_signatures": sigs.n_obs, "n_drugs": len(D_all), "tier_counts": meta.tier.value_counts().to_dict()}

    sig = pd.read_csv(proc / "maturation_signature.csv").set_index("gene")
    sig_full = pd.read_csv(proc / "maturation_signature_full.csv").set_index("gene")

    # ---------- generic maturation score of each drug (patient-independent)
    def generic(signature, genes_mask):
        genes = [g for g in signature.index if g in set(sigs.var_names[genes_mask.values])]
        w = signature.loc[genes, "spearman"].values
        cos = S.cosine(D_all[genes].values, w[None, :])[0]
        up = np.where(w > 0)[0]; dn = np.where(w < 0)[0]
        con = S.connectivity(D_all[genes].values, up, dn)
        return pd.DataFrame({"cosine": cos, "connectivity": con}, index=D_all.index), len(genes)

    G, n_g = generic(sig, space[p3["gene_space"]])
    G_lm, n_lm = generic(sig, space["landmark"])
    G_full, _ = generic(sig_full, space[p3["gene_space"]])
    G = G.join(meta)
    G["cosine_landmark"] = G_lm.cosine; G["cosine_full_signature"] = G_full.cosine
    G["rank_cosine"] = G.cosine.rank(ascending=False).astype(int)
    st["generic"] = {"n_genes": n_g, "n_genes_landmark": n_lm,
                     "cosine_vs_connectivity": float(spearmanr(G.cosine, G.connectivity)[0]),
                     "primary_vs_landmark": float(spearmanr(G.cosine, G.cosine_landmark)[0]),
                     "primary_vs_full_signature": float(spearmanr(G.cosine, G.cosine_full_signature)[0]),
                     "cosine_vs_tas": float(spearmanr(G.cosine, G.mean_tas, nan_policy="omit")[0])}

    # ---------- positive controls (evaluated on AML/monocytic-line consensus, all compounds as background)
    pc = {t: {cmap_to_drug.get(c, c) for c in p3["positive_controls"][t]} for t in ("tier1", "tier2")}
    bg = G[G.tier == "AML/monocytic"]
    ctrl_rows = []
    for t, names in pc.items():
        for d in sorted(names):
            if d in G.index:
                r = G.loc[d]
                in_bg = d in bg.index
                ctrl_rows.append({"tier": t, "drug": d, "line_tier": r.tier, "cells": r.cells, "cosine": r.cosine,
                                  "connectivity": r.connectivity,
                                  "percentile_in_myeloid_screen": float((bg.cosine < r.cosine).mean() * 100) if in_bg else np.nan,
                                  "rank_in_myeloid_screen": int((bg.cosine > r.cosine).sum() + 1) if in_bg else np.nan})
            else:
                ctrl_rows.append({"tier": t, "drug": d, "line_tier": "absent"})
    st["controls"] = pd.DataFrame(ctrl_rows)
    st["control_tests"] = {}
    for t in ("tier1", "tier2"):
        for metric in ("cosine", "connectivity"):
            pos = pc[t] & set(bg.index)
            st["control_tests"][f"{t}_{metric}"] = {
                "n_controls": len(pos), "n_background": len(bg),
                "auroc": S.auroc(bg[metric], pos),
                "perm_p": S.auroc_perm_p(bg[metric], pos, p3["n_permutations"], seed)}
    st["top_generic"] = G[G.tier == "AML/monocytic"].nlargest(25, "cosine")
    st["bottom_generic"] = G[G.tier == "AML/monocytic"].nsmallest(10, "cosine")

    # ---------- per-patient reversion scores
    genes_ok = set(sigs.var_names[space[p3["gene_space"]].values])
    match = pd.read_csv(proc / "drug_matching.csv").dropna(subset=["lincs_names"])
    beat_drugs = [d for d in match.beat_name if d in D_all.index]
    out = {}
    for cohort, file, drugs in [("vangalen", "stuck_signatures_vangalen.parquet", list(D_all.index)),
                                ("beataml", "stuck_signatures_beataml.parquet", beat_drugs)]:
        stuck = pd.read_parquet(proc / file)
        push = -stuck[[g for g in stuck.columns if g in genes_ok]]
        cos, con = S.score_patients(D_all.loc[drugs], push, p3["connectivity_set_size"])
        cos.to_parquet(proc / f"reversion_cosine_{cohort}.parquet")
        con.to_parquet(proc / f"reversion_connectivity_{cohort}.parquet")
        out[cohort] = (cos, con)
        st[f"{cohort}_n_genes"] = push.shape[1]
        st[f"{cohort}_method_agreement"] = float(np.median([spearmanr(cos.loc[p], con.loc[p])[0] for p in cos.index]))
        # how patient-specific are rankings? mean pairwise Spearman between patients' cosine rankings
        R = cos.T.rank()
        c = np.corrcoef(R.values.T)
        st[f"{cohort}_between_patient_rank_corr"] = float(c[np.triu_indices_from(c, 1)].mean())
    # ranked lists for van Galen patients (used by the web app)
    cos_vg, con_vg = out["vangalen"]
    ranked = []
    for p in cos_vg.index:
        df = pd.DataFrame({"cosine": cos_vg.loc[p], "connectivity": con_vg.loc[p]}).join(meta[["tier", "cells", "n_sigs"]])
        df["in_beataml_panel"] = df.index.isin(beat_drugs)
        df = df.sort_values("cosine", ascending=False)
        df["rank"] = np.arange(1, len(df) + 1)
        ranked.append(df.assign(patient=p).reset_index(names="drug"))
    ranked = pd.concat(ranked, ignore_index=True)
    ranked.to_parquet(proc / "drug_rankings_vangalen.parquet")
    G.to_parquet(proc / "drug_generic_scores.parquet")
    st["controls"].to_csv(proc / "phase3_positive_controls.csv", index=False)
    D_all.to_parquet(proc / "drug_consensus_signatures.parquet")
    meta.to_parquet(proc / "drug_signature_meta.parquet")
    st["ranked_vg"] = ranked
    ser = {k: v for k, v in st.items() if isinstance(v, (dict, int, float))}
    (proc / "phase3_stats.json").write_text(json.dumps(ser, indent=2, default=str))
    write_report(cfg, st, G, out)
    log.info("phase 3 done")
