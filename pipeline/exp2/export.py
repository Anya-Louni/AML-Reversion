"""Experiment 2 outputs for the web app: hierarchy-based drug ranking for the van Galen patients and a
validation summary. Uses only the development-fitted primary model (cell-type sensitivities) and the saved
single test evaluation; nothing here refits on or re-scores the test set.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from pipeline.common import path
from pipeline.phase2_patients.deconvolution import label_cells

NORMAL_TO_LIKE = {"HSC": "HSC-like", "Prog": "Prog-like", "GMP": "GMP-like", "ProMono": "ProMono-like", "Mono": "Mono-like", "cDC": "cDC-like"}


def sample_fractions(obs: pd.DataFrame, like_types, normal_groups: dict) -> pd.Series:
    """Cell-type fractions of one single-cell sample in the same 10 categories as the bulk deconvolution:
    malignant like-types plus normal myeloid cells counted toward the matching like-type (as in bulk), and
    the non-myeloid normal groups."""
    lab = label_cells(obs, like_types, normal_groups)
    norm_my = (obs.PredictionRefined == "normal") & obs.CellType.isin(list(NORMAL_TO_LIKE))
    lab[norm_my.values] = obs.CellType[norm_my].astype(str).map(NORMAL_TO_LIKE).values
    cats = list(like_types) + list(normal_groups)
    return lab.dropna().value_counts(normalize=True).reindex(cats).fillna(0.0)


def hierarchy_predictions(fractions: pd.Series, S: pd.DataFrame, means: pd.Series) -> pd.DataFrame:
    """Predicted response = sum_k fraction(k) * sensitivity(k, d); residual = response - development drug mean."""
    f = fractions.reindex(S.columns).fillna(0.0)
    resp = S.values @ f.values
    return pd.DataFrame({"pred_response": resp, "pred_residual": resp - means.reindex(S.index).values}, index=S.index)


def web_payload(cfg, vg_obs: pd.DataFrame, patients: list[str]) -> dict:
    proc = path(cfg, "processed", "x").parent
    S = pd.read_csv(proc / "exp2_celltype_sensitivity.csv", index_col=0)
    order = list(cfg["phase2"]["malignant_like_types"]) + list(cfg["phase2"]["deconvolution"]["normal_groups"])
    S = S[[c for c in order if c in S.columns]]   # maturation order, then non-myeloid groups
    from pipeline.exp2.step3_run import assemble
    means = assemble(cfg, proc)["means"]
    p2 = cfg["phase2"]
    out = {"cell_types": list(S.columns), "sensitivity": {d: [round(float(v), 1) for v in S.loc[d]] for d in S.index},
           "drug_means": {d: round(float(means[d]), 1) for d in S.index}, "patients": {}}
    for p in patients:
        obs = vg_obs[(vg_obs.patient == p) & (vg_obs.day == 0) & (vg_obs.sample_type == "AML")]
        f = sample_fractions(obs, p2["malignant_like_types"], p2["deconvolution"]["normal_groups"])
        pr = hierarchy_predictions(f.rename(lambda c: c), S, means).sort_values("pred_residual", ascending=False)
        out["patients"][p] = {"fractions": [round(float(f.get(c, 0.0)), 4) for c in S.columns],
                              "ranking": [{"drug": d, "pred_residual": round(float(r.pred_residual), 2)} for d, r in pr.iterrows()]}
    return out


def validation_payload(cfg) -> dict:
    proc = path(cfg, "processed", "x").parent
    s1 = json.loads((proc / "exp2_step1_stats.json").read_text())
    r = json.loads((proc / "exp2_test_results.json").read_text())
    return {"ceiling": {"central": s1["central_ceiling"], "bias_corrected": s1["split_replicate_debiased"]["ceiling"],
                        "test_retest": s1["ceiling_test_retest"], "variance_model": s1["reliability_replicates"]["ceiling_corr"],
                        "shares": s1["shares"]},
            "primary": r["primary"], "models": r["table"], "secondary": r["secondary"],
            "n_test_patients": r["n_test_patients"], "n_test_samples": r["n_test_samples"],
            "phase4_descriptive": r["phase4_descriptive"], "venetoclax_sensitivity": r["venetoclax_sensitivity"]}
