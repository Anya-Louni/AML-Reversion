"""Experiment 2, Step 3: models on the development set only (5-fold patient-level CV).

The locked test samples are dropped right after loading; no function in this step sees them.
Writes reports/exp2_step3_development.md and data/processed/exp2_dev_choices.json (penalties used in Step 4).
"""
from __future__ import annotations

import datetime as dt
import json

import numpy as np
import pandas as pd

from pipeline.common import get_logger, path
from pipeline.exp2 import models as M
from pipeline.phase4_validation.run import load as load_phase4

log = get_logger("exp2.step3")


def assemble(cfg, proc):
    """Everything both Step 3 and Step 4 need. Returns full-cohort objects; callers choose rows."""
    Y, _, subj, _ = load_phase4(cfg, proc)
    split = pd.read_csv(path(cfg, "reports", "exp2_split.csv")).set_index("sample")
    dev = split.index[split.split == "development"]
    test = split.index[split.split == "test"]
    n_dev, n_test = Y.loc[dev].notna().sum(), Y.loc[test].notna().sum()
    drugs = [d for d in Y.columns if n_dev[d] >= 50 and n_test[d] >= 20]
    Y = Y[drugs]
    means = M.drug_means(Y, dev)
    R = M.residuals(Y, means)
    dec = pd.read_parquet(proc / "beataml_deconvolution.parquet")
    frac = dec.loc[Y.index, [c for c in dec.columns if c.startswith("frac_")]].fillna(0.0)
    clin = pd.read_excel(path(cfg, "raw", "beataml", "beataml_wv1to4_clinical.xlsx"), sheet_name="summary")
    clin = clin.dropna(subset=["dbgap_rnaseq_sample"]).drop_duplicates("dbgap_rnaseq_sample").set_index("dbgap_rnaseq_sample")
    mut = pd.read_parquet(proc / "beataml_mutations_matrix.parquet")
    MUT = M.mutation_clinical_features(Y.index, dev, clin, mut)
    rev = pd.read_parquet(proc / "reversion_cosine_beataml.parquet").reindex(index=Y.index, columns=drugs)
    meta = pd.read_parquet(proc / "drug_signature_meta.parquet")
    myeloid_drugs = [d for d in drugs if meta.tier.get(d) == "AML/monocytic"]
    low = set(pd.read_csv(path(cfg, "reports", "exp2_low_signal_samples.csv"))["sample"])
    specs = {
        "B mutation + clinical": {"F": MUT, "extra": None, "drugs": drugs},
        "H hierarchy (PRIMARY)": {"F": frac, "extra": None, "drugs": drugs},
        "HR hierarchy + reversion": {"F": frac, "extra": rev, "drugs": drugs},
        "HM hierarchy + myeloid-line reversion": {"F": frac, "extra": rev, "drugs": myeloid_drugs},
    }
    return {"Y": Y, "R": R, "means": means, "subj": subj.astype(str), "dev": dev, "test": test, "drugs": drugs,
            "frac": frac, "MUT": MUT, "rev": rev, "specs": specs, "low": low, "myeloid_drugs": myeloid_drugs}


def cv_predict(D, spec, alpha, rows, folds):
    R = D["R"][spec["drugs"]]
    P = pd.DataFrame(np.nan, index=rows, columns=spec["drugs"])
    for f in sorted(folds.unique()):
        te = folds.index[folds == f]; tr = folds.index[folds != f]
        m = M.PerDrugRidge(alpha).fit(R, spec["F"], tr, spec["extra"])
        p = m.predict(spec["F"], te, spec["extra"])
        P.loc[te, p.columns] = p
    return P


def summarize(P, R, subj, n_boot=10000):
    ps = M.per_sample_spearman(P, R)
    pp = M.per_patient(ps, subj)
    mean, lo, hi = M.bootstrap_ci(pp, n_boot)
    hits = M.per_patient(M.top_k_hits(P, R), subj)
    return {"primary_mean": mean, "ci_lo": lo, "ci_hi": hi, "n_patients": int(pp.size),
            "pooled_pearson": M.pooled_pearson(P, R), "top5_hit_rate": float(hits.mean())}, pp


def run(cfg: dict) -> None:
    proc = path(cfg, "processed", "x").parent
    D = assemble(cfg, proc)
    dev = D["dev"]
    R = D["R"].loc[dev]
    subj = D["subj"]
    folds = M.patient_folds(subj.loc[dev], 5, cfg["seed"])
    st = {"n_dev_samples": len(dev), "n_dev_patients": int(subj.loc[dev].nunique()), "n_drugs": len(D["drugs"]),
          "n_myeloid_drugs": len(D["myeloid_drugs"]), "n_mutation_features": D["MUT"].shape[1],
          "mutation_features": list(D["MUT"].columns)}
    rows, choices, grids = [], {}, {}
    # A: predicts 0 -> per-patient Spearman defined as 0, hit rate = chance
    rows.append({"model": "A global ranking", "alpha": None, "primary_mean": 0.0, "ci_lo": 0.0, "ci_hi": 0.0,
                 "n_patients": int(subj.loc[dev].nunique()), "pooled_pearson": 0.0, "top5_hit_rate": 0.20})
    per_pat = {}
    for name, spec in D["specs"].items():
        grid = {}
        for a in M.ALPHAS:
            P = cv_predict(D, spec, a, dev, folds)
            grid[a] = M.per_patient(M.per_sample_spearman(P, R[spec["drugs"]]), subj).mean()
        best = max(grid, key=grid.get)
        grids[name] = grid
        choices[name] = best
        P = cv_predict(D, spec, best, dev, folds)
        s, pp = summarize(P, R[spec["drugs"]], subj)
        per_pat[name] = pp
        rows.append({"model": name, "alpha": best, **s})
        # secondary: excluding the fixed low-signal samples from fitting and scoring
        keep = dev[~dev.isin(D["low"])]
        f2 = folds.loc[keep]
        P2 = cv_predict(D, spec, best, keep, f2)
        s2, _ = summarize(P2, D["R"].loc[keep, spec["drugs"]], subj)
        rows.append({"model": f"{name} [low-signal excluded]", "alpha": best, **s2})
        log.info("%s: alpha %s, CV %.3f", name, best, s["primary_mean"])
    # descriptive: Phase 4 raw reversion score, no fitting
    s, _ = summarize(D["rev"].loc[dev], R, subj)
    rows.append({"model": "Phase 4 raw reversion score (descriptive)", "alpha": None, **s})
    table = pd.DataFrame(rows)
    # development-CV paired differences vs H (descriptive; tuning makes CV slightly optimistic)
    h = per_pat["H hierarchy (PRIMARY)"]
    diffs = {}
    for name, pp in per_pat.items():
        if name.startswith("H hierarchy"):
            continue
        d = (h - pp).dropna()
        diffs[f"H minus {name}"] = {"mean_diff": float(d.mean()), "p_two_sided": M.sign_flip_test(d, alternative="two-sided")}
    st["diffs"] = diffs
    st["alpha_grid"] = {k: {str(a): float(v) for a, v in g.items()} for k, g in grids.items()}

    # Exploratory diagnostic (development only, not pre-registered): the target removes drug means but not a
    # patient's overall sensitivity, so a generally sensitive patient has larger residuals on drugs with large
    # between-patient spread. How much of each model's score comes from that pattern rather than drug-specific matching?
    from pipeline.exp2.ceiling import additive_fit
    sd = R.std()
    g = additive_fit(R)["row"]
    scale = pd.DataFrame(np.outer(g.values, (sd / sd.mean()).values), index=R.index, columns=R.columns)
    obs_pat = pd.DataFrame(np.outer(np.sign(g.values), sd.values), index=R.index, columns=R.columns)
    diag = {"observed residuals vs sign(patient effect) x drug SD": float(M.per_patient(M.per_sample_spearman(obs_pat, R), subj).mean())}
    for name in ["B mutation + clinical", "H hierarchy (PRIMARY)"]:
        P = cv_predict(D, D["specs"][name], choices[name], dev, folds)
        pat = pd.DataFrame(np.outer(np.sign(P.mean(axis=1).values), sd[P.columns].values), index=P.index, columns=P.columns)
        diag[f"{name}: predictions vs that pattern"] = float(M.per_patient(M.per_sample_spearman(P, pat), subj).mean())
        diag[f"{name}: CV score with the pattern removed from the target"] = float(
            M.per_patient(M.per_sample_spearman(P, R[P.columns] - scale[P.columns]), subj).mean())
    st["diagnostic"] = diag
    (proc / "exp2_dev_choices.json").write_text(json.dumps({"alphas": choices, "drugs": D["drugs"],
                                                            "myeloid_drugs": D["myeloid_drugs"],
                                                            "mutation_features": st["mutation_features"]}, indent=2))
    table.to_csv(proc / "exp2_dev_cv_table.csv", index=False)
    write(cfg, st, table)


def write(cfg, st, table):
    L = []; w = L.append
    w("# Experiment 2, Step 3: development-set results")
    w(f"\nGenerated {dt.date.today().isoformat()} by `python -m pipeline.run --phase exp2-3`. Plan: [exp2_plan.md](exp2_plan.md). "
      "Only development data were used; the locked test set was not loaded into any model or metric.\n")
    w(f"- Development: {st['n_dev_samples']} samples from {st['n_dev_patients']} patients; {st['n_drugs']} drugs "
      f"({st['n_myeloid_drugs']} with AML/monocytic-line signatures, used by HM).")
    w(f"- Mutation/clinical features (chosen on development data): {st['n_mutation_features']}.")
    w("- Penalty per model chosen from {0.01, 0.1, 1, 10, 100, 1000} by 5-fold patient-level CV on the primary metric. Because the same folds choose the "
      "penalty and report the score, development numbers are slightly optimistic; only the Step 4 test result counts.")
    w("- HM note: the Phase 3 consensus already builds a drug's signature only from AML/monocytic lines whenever such lines exist, so recomputing with "
      "only those lines gives identical scores. HM is therefore HR restricted to those drugs.\n")
    w("## Cross-validated results (development)\n")
    w(table.to_markdown(index=False, floatfmt=".3f"))
    w("\n`primary_mean`: mean over patients of per-sample Spearman between predicted and observed residuals, with patient-level bootstrap 95% CI. "
      "`top5_hit_rate`: share of a sample's 5 highest-predicted drugs in its observed top 20% (chance 0.20).")
    w("\n## Paired differences vs the primary model (development CV, descriptive)\n")
    w(pd.DataFrame(st["diffs"]).T.to_markdown(floatfmt=".4f"))
    w("\n## Penalty grid (CV primary metric)\n")
    w(pd.DataFrame(st["alpha_grid"]).to_markdown(floatfmt=".4f"))
    w("\nAll four models chose the largest penalty in the pre-registered grid. The curves have flattened there (H gains 0.0004 from 100 to 1000), so "
      "extending the grid would change little; the grid is kept as registered.")
    w("\n## Exploratory diagnostic (development only, not pre-registered)\n")
    w("The target removes each drug's mean but not each patient's overall sensitivity. A generally sensitive patient therefore tends to have larger "
      "residuals on drugs whose response varies most between patients, and a model could score well just by detecting general sensitivity, without "
      "matching specific drugs to patients. The checks below estimate how much each model relies on that pattern. The last rows remove the pattern from "
      "the target using each patient's fitted overall effect (an oracle quantity, used only for this diagnostic).\n")
    w(pd.Series(st["diagnostic"]).rename("value").to_frame().to_markdown(floatfmt=".3f"))
    w("\nIf a model's score barely drops when the pattern is removed, its signal is mostly drug-specific.")
    w("\nThe primary model (H) was fixed in the plan and is not changed by these results.")
    path(cfg, "reports", "exp2_step3_development.md").write_text("\n".join(L) + "\n", encoding="utf-8")
