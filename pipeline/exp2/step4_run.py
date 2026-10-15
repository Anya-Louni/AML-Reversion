"""Experiment 2, Step 4: the single evaluation on the locked test set (plan: reports/exp2_plan.md).

Runs only with EXP2_CONFIRM_TEST=yes and only if the test set has not been evaluated before.
Penalties and drug lists are read from the committed reports/exp2_dev_choices.json.
"""
from __future__ import annotations

import datetime as dt
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pipeline.common import get_logger, path
from pipeline.exp2 import models as M
from pipeline.exp2.split import assert_test_not_evaluated, mark_test_evaluated
from pipeline.exp2.step3_run import assemble

log = get_logger("exp2.step4")
N_BOOT = N_PERM = 10000
N_LABEL_PERM = 2000
CEILINGS = {"pre-specified": 0.67, "bias-corrected": 0.77}
PRIMARY = "H hierarchy (PRIMARY)"


def evaluate(P, R, subj):
    ps = M.per_sample_spearman(P, R)
    pp = M.per_patient(ps, subj)
    mean, lo, hi = M.bootstrap_ci(pp, N_BOOT)
    hits = M.per_patient(M.top_k_hits(P, R), subj)
    return {"primary_mean": mean, "ci_lo": lo, "ci_hi": hi, "n_patients": int(pp.size), "n_samples": int(ps.size),
            "pooled_pearson": M.pooled_pearson(P, R), "top5_hit_rate": float(hits.mean())}, pp


def run(cfg: dict) -> None:
    if os.environ.get("EXP2_CONFIRM_TEST") != "yes":
        raise RuntimeError("set EXP2_CONFIRM_TEST=yes to run the single test evaluation")
    proc = path(cfg, "processed", "x").parent
    assert_test_not_evaluated(proc)
    choices = json.loads(path(cfg, "reports", "exp2_dev_choices.json").read_text())
    D = assemble(cfg, proc)
    if D["drugs"] != choices["drugs"] or D["myeloid_drugs"] != choices["myeloid_drugs"] or list(D["MUT"].columns) != choices["mutation_features"]:
        raise RuntimeError("data assembly differs from the locked development choices")
    dev, test, subj = D["dev"], D["test"], D["subj"]
    R = D["R"]
    res, pp, fitted = {}, {}, {}

    # ---- all models: fit on the whole development set, predict the test set once
    for name, spec in D["specs"].items():
        m = M.PerDrugRidge(choices["alphas"][name]).fit(R[spec["drugs"]], spec["F"], dev, spec["extra"])
        fitted[name] = m
        P = m.predict(spec["F"], test, spec["extra"])
        res[name], pp[name] = evaluate(P, R.loc[test, spec["drugs"]], subj)
        res[name]["alpha"] = choices["alphas"][name]
        if name == PRIMARY:
            P_primary = P
    test_patients = subj.loc[test].unique()
    res["A global ranking"] = {"primary_mean": 0.0, "ci_lo": 0.0, "ci_hi": 0.0, "n_patients": res[PRIMARY]["n_patients"],
                               "n_samples": res[PRIMARY]["n_samples"], "pooled_pearson": 0.0, "top5_hit_rate": 0.20, "alpha": None}

    # ---- primary test: H vs A (A = 0 for every patient), paired sign-flip, one-sided
    h = pp[PRIMARY]
    p_primary = M.sign_flip_test(h - 0.0, N_PERM, cfg["seed"], "greater")
    supported = bool(p_primary < 0.05 and res[PRIMARY]["ci_lo"] > 0)

    # ---- secondary comparisons (two-sided)
    sec = {}
    sec["H vs B"] = {"mean_diff": float((h - pp["B mutation + clinical"]).mean()),
                     "p": M.sign_flip_test((h - pp["B mutation + clinical"]).dropna(), N_PERM, cfg["seed"], "two-sided")}
    d = (pp["HR hierarchy + reversion"] - h).dropna()
    sec["HR vs H"] = {"mean_diff": float(d.mean()), "p": M.sign_flip_test(d, N_PERM, cfg["seed"], "two-sided")}
    # HM uses a drug subset; compare with H restricted to the same drugs
    my = D["myeloid_drugs"]
    _, pp_h_my = evaluate(P_primary[my], R.loc[test, my], subj)
    d = (pp["HM hierarchy + myeloid-line reversion"] - pp_h_my).dropna()
    sec["HM vs H (same drug subset)"] = {"mean_diff": float(d.mean()), "p": M.sign_flip_test(d, N_PERM, cfg["seed"], "two-sided"),
                                         "H_on_subset": float(pp_h_my.mean())}
    p_label = M.label_permutation_test(P_primary, R.loc[test, P_primary.columns], subj, N_LABEL_PERM, cfg["seed"])

    # ---- amendment 1a: low-signal samples excluded from development fitting and test scoring
    low = D["low"]
    dev2, test2 = dev[~dev.isin(low)], test[~test.isin(low)]
    low_res = {}
    for name, spec in D["specs"].items():
        m = M.PerDrugRidge(choices["alphas"][name]).fit(R[spec["drugs"]], spec["F"], dev2, spec["extra"])
        s_, pp_ = evaluate(m.predict(spec["F"], test2, spec["extra"]), R.loc[test2, spec["drugs"]], subj)
        if name == PRIMARY:
            s_["p_vs_A"] = M.sign_flip_test(pp_, N_PERM, cfg["seed"], "greater")
        low_res[name] = s_

    # ---- amendment 1b: Phase 4 raw reversion score, descriptive
    p4, pp_p4 = evaluate(D["rev"].loc[test], R.loc[test], subj)

    # ---- sensitivity heatmap from the development-fitted primary model
    S = fitted[PRIMARY].pure_type_sensitivity(D["frac"], D["means"])
    S.columns = [c.replace("frac_", "") for c in S.columns]
    ven = S.loc["Venetoclax"] if "Venetoclax" in S.index else None

    # ---- lock, then save and report
    mark_test_evaluated(proc)
    table = pd.DataFrame(res).T.loc[["A global ranking", "B mutation + clinical", PRIMARY, "HR hierarchy + reversion",
                                     "HM hierarchy + myeloid-line reversion"]]
    out = {"primary": {"mean": res[PRIMARY]["primary_mean"], "ci": [res[PRIMARY]["ci_lo"], res[PRIMARY]["ci_hi"]], "p_vs_A": p_primary,
                       "supported": supported, "p_label_permutation": p_label,
                       "fraction_of_ceiling": {k: res[PRIMARY]["primary_mean"] / v for k, v in CEILINGS.items()}},
           "table": table.reset_index(names="model").to_dict(orient="records"), "secondary": sec,
           "low_signal_excluded": low_res, "phase4_descriptive": p4,
           "n_test_samples": len(test), "n_test_patients": int(len(test_patients)),
           "n_test_low_signal_removed": int(len(test) - len(test2)),
           "venetoclax_sensitivity": ven.round(2).to_dict() if ven is not None else None}
    (proc / "exp2_test_results.json").write_text(json.dumps(out, indent=2, default=float))
    S.to_csv(proc / "exp2_celltype_sensitivity.csv")
    F = figures(cfg, P_primary, R.loc[test], pp, h, S)
    write(cfg, out, table, F, S)
    log.info("test evaluated once; primary %.3f [%.3f, %.3f], p = %.4f, supported = %s",
             res[PRIMARY]["primary_mean"], res[PRIMARY]["ci_lo"], res[PRIMARY]["ci_hi"], p_primary, supported)


def figures(cfg, P, R, pp, h, S):
    F = {}
    cols = P.columns
    a, b = P[cols].values.ravel(), R[cols].values.ravel()
    ok = ~np.isnan(a) & ~np.isnan(b)
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    ax[0].hexbin(a[ok], b[ok], gridsize=55, cmap="Blues", mincnt=1)
    ax[0].axhline(0, color="k", lw=0.5); ax[0].axvline(0, color="k", lw=0.5)
    ax[0].set_xlabel("predicted residual (H)"); ax[0].set_ylabel("observed residual (sensitivity - dev drug mean)")
    ax[0].set_title(f"Locked test set: {ok.sum():,} patient-drug pairs")
    data = [pp["B mutation + clinical"].values, h.values]
    ax[1].violinplot(data, positions=[1, 2], showmedians=True)
    ax[1].axhline(0, color="#C44E52", ls="--", lw=1, label="A global ranking (0)")
    ax[1].set_xticks([1, 2], ["B mutation + clinical", "H hierarchy (primary)"])
    ax[1].set_ylabel("per-patient Spearman"); ax[1].set_title("Per-patient correlation on the test set"); ax[1].legend(frameon=False, fontsize=8)
    fig.tight_layout(); F["test"] = path(cfg, "figures", "exp2_test_results.png"); fig.savefig(F["test"], dpi=150); plt.close(fig)

    order = ["HSC-like", "Prog-like", "GMP-like", "ProMono-like", "Mono-like", "cDC-like", "T_NK", "B", "Ery", "pDC"]
    Z = S[[c for c in order if c in S.columns]]
    Zc = Z.sub(Z.mean(axis=1), axis=0)
    Zc = Zc.loc[Zc.var(axis=1).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(8, 18))
    v = np.nanpercentile(np.abs(Zc.values), 98)
    im = ax.imshow(Zc.values, cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
    ax.set_xticks(range(Zc.shape[1]), Zc.columns, rotation=60, fontsize=8)
    ax.set_yticks(range(len(Zc)), Zc.index, fontsize=5)
    plt.colorbar(im, ax=ax, fraction=0.04, label="sensitivity relative to drug's row mean (-AUC units)")
    ax.set_title("Learned cell-type sensitivity (primary model, development fit)\nred = more sensitive", fontsize=10)
    fig.tight_layout(); F["heat"] = path(cfg, "figures", "exp2_celltype_sensitivity.png"); fig.savefig(F["heat"], dpi=150); plt.close(fig)
    return F


def write(cfg, out, table, F, S):
    pr, sec = out["primary"], out["secondary"]
    L = []; w = L.append
    w("# Experiment 2, Step 4: locked test set results")
    w(f"\nGenerated {dt.date.today().isoformat()} by the single test evaluation (`EXP2_CONFIRM_TEST=yes python -m pipeline.run --phase exp2-4`). "
      "Plan: [exp2_plan.md](exp2_plan.md), including Amendment 1. Development results: [exp2_step3_development.md](exp2_step3_development.md). "
      "Experiment 2 is a follow-up designed after the first pre-registered validation (Phase 4) found no patient-specific value.\n")
    w("## Primary result\n")
    w(f"- Test set: {out['n_test_samples']} samples from {out['n_test_patients']} patients (Beat AML waves 3+4 only), {len(S)} drugs scored by the primary model.")
    w(f"- **Primary model H (hierarchy-weighted): mean per-patient Spearman {pr['mean']:.3f}, 95% CI {pr['ci'][0]:.3f} to {pr['ci'][1]:.3f}.**")
    w(f"- **H vs A (global ranking), paired sign-flip permutation, one-sided: p = {pr['p_vs_A']:.4f}.**")
    w(f"- **Decision (pre-registered rule: p < 0.05 and CI excluding 0): patient-specific value {'SUPPORTED' if pr['supported'] else 'NOT supported'}.**")
    w(f"- Robustness, within-patient drug-label permutation ({N_LABEL_PERM:,} permutations): p = {pr['p_label_permutation']:.4f}.")
    w(f"- Relative to the Step 1 ceiling: {pr['fraction_of_ceiling']['pre-specified']*100:.0f}% of the pre-specified estimate (0.67) and "
      f"{pr['fraction_of_ceiling']['bias-corrected']*100:.0f}% of the bias-corrected estimate (0.77).")

    w("\n## All models on the test set\n")
    w(table.to_markdown(floatfmt=".3f"))
    w("\n`primary_mean`: mean over patients of per-sample Spearman (samples with at least 20 drugs), patient-level bootstrap CI (10,000). "
      "`top5_hit_rate`: share of the 5 highest-predicted drugs in the observed top 20% (chance 0.20).")

    w("\n## Secondary comparisons (two-sided paired sign-flip tests over patients)\n")
    w(pd.DataFrame(sec).T.to_markdown(floatfmt=".4f"))
    w("\nHM covers only drugs with AML/monocytic-line signatures, so it is compared with H restricted to the same drugs.")

    w("\n## Amendment 1 analyses\n")
    w(f"**Low-signal samples excluded** ({out['n_test_low_signal_removed']} test samples removed; development fitting also without them):\n")
    w(pd.DataFrame(out["low_signal_excluded"]).T.to_markdown(floatfmt=".4f"))
    p4 = out["phase4_descriptive"]
    w(f"\n**Phase 4 raw reversion score (descriptive, no fitting):** mean per-patient Spearman {p4['primary_mean']:.3f} "
      f"(95% CI {p4['ci_lo']:.3f} to {p4['ci_hi']:.3f}), top-5 hit rate {p4['top5_hit_rate']:.3f}.")

    w("\n## Learned cell-type sensitivities and the venetoclax check\n")
    w("Sensitivity(k, d) = development drug mean + the primary model's prediction for a sample made only of cell type k, in sensitivity units "
      "(minus AUC; higher = more sensitive). Fitted on development data only, so this does not depend on the test set.")
    ven = out["venetoclax_sensitivity"]
    if ven:
        stem = np.mean([ven[k] for k in ("HSC-like", "Prog-like") if k in ven])
        mono = ven.get("Mono-like")
        w(f"\nVenetoclax: {', '.join(f'{k} {v:.1f}' for k, v in ven.items())}.")
        w(f"\nStem/progenitor-like (mean of HSC-like and Prog-like) {stem:.1f} vs Mono-like {mono:.1f}: "
          f"{'higher sensitivity in stem/progenitor-like cells, as expected from Pei et al. 2020' if stem > mono else 'NOT in the expected direction'}.")
    w(f"\n![test results](figures/{F['test'].name})")
    w(f"\n![cell-type sensitivity](figures/{F['heat'].name})")
    w("\n## Interpretation\n")
    w("Written after the single evaluation; see the summary in the project reports.")
    path(cfg, "reports", "exp2_results.md").write_text("\n".join(L) + "\n", encoding="utf-8")
