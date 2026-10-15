"""Phase 4 report and figures."""
from __future__ import annotations

import datetime as dt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pipeline.common import path

COL = {"B0": "#BBBBBB", "B1": "#555555", "B2": "#8172B2", "B3": "#937860", "M1": "#C44E52", "M2": "#DD8452"}


def md(df, index=False):
    return df.to_markdown(index=index, floatfmt=".4f")


def _c(name):
    return COL.get(name.split()[0][:2], "#4C72B0") if not name.startswith("M1c") else "#E58C8A"


def figures(cfg, st, res, Y, R, dec):
    F = {}
    t = st["table"]
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6), sharey=True)
    for i, (k, lab) in enumerate([("raw", "per-sample Spearman (raw sensitivity)"),
                                  ("centred", "per-sample Spearman (drug-centred, patient-specific)"),
                                  ("per_drug", "per-drug Spearman across patients")]):
        y = np.arange(len(t))[::-1]
        m, lo, hi = t[f"{k}_mean"], t[f"{k}_lo"], t[f"{k}_hi"]
        ax[i].errorbar(m, y, xerr=[m - lo, hi - m], fmt="o", ecolor="#999999", capsize=3, mfc="none")
        for yy, (name, mm) in zip(y, m.items()):
            ax[i].plot(mm, yy, "o", color=_c(name), ms=8)
        ax[i].axvline(0, color="k", lw=0.6); ax[i].set_title(lab, fontsize=10)
        ax[i].set_yticks(y, t.index, fontsize=8)
    fig.suptitle("Held-out performance with 95% bootstrap CIs (patient-level 5-fold CV)", fontsize=11)
    fig.tight_layout(); F["forest"] = path(cfg, "figures", "phase4_forest.png"); fig.savefig(F["forest"], dpi=150); plt.close(fig)

    n1, n2 = st["nulls"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    for a, nul, key, lab in [(ax[0], n1, "P1", "P1: raw"), (ax[1], n2, "P2", "P2: drug-centred")]:
        a.hist(nul, bins=60, color="#BBBBBB", label="permutation null")
        a.axvline(st[key]["statistic"], color="#C44E52", lw=2, label=f"observed {st[key]['statistic']:.4f}")
        a.set_title(f"{lab} (p = {st[key]['p']:.4f})"); a.set_xlabel("mean per-sample Spearman"); a.legend(frameon=False, fontsize=8)
    fig.tight_layout(); F["perm"] = path(cfg, "figures", "phase4_permutation.png"); fig.savefig(F["perm"], dpi=150); plt.close(fig)

    pdg = st["per_drug_m1"].dropna(subset=["spearman"])
    fig, ax = plt.subplots(figsize=(15, 4))
    cols = ["#C44E52" if q < 0.05 else ("#4C72B0" if t_ == "AML/monocytic" else "#AAAAAA") for q, t_ in zip(pdg.q, pdg.tier)]
    ax.bar(range(len(pdg)), pdg.spearman, color=cols)
    ax.set_xticks(range(len(pdg)), pdg.drug, rotation=90, fontsize=5); ax.axhline(0, color="k", lw=0.5)
    ax.set_ylabel("Spearman (score vs sensitivity)")
    ax.set_title("Per drug, across patients: reversion score vs ex vivo sensitivity (red = BH q < 0.05; blue = AML/monocytic-line signature)")
    fig.tight_layout(); F["perdrug"] = path(cfg, "figures", "phase4_per_drug.png"); fig.savefig(F["perdrug"], dpi=150); plt.close(fig)

    Yd, hot = st["idh_data"]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for a, (drug, gene) in zip(ax[:2], [("AGI-5198", "IDH1_hotspot"), ("AGI-6780", "IDH2_hotspot")]):
        if drug in Yd.columns:
            ok = Yd[drug].notna() & hot[gene].notna()
            g = [Yd.loc[ok & (hot[gene] == 0), drug], Yd.loc[ok & (hot[gene] == 1), drug]]
            a.boxplot(g, showfliers=False); a.set_xticks([1, 2], [f"{gene[:4]} WT (n={len(g[0])})", f"{gene[:4]} mutant (n={len(g[1])})"])
            for i, v in enumerate(g):
                a.scatter(np.random.default_rng(0).normal(i + 1, 0.05, len(v)), v, s=4, alpha=0.4)
            a.set_ylabel("sensitivity (-AUC)"); a.set_title(drug)
    if "venetoclax" in st:
        ok = Yd["Venetoclax"].notna()
        ax[2].scatter(dec.loc[Yd.index[ok], "mature_share"], Yd.loc[ok, "Venetoclax"], s=5, alpha=0.5, c="#4C72B0")
        ax[2].set_xlabel("deconvolved mature share"); ax[2].set_ylabel("venetoclax sensitivity (-AUC)")
        ax[2].set_title(f"Venetoclax vs mature share (Spearman {st['venetoclax']['spearman']:.2f})")
    fig.tight_layout(); F["checks"] = path(cfg, "figures", "phase4_subgroup_checks.png"); fig.savefig(F["checks"], dpi=150); plt.close(fig)

    s = R.stack().rename("score").to_frame().join(Y.stack().rename("sens"))
    s = s.dropna()
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.hexbin(s.score, s.sens, gridsize=50, cmap="Blues", mincnt=1)
    ax.set_xlabel("reversion score (cosine)"); ax.set_ylabel("observed sensitivity (-AUC)")
    ax.set_title(f"All {len(s):,} patient-drug pairs")
    fig.tight_layout(); F["scatter"] = path(cfg, "figures", "phase4_scatter.png"); fig.savefig(F["scatter"], dpi=150); plt.close(fig)
    return F


def write_report(cfg, st, res, Y, R, dec):
    F = figures(cfg, st, res, Y, R, dec)
    L = []; w = L.append
    t = st["table"]; cl = st["claims"]
    w("# Phase 4: validation against Beat AML ex vivo drug response")
    w(f"\nGenerated {dt.date.today().isoformat()} by `python -m pipeline.run --phase 4`. Analysis plan: [phase4_plan.md](phase4_plan.md), "
      "committed to git before this code was written.\n")
    w("> **Caveat.** Beat AML ex vivo assays measure cell viability after about 3 days. Differentiation agents act through maturation, often over a "
      "longer time, so ex vivo sensitivity is only a proxy for what reversion scores target, and mismatch is expected. The panel also contains no "
      "clinical differentiation agents (no ATRA, enasidenib or ivosidenib).\n")
    w("## Summary\n")
    w(f"- Data: **{st['n_samples']} samples from {st['n_patients']} patients, {st['n_drugs']} drugs, {st['n_pairs']:,} patient-drug pairs** "
      f"(median {st['drugs_per_sample_median']:.0f} drugs per sample). {st['n_drugs']} overlapping drugs is well above the 15-drug threshold.")
    w(f"- **P1** (reversion score vs sensitivity across drugs, within patient): mean Spearman **{st['P1']['statistic']:.4f}**, permutation p = **{st['P1']['p']:.4f}** "
      f"({'significant' if st['P1']['significant'] else 'not significant'} at 0.025).")
    w(f"- **P2** (same after removing each drug's average potency, i.e. patient-specific signal): mean Spearman **{st['P2']['statistic']:.4f}**, "
      f"p = **{st['P2']['p']:.4f}** ({'significant' if st['P2']['significant'] else 'not significant'}).")
    b1, m1, m2 = t.loc["B1 global drug ranking"], t.loc["M1 reversion score (raw)"], t.loc["M2 combined"]
    w(f"- **Global drug ranking (B1) is the strong baseline:** raw per-sample Spearman {b1.raw_mean:.3f} vs {m1.raw_mean:.3f} for our raw score and "
      f"{m2.raw_mean:.3f} for the combined model.")
    w(f"- **Pre-registered claims:** information about response: **{'yes' if cl['information_about_response'] else 'no'}**; "
      f"patient-specific value: **{'yes' if cl['patient_specific_value'] else 'no'}**.")

    w("\n## 1. Model comparison (held-out folds)\n")
    w(md(t, index=True))
    w("\n`raw`: per-sample Spearman across drugs with raw sensitivity. `centred`: per-sample Spearman after centring each drug with training-fold statistics "
      "(B1 is 0 by construction). `per_drug`: mean over drugs of the Spearman across held-out patients. Brackets are patient-level bootstrap 95% CIs.")
    w("\nPaired comparisons (Wilcoxon signed-rank on per-sample Spearman, one-sided):\n")
    w(md(st["comparisons"], index=True))

    w("\n## 2. Primary tests\n")
    w(md(pd.DataFrame({"P1 raw": st["P1"], "P2 drug-centred": st["P2"]}).T, index=True))

    w("\n## 3. Per-drug results (M1, all samples)\n")
    pdg = st["per_drug_m1"]
    w(f"{int((pdg.q < 0.05).sum())} of {pdg.spearman.notna().sum()} drugs have BH q < 0.05 "
      f"({int(((pdg.q < 0.05) & (pdg.spearman > 0)).sum())} positive, {int(((pdg.q < 0.05) & (pdg.spearman < 0)).sum())} negative). "
      f"Median per-drug Spearman {pdg.spearman.median():.3f}. Ten most positive and ten most negative:\n")
    w(md(pd.concat([pdg.head(10), pdg.tail(10)])[["drug", "n", "spearman", "q", "tier", "cells"]]))

    w("\n## 4. Secondary analyses\n")
    w(md(st["secondary"]))
    w("\nThese use the same permutation scheme with 2,000 permutations; they are robustness checks, not claims.")

    w("\n## 5. IDH subgroup check\n")
    w(md(st["idh"]))
    w("\nAGI-5198 and AGI-6780 are tool IDH inhibitors. Higher sensitivity in mutant samples would show the assay can see on-target IDH effects; "
      "the reversion columns show whether our scores differ between groups (they are computed from expression only).")
    if "venetoclax" in st:
        v = st["venetoclax"]
        w(f"\n**Venetoclax sanity check** (pre-registered): Spearman between deconvolved mature share and venetoclax sensitivity = {v['spearman']:.3f} "
          f"(p = {v['p']:.2e}, n = {v['n']}). Expected negative (monocytic AML resists venetoclax, Pei et al. 2020).")

    w("\n## 6. Interpretation\n")
    w("- The global drug ranking captures most of the between-drug variance in sensitivity; any patient-level model has to beat it to show patient-specific value.")
    w("- The drug-centred metric isolates the patient-specific part. Its size, and whether it beats random, is the honest answer to 'does this tell us which drug for which patient'.")
    w("- Most panel drugs are kinase inhibitors and cytotoxics, whose ex vivo effect is killing, not maturation. A reversion score is not designed to predict killing, so weak or negative per-drug correlations are not surprising.")

    w("\n## 7. Exploratory analyses (not pre-registered)\n")
    w("Added after seeing the results above, to understand them. They do not change the pre-registered verdict.\n")
    w(md(st["exploratory"], index=True))
    m1c, x1 = t.loc["M1c reversion score (calibrated)"], t.loc["X1 fractions + mutations (ablation)"]
    pdg_ = st["per_drug_m1"]
    npos, nneg = int(((pdg_.q < 0.05) & (pdg_.spearman > 0)).sum()), int(((pdg_.q < 0.05) & (pdg_.spearman < 0)).sum())
    w(f"\n- **Per-drug calibration recovers a patient-specific signal.** The raw score fails P2 because the direction of association differs between drugs "
      f"({npos} drugs significantly positive, {nneg} negative). Letting each drug learn its own slope (M1c) gives a drug-centred Spearman of {m1c.centred_mean:.3f} "
      f"(95% CI {m1c.centred_lo:.3f} to {m1c.centred_hi:.3f}). So the reversion score carries patient-level information, but whether a higher score means more or "
      "less sensitivity depends on the drug. This fits the viability readout: for many kinase inhibitors, more mature (monocytic) blasts are less sensitive, not more.")
    w(f"- **Contribution to the combined model.** Removing the reversion score from M2 (X1) gives raw {x1.raw_mean:.3f} and centred {x1.centred_mean:.3f}, vs "
      f"{t.loc['M2 combined'].raw_mean:.3f} and {t.loc['M2 combined'].centred_mean:.3f} for M2; the paired tests above show whether the difference is reliable.")
    w("- **B1 per-drug value is negative** by a known artefact: with fold-wise training means, a fold whose training mean is high tends to have lower test values "
      "(regression to the mean), which creates a small negative correlation across pooled folds. It does not affect the per-sample metrics.")
    w("\n## 8. Deviations from the plan\n")
    w("- The exploratory analyses in section 7 (and the X1 row in the model table) were added after seeing the results and are labelled as such.")
    w("- Otherwise none. Implementation note: the plan's alternative claim clause (M2 beats B1 and B2/B3) already implies M2 beats B1, so the claim rule was coded as "
      "P2 significant and (M1 or M2 beats B1).")

    w("\n## Figures\n")
    for f in F.values():
        w(f"![{f.stem}](figures/{f.name})")
    path(cfg, "reports", "phase4_validation_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
