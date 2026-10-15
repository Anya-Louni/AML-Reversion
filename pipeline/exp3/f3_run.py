"""Experiment 3, F3 runner. Plan: reports/exp3_f3_plan.md. Writes reports/exp3_f3_results.md."""
from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

import anndata as ad
import joblib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pipeline.common import get_logger, path
from pipeline.exp3 import f3 as F3

log = get_logger("exp3.f3")
PDX = {"control": "GSM7118799_DF68555_Control", "treated": "GSM7118800_DF68555_500VTP"}
VG_PATIENTS = ["AML328", "AML329", "AML707B"]


def arm_pdx(cfg, proc, mapper, ref_Z, raw: Path) -> pd.DataFrame:
    """PDX AML cells, control vs VTP50469 at 16 h. One PDX, so cell-level folds only."""
    frames = []
    for label, prefix in PDX.items():
        X, cells, genes = F3.load_10x(raw, prefix)
        X, genes = F3.collapse_duplicate_genes(X, genes)
        d = F3.place(mapper, X, genes, ref_Z).join(F3.marker_frame(X, genes))
        d["treated"] = label == "treated"
        d["sample"] = label
        frames.append(d)
    df = pd.concat(frames, ignore_index=True)
    # cell type is not annotated in this dataset; assign by nearest reference cell type on the landscape
    df["cell_type"] = assign_type(df, mapper, proc)
    return df


def assign_type(df: pd.DataFrame, mapper, proc) -> np.ndarray:
    """Nearest-reference-cell type for each mapped cell, using landscape coordinates."""
    from sklearn.neighbors import KNeighborsClassifier
    ref = ad.read_h5ad(proc / "landscape_reference.h5ad")
    k = KNeighborsClassifier(15).fit(ref.obsm["X_landscape"], ref.obs.CellType.astype(str).values)
    return k.predict(df[["x", "y"]].values)


def arm_vangalen(cfg, proc, mapper, ref_Z) -> pd.DataFrame:
    """Malignant cells at day 0 vs follow-up, for the 3 patients with enough of both."""
    vg = ad.read_h5ad(proc / "vangalen_raw.h5ad")
    o = vg.obs
    keep = (o.sample_type == "AML") & (o.PredictionRefined == "malignant") & o.patient.isin(VG_PATIENTS) \
        & o.CellType.isin(cfg["phase2"]["malignant_like_types"])
    sub = vg[keep.values]
    counts = sub.X
    d = F3.place(mapper, counts, sub.var_names, ref_Z).join(F3.marker_frame(counts, sub.var_names))
    d["patient"] = sub.obs.patient.astype(str).values
    d["day"] = sub.obs.day.values
    d["cell_type"] = sub.obs.CellType.astype(str).values
    d["treated"] = d.day > 0
    # keep only timepoints with enough malignant cells on both sides, as pre-registered
    ok = []
    for p, g in d.groupby("patient"):
        if (~g.treated).sum() >= 50 and g.treated.sum() >= 50:
            ok.append(g)
    return pd.concat(ok, ignore_index=True)


def run(cfg: dict) -> None:
    proc = path(cfg, "processed", "x").parent
    raw = Path("data/raw/exp3")
    mapper = joblib.load(proc / "landscape_mapper.joblib")
    ref = ad.read_h5ad(proc / "landscape_reference.h5ad")
    ref_Z = mapper.embed(mapper.align(ref.X, ref.var_names)[0])
    thr = json.loads((proc / "phase2_stats.json").read_text())["threshold"]["threshold"]

    arms = {}
    log.info("arm A: PDX menin inhibitor")
    A = arm_pdx(cfg, proc, mapper, ref_Z, raw)
    A_folds = np.array([i % 5 for i in range(len(A))])
    rng = np.random.default_rng(cfg["seed"]); rng.shuffle(A_folds)
    arms["A PDX, menin inhibitor 16 h"] = (A, A_folds, "cells (5-fold); 1 PDX, no held-out sample possible")
    log.info("arm B: van Galen chemotherapy")
    B = arm_vangalen(cfg, proc, mapper, ref_Z)
    B_folds = pd.factorize(B.patient)[0]
    arms["B van Galen, chemotherapy"] = (B, B_folds, f"leave-one-patient-out ({B.patient.nunique()} patients)")

    rows, direction, importances = [], {}, {}
    for name, (df, folds, split) in arms.items():
        for mid, feats in F3.FEATURE_SETS.items():
            r = F3.fit_eval(df, feats, folds, cfg["seed"])
            r.pop("pred")
            rows.append({"arm": name, "split": split, "model": mid, **r})
            log.info("%s | %s: AUROC %.3f", name, mid, r["auroc"])
        direction[name] = F3.maturity_direction(df, thr)
        importances[name] = F3.permutation_importance_auroc(df, F3.FEATURE_SETS["F full state"], folds, 3, cfg["seed"])

    res = pd.DataFrame(rows)
    piv = res.pivot(index="arm", columns="model", values="auroc")
    verdict = {}
    for arm in piv.index:
        full = piv.loc[arm, "F full state"]
        verdict[arm] = {"full": float(full),
                        "over_S1": float(full - piv.loc[arm, "S1 maturity score only"]),
                        "over_S2": float(full - piv.loc[arm, "S2 cell type only"])}
        verdict[arm]["beats_baselines"] = bool(verdict[arm]["over_S1"] > 0.02 and verdict[arm]["over_S2"] > 0.02)
    overall = all(v["beats_baselines"] for v in verdict.values())

    # link to our own simulation, for the van Galen patients
    sim_link = simulation_link(proc, B, thr)
    st = {"results": res.to_dict(orient="records"), "verdict": verdict, "overall": overall,
          "direction": direction, "threshold": thr,
          "importances": {k: v.round(4).to_dict() for k, v in importances.items()},
          "simulation_link": sim_link,
          "cells": {k: {"untreated": int((~v[0].treated).sum()), "treated": int(v[0].treated.sum())} for k, v in arms.items()}}
    res.to_csv(proc / "exp3_f3_results.csv", index=False)
    (proc / "exp3_f3_stats.json").write_text(json.dumps(st, indent=2, default=str))
    figure(cfg, arms, piv, direction, thr)
    write(cfg, st, res, piv, verdict, overall, direction, importances, sim_link)
    log.info("F3 done: beats baselines in both arms = %s", overall)


def simulation_link(proc, B: pd.DataFrame, thr: float) -> dict:
    """Does Phase 5's simulated resistant fraction track measured persistence? Weak by construction."""
    f = proc / "phase5_simulation_summary.csv"
    if not f.exists():
        return {"available": False}
    sim = pd.read_csv(f)
    sim = sim[sim.patient.isin(B.patient.unique())]
    if sim.empty:
        return {"available": False}
    simf = sim.groupby("patient").apply(lambda g: (g.resistant / g.n).mean(), include_groups=False)
    meas = B[B.treated].groupby("patient").apply(lambda g: float((g.score < thr).mean()), include_groups=False)
    both = pd.DataFrame({"simulated_resistant_fraction": simf, "measured_stuck_fraction_after": meas}).dropna()
    return {"available": True, "n_patients": len(both), "table": both.round(3).to_dict(orient="index"),
            "note": "Phase 5 simulates differentiation agents; these patients received chemotherapy, so no relationship is expected."}


def figure(cfg, arms, piv, direction, thr):
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.3))
    P = piv[["S1 maturity score only", "S2 cell type only", "F full state"]]
    x = np.arange(len(P)); w = 0.26
    for i, c in enumerate(P.columns):
        ax[0].bar(x + i * w - w, P[c], w, label=c)
    ax[0].axhline(0.5, color="k", ls="--", lw=0.8)
    ax[0].set_xticks(x, [i.split(",")[0] for i in P.index], fontsize=8)
    ax[0].set_ylabel("held-out AUROC"); ax[0].set_title("Separating treated from untreated cells")
    ax[0].legend(fontsize=7, frameon=False)
    for a, (name, (df, _, _)) in zip(ax[1:], arms.items()):
        for lab, sel in [("untreated", ~df.treated.astype(bool)), ("treated", df.treated.astype(bool))]:
            a.hist(df.score[sel], bins=40, alpha=0.6, density=True, label=lab)
        a.axvline(thr, color="#C44E52", ls="--", lw=1, label="maturity threshold")
        a.set_xlabel("maturity score"); a.set_title(name.split(",")[0] + ": where cells sit", fontsize=9)
        a.legend(fontsize=7, frameon=False)
    fig.tight_layout(); f = path(cfg, "figures", "exp3_f3_resistance.png"); fig.savefig(f, dpi=150); plt.close(fig)


def write(cfg, st, res, piv, verdict, overall, direction, importances, sim_link):
    L = []; w = L.append
    w("# Experiment 3, F3: can we predict which cells stay stuck?")
    w(f"\nGenerated {dt.date.today().isoformat()} by `python -m pipeline.run --phase exp3-f3`. "
      "Plan: [exp3_f3_plan.md](exp3_f3_plan.md), committed before any model was fitted. "
      "Constraints: [exp3_f3_feasibility.md](exp3_f3_feasibility.md).\n")
    w("## Read this first\n")
    w("**1. The question as posed cannot be answered.** Single-cell sequencing destroys the cell it measures, so no cell is observed both before and "
      "after treatment. Whether a *given* cell resists is not measurable in any open dataset. We substitute a state-level question: which untreated "
      "cell states are over-represented after treatment. Everything below is about states, not cells, and the whole of F3 is exploratory.")
    w("\n**2. Treatment and batch are confounded.** Each arm compares separate sequencing runs, so a classifier that separates treated from untreated "
      "cells is partly detecting batch. The AUROC numbers below therefore **do not support a biological claim on their own**. The maturity-direction "
      "result further down is the more trustworthy readout, because it is a prediction about one specific axis rather than about separability.")
    w(f"\n**3. The data is thin:** {st['cells']['A PDX, menin inhibitor 16 h']['untreated']} and "
      f"{st['cells']['A PDX, menin inhibitor 16 h']['treated']} cells from **one PDX and one drug**, plus **three patients** on chemotherapy.")

    w("\n## Classifier results\n")
    w(res[["arm", "split", "model", "auroc", "auroc_per_fold_mean", "n_cells"]].to_markdown(index=False, floatfmt=".4f"))
    w("\n### Against the pre-registered bar\n")
    w(pd.DataFrame(verdict).T.to_markdown(floatfmt=".4f"))
    w(f"\nPre-registered rule: the full state model must beat **both** baselines by more than 0.02 AUROC in **both** arms. "
      f"**Result: {'met' if overall else 'not met'}.**")

    w("\n## Maturity direction (the more trustworthy readout)\n")
    w(pd.DataFrame(direction).T.to_markdown(floatfmt=".4f"))
    w(f"\n`stuck_frac` is the share of cells below the Phase 2 maturity threshold ({st['threshold']:.3f}), which was fixed on normal marrow long before "
      "Experiment 3. A differentiation agent should lower it. A rise means the cells that persist are the less mature ones, which is the resistance "
      "pattern this project predicts.")

    w("\n## What marks the persisting cells\n")
    for name, imp in importances.items():
        w(f"\n**{name}** (drop in held-out AUROC when each feature is shuffled):\n")
        w(imp.to_frame("auroc_drop").to_markdown(floatfmt=".4f"))

    w("\n## Link to our own simulation\n")
    if sim_link.get("available"):
        w(pd.DataFrame(sim_link["table"]).T.to_markdown(floatfmt=".3f"))
        w(f"\n{sim_link['note']} With {sim_link['n_patients']} patients no correlation is computable; the numbers are shown for transparency only.")
    else:
        w("Not available.")
    w("\n## Interpretation\n")
    a_key, b_key = "A PDX, menin inhibitor 16 h", "B van Galen, chemotherapy"
    va, vb = verdict[a_key], verdict[b_key]
    da, db = direction[a_key], direction[b_key]
    rb = res[(res.arm == b_key)].set_index("model")
    w(f"**The pre-registered bar was not met.** Full cell state beat both baselines in the PDX arm ({va['full']:.3f}, "
      f"+{va['over_S1']:.3f} over maturity alone and +{va['over_S2']:.3f} over cell type alone) but not in the patient arm.")
    w(f"\n**1. Nothing transfers between patients.** Every van Galen model scores *below* chance when predictions are pooled across patients "
      f"({rb.loc['S1 maturity score only', 'auroc']:.2f}, {rb.loc['S2 cell type only', 'auroc']:.2f} and {rb.loc['F full state', 'auroc']:.2f}), which "
      "means a rule learned on two patients is actively misleading on the third. Scoring each held-out patient separately gives a different picture "
      f"({rb.loc['S1 maturity score only', 'auroc_per_fold_mean']:.2f}, {rb.loc['S2 cell type only', 'auroc_per_fold_mean']:.2f} and "
      f"{rb.loc['F full state', 'auroc_per_fold_mean']:.2f}), and cell type is the clearest case: {rb.loc['S2 cell type only', 'auroc_per_fold_mean']:.2f} "
      f"within patients against {rb.loc['S2 cell type only', 'auroc']:.2f} pooled. The gap is the signature of an effect whose direction differs from patient "
      "to patient, so the pooled ranking is scrambled. With three patients we cannot say whether that is real heterogeneity or noise; either way nothing transfers.")
    w(f"\n**2. The PDX result cannot be read as biology.** AUROC {va['full']:.3f} looks strong, but control and treated cells are separate sequencing "
      "runs from one animal, so batch and treatment are perfectly confounded. A classifier that good is exactly what a batch effect would produce. "
      "This number is reported for completeness and should not be cited as evidence that resistant states are predictable.")
    w(f"\n**3. The one encouraging result is the maturity shift.** Under the menin inhibitor, mean maturity rises by {da['delta_mean_score']:+.3f} "
      "— the direction a differentiation-inducing drug should produce, on an axis fixed long before this experiment. It is small, and the batch "
      f"confound applies here too, but a pure batch effect has no reason to favour this particular direction. The share of cells below the maturity "
      f"threshold barely moves ({da['delta_stuck_frac']:+.3f}), so at 16 h the population shifts without many cells crossing over.")
    w(f"\n**4. Chemotherapy does not mature cells.** In the patient arm, mean maturity changes by {db['delta_mean_score']:+.3f} and the stuck fraction by "
      f"{db['delta_stuck_frac']:+.3f}, both negligible. The malignant cells that survive are not more mature, they are simply fewer. This matches what "
      "F0 predicted from the cell counts and confirms that these samples cannot stand in for differentiation therapy.")
    w("\n**5. Features.** In the PDX arm the strongest single feature is PROM1 (CD133), a stem/progenitor marker, followed by LYZ, a monocytic one, which "
      "is at least consistent with a stem-to-mature axis separating the two populations. In the patient arm almost every importance is negative — "
      "shuffling a feature *improves* held-out AUROC — which is what one expects when there is no transferable signal to destroy.")
    w("\n**Against published biology.** Resistance to menin inhibitors is well documented to involve MEN1 mutations at the drug-binding interface "
      "(PubMed lists dozens of reports). **We cannot test that here:** these data are expression only, with no genotyping of the treated cells, and 16 h "
      "is far too early for an outgrowth of resistant clones. Any claim that our features correspond to that mechanism would be unsupported, so we make none.")
    w("\n**Bottom line.** F3 does not demonstrate that resistant cell states are predictable. It establishes that the question cannot be asked properly of "
      "open data, that the one dataset able to address it is confounded, and that patient-level signal does not transfer across three patients. The single "
      "result worth carrying forward is the small, correctly directed maturity shift under a menin inhibitor.")

    w("\n## Figure\n")
    w("![resistance](figures/exp3_f3_resistance.png)")
    path(cfg, "reports", "exp3_f3_results.md").write_text("\n".join(L) + "\n", encoding="utf-8")
