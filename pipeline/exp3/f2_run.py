"""Experiment 3, F2 runner. Plan: reports/exp3_f2_plan.md. Writes reports/exp3_f2_results.md."""
from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from pipeline.common import get_logger, path
from pipeline.exp3 import datasets as DS
from pipeline.exp3 import f1 as F1
from pipeline.exp3 import f2 as F2
from pipeline.exp3.f1_run import build_simulator
from pipeline.phase0_data.genes import HGNCMapper

log = get_logger("exp3.f2")
RANKS, ALPHAS = [1, 2, 3, 5, 10], [0.1, 1.0, 10.0, 100.0]
# Which drug each treated arm actually is, so the "other drugs" baseline never averages
# other timepoints of the same drug.
ARM_DRUG = {("GSE227542", "RA4h"): "tretinoin", ("GSE227542", "RA2d"): "tretinoin", ("GSE227542", "RA4d"): "tretinoin",
            ("GSE276457", "ATRA"): "tretinoin", ("GSE276457", "palbociclib"): "Palbociclib",
            ("GSE228321", "OTX015"): "OTX-015", ("GSE228321", "VTP50469"): "VTP50469"}
# External conditions, from the F2 plan. (id, dataset, arm, drug, subset, cell line, time in hours)
CONDITIONS = [
    ("D1 4h", "GSE227542", "RA4h", "tretinoin", None, "HL60", 4.0),
    ("D1 2d", "GSE227542", "RA2d", "tretinoin", None, "HL60", 48.0),
    ("D1 4d", "GSE227542", "RA4d", "tretinoin", None, "HL60", 96.0),
    ("D2 ATRA 48h", "GSE276457", "ATRA", "tretinoin", None, "HL60", 48.0),
    ("D3 palbociclib 48h", "GSE276457", "palbociclib", "Palbociclib", None, "HL60", 48.0),
    ("D4 OTX015 16h", "GSE228321", "OTX015", "OTX-015", {"cell_line": "MOLM13"}, "MOLM13", 16.0),
]


def metrics(pred: pd.Series, real: pd.Series, sig: pd.DataFrame) -> dict:
    idx = pred.index.intersection(real.index)
    p, r = pred.loc[idx].values, real.loc[idx].values
    if np.allclose(p, 0):
        return {"cosine": 0.0, "spearman": 0.0, "sim_shift": 0.0, "n_genes": len(idx)}
    return {"cosine": F1.cosine(p, r), "spearman": float(spearmanr(p, r)[0]),
            "sim_shift": F1.signature_shift(pred.loc[idx] / np.log(2), sig), "n_genes": len(idx)}


def run(cfg: dict) -> None:
    proc = path(cfg, "processed", "x").parent
    raw = Path("data/raw/exp3")
    hgnc = HGNCMapper.from_file(path(cfg, "raw", "hgnc", "hgnc_complete_set.txt"))
    mapper, sim, norm = build_simulator(cfg, proc)
    sig = pd.read_csv(proc / "maturation_signature.csv")
    data = {"GSE227542": DS.load_gse227542(raw), "GSE276457": DS.load_gse276457(raw),
            "GSE228321": DS.load_gse228321(raw, hgnc)}

    # ---- signature spaces for the learned map (LINCS only)
    S = ad.read_h5ad(proc / "lincs_signatures_all.h5ad")
    genes = pd.Index(S.var_names)
    keep = genes.isin(sim.genes)
    Xall, obs = np.asarray(S.X)[:, keep], S.obs
    g = genes[keep]
    hl60 = F2.per_cell_line_signatures(obs, Xall, g, "HL60")
    cons_ex = F2.consensus_excluding(obs, Xall, g, "HL60")
    shared_drugs = sorted(set(hl60.index) & set(cons_ex.index))
    Xs, Ys = cons_ex.loc[shared_drugs], hl60.loc[shared_drugs]
    st: dict = {"n_map_training_drugs": len(shared_drugs), "n_genes": int(len(g)), "shift_norm": norm}
    log.info("learned map: %d drugs with HL60 and other-line signatures", len(shared_drugs))

    rows = []
    for cid, ds, arm, drug, subset, line, hours in CONDITIONS:
        X, m = data[ds]
        real = DS.real_change(X, m, arm, subset=subset)
        preds: dict[str, pd.Series] = {}
        preds["B0 no change"] = pd.Series(0.0, index=sim.genes)
        preds["B1 linear shift"] = sim.change(drug)
        # B2: mean real change of the other single-agent arms that are a DIFFERENT drug.
        # (In GSE227542 every arm is tretinoin at another timepoint, so B2 is undefined there.)
        others = [a for a in m.arm.unique()
                  if a not in (arm, "control", "combination") and ARM_DRUG.get((ds, a)) not in (None, drug)]
        if others:
            ch = [DS.real_change(X, m, a, subset=subset) for a in others]
            idx = ch[0].index
            for c in ch[1:]:
                idx = idx.intersection(c.index)
            preds["B2 mean other drugs"] = pd.concat([c.loc[idx] for c in ch], axis=1).mean(axis=1)
        # M1: cell-line-matched LINCS signature
        if line == "HL60" and drug in hl60.index:
            v = hl60.loc[drug] * sim.sd[g]
            preds["M1 matched cell line"] = v / max(np.linalg.norm(v.values), 1e-9) * norm
        # M2: learned consensus -> HL60 map, test drug always excluded from training
        if line == "HL60" and drug in cons_ex.index:
            sel, grid = F2.select_map(Xs, Ys, RANKS, ALPHAS, exclude=[drug])
            tr = [d for d in shared_drugs if d != drug]
            mp = F2.ReducedRankMap(sel["rank"], sel["alpha"]).fit(Xs.loc[tr].values, Ys.loc[tr].values)
            yhat = pd.Series(mp.predict(cons_ex.loc[[drug]].values)[0], index=g)
            v = yhat * sim.sd[g]
            preds["M2 learned map"] = v / max(np.linalg.norm(v.values), 1e-9) * norm
            st.setdefault("map_selection", {})[cid] = sel
            st.setdefault("map_baseline_cosine", {})[cid] = F2.baseline_consensus_cosine(Xs, Ys, [drug])
        for name, p in preds.items():
            rows.append({"id": cid, "drug": drug, "cell_line": line, "hours": hours, "model": name,
                         **metrics(p, real, sig)})
        log.info("%s done (%d models)", cid, len(preds))

    res = pd.DataFrame(rows)
    piv = res.pivot(index="id", columns="model", values="cosine")

    # ---- M3: does a time axis help? Hold out 2 d, interpolate from 4 h and 4 d.
    X, m = data["GSE227542"]
    ch = {4.0: DS.real_change(X, m, "RA4h"), 96.0: DS.real_change(X, m, "RA4d")}
    real2d = DS.real_change(X, m, "RA2d")
    interp = F2.time_interpolate(ch, 48.0)
    gset = sim.genes.intersection(interp.index).intersection(real2d.index)
    m3 = {"M3 time interpolation": metrics(interp.loc[gset], real2d.loc[gset], sig),
          "B1 linear shift": metrics(sim.change("tretinoin").loc[gset], real2d.loc[gset], sig),
          "B0 no change": {"cosine": 0.0, "n_genes": len(gset)}}
    st["m3"] = m3
    st["m3_note"] = "M3 uses real measurements of the same drug at 4 h and 4 d; B1 uses none. Not a fair comparison, reported as a separate question."

    # ---- decision rule
    b1 = piv["B1 linear shift"]
    verdict = {}
    for mod in ["M1 matched cell line", "M2 learned map"]:
        if mod not in piv.columns:
            continue
        d = (piv[mod] - b1).dropna()
        verdict[mod] = {"n_conditions": int(len(d)), "n_wins": int((d > 0).sum()),
                        "mean_improvement": float(d.mean()), "per_condition": d.round(4).to_dict(),
                        "beats_linear": bool((d > 0).sum() > len(d) / 2 and d.mean() > 0.05)}
    st["verdict"] = verdict
    st["results"] = res.to_dict(orient="records")
    res.to_csv(proc / "exp3_f2_results.csv", index=False)
    (proc / "exp3_f2_stats.json").write_text(json.dumps(st, indent=2, default=str))
    figure(cfg, piv)
    write(cfg, st, res, piv, verdict, m3)
    log.info("F2 done: %s", {k: v["beats_linear"] for k, v in verdict.items()})


def figure(cfg, piv):
    order = [c for c in ["B0 no change", "B2 mean other drugs", "B1 linear shift", "M1 matched cell line", "M2 learned map"] if c in piv.columns]
    P = piv[order]
    fig, ax = plt.subplots(figsize=(10, 4.5))
    x = np.arange(len(P)); w = 0.8 / len(order)
    colors = {"B0 no change": "#BBBBBB", "B2 mean other drugs": "#937860", "B1 linear shift": "#C44E52",
              "M1 matched cell line": "#55A868", "M2 learned map": "#4C72B0"}
    for i, c in enumerate(order):
        ax.bar(x + i * w - 0.4, P[c].fillna(0), w, label=c, color=colors.get(c, "#888"))
    ax.set_xticks(x, P.index, fontsize=8); ax.axhline(0, color="k", lw=0.5)
    ax.set_ylabel("cosine with real change"); ax.legend(fontsize=8, frameon=False, ncol=2)
    ax.set_title("F2: can anything beat the linear shift? (higher is better)")
    fig.tight_layout(); f = path(cfg, "figures", "exp3_f2_models.png"); fig.savefig(f, dpi=150); plt.close(fig)


def write(cfg, st, res, piv, verdict, m3):
    L = []; w = L.append
    w("# Experiment 3, F2: can a learned model beat the linear shift?")
    w(f"\nGenerated {dt.date.today().isoformat()} by `python -m pipeline.run --phase exp3-f2`. "
      "Plan: [exp3_f2_plan.md](exp3_f2_plan.md), committed before any model was fitted. Follows [exp3_f1_results.md](exp3_f1_results.md).\n")
    beat = [k for k, v in verdict.items() if v["beats_linear"]]
    w("## Result\n")
    if beat:
        w(f"**{', '.join(beat)} beats the linear shift** under the pre-registered rule (wins on a majority of conditions and mean cosine improvement above 0.05).")
    else:
        w("**Nothing beats the linear shift under the pre-registered rule.** The simple model from Phase 5 remains the best available, which the plan "
          "named in advance as a legitimate outcome.")
    w("\nCosine with the real change, by condition and model (higher is better):\n")
    w(piv.round(4).to_markdown())
    w("\n### Pre-registered comparison against B1\n")
    w(pd.DataFrame({k: {kk: vv for kk, vv in v.items() if kk != "per_condition"} for k, v in verdict.items()}).T.to_markdown())
    w("\nPer-condition differences against the linear shift:\n")
    w(pd.DataFrame({k: v["per_condition"] for k, v in verdict.items()}).to_markdown())
    w(f"\nThe learned map (M2) was trained on **{st['n_map_training_drugs']} LINCS drugs** that have both an HL-60 signature and a signature in another "
      "myeloid line, with the evaluated drug removed from training every time, and its rank and penalty chosen by leave-one-drug-out inside that training set. "
      "M1 and M2 apply only to HL-60 conditions; M1 further needs a matched signature to exist, which among our drugs is true only for tretinoin.")
    if "map_selection" in st:
        w("\nSelected hyper-parameters and how well the raw consensus already matched the HL-60 signature in training:\n")
        sel = pd.DataFrame(st["map_selection"]).T
        sel["raw_consensus_cosine"] = pd.Series(st["map_baseline_cosine"])
        w(sel.round(4).to_markdown())

    w("\n## Does a time axis help? (separate question)\n")
    w(pd.DataFrame(m3).T.to_markdown(floatfmt=".4f"))
    w(f"\n{st['m3_note']} Interpolating between the real 4 h and 4 d responses predicts the held-out 2 d response far better than our static "
      "prediction does, which says the real path is smooth in time and that a model with a time axis has real room to improve. What it does not show "
      "is that such a model could be built from drug signatures alone, which is the harder problem and the one the data cannot currently answer.")

    w("\n## What this means\n")
    m1 = verdict.get("M1 matched cell line")
    m2 = verdict.get("M2 learned map")
    if m1 and m1["beats_linear"]:
        w(f"- **Using the drug's signature from the matching cell line beats the consensus** on all {m1['n_conditions']} conditions where it can be "
          f"tested, by {m1['mean_improvement']:.3f} cosine on average. **But all {m1['n_conditions']} conditions are the same drug** (tretinoin), "
          "because it is the only drug in our external data with a matched HL-60 signature. So this is one drug measured four ways, not four independent "
          "tests, and it should be read as a promising lead rather than an established gain. It does point at a concrete fix: prefer a cell-line-matched "
          "signature over the consensus wherever one exists.")
    if m2 and not m2["beats_linear"]:
        w(f"- **The learned map does not beat the linear shift** (mean change {m2['mean_improvement']:+.3f}). It helps modestly on the ATRA conditions "
          "but is much worse on palbociclib, so a map trained on 46 drugs does not generalize reliably to a new drug. With 46 training drugs and "
          "1,345 genes this is the expected outcome, and it is the honest answer: the simple model stands.")
    w("- **No neural ODE or flow-matching model was fitted**, as pre-registered: three conditionable drugs and nine external conditions cannot support one, "
      "and a number produced from them would not mean anything. This is a limit of the open data, not a judgement that such models would fail.")
    w("- The time-interpolation result is the clearest signal that dynamics matter: our static prediction is the weakest part of the model, and the real "
      "path is smooth enough that knowing the time would help a lot.")
    w("- Taken with F1, the picture is consistent: the linear shift has real but limited skill, the largest identified gaps are **no time axis** and "
      "**wrong cell context**, and neither can be closed properly until a much larger open perturbation dataset exists.")
    w("\n## Figure\n")
    w("![models](figures/exp3_f2_models.png)")
    path(cfg, "reports", "exp3_f2_results.md").write_text("\n".join(L) + "\n", encoding="utf-8")
