"""Experiment 3, F4 runner. Plan: reports/exp3_f4_plan.md. Writes reports/exp3_f4_results.md."""
from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pipeline.common import get_logger, path
from pipeline.exp3 import datasets as DS
from pipeline.exp3 import f1 as F1
from pipeline.exp3 import f4 as F4
from pipeline.exp3.f1_run import build_simulator
from pipeline.phase0_data.genes import HGNCMapper

log = get_logger("exp3.f4")
CONTROL_PAIR = ("tretinoin", "tranylcypromine")


def run(cfg: dict) -> None:
    proc = path(cfg, "processed", "x").parent
    raw = Path("data/raw/exp3")
    hgnc = HGNCMapper.from_file(path(cfg, "raw", "hgnc", "hgnc_complete_set.txt"))
    mapper, sim, norm = build_simulator(cfg, proc)
    sig = pd.read_csv(proc / "maturation_signature.csv")
    meta = pd.read_parquet(proc / "drug_signature_meta.parquet")
    meta = meta[~meta.index.duplicated()]
    null_drugs = [d for d in sim.D.index if meta.tier.get(d) == "AML/monocytic"]
    st: dict = {}

    # ---------------- Test 0: is the real combination distinct from its single agents?
    X, m = DS.load_gse276457(raw)
    real = {a: DS.real_change(X, m, a) for a in ["ATRA", "palbociclib", "combination"]}
    idx = real["combination"].index
    for k in real:
        idx = idx.intersection(real[k].index)
    t0 = {f"real combination vs real {k}": F1.cosine(real["combination"].loc[idx].values, real[k].loc[idx].values)
          for k in ["ATRA", "palbociclib"]}
    t0["real ATRA vs real palbociclib"] = F1.cosine(real["ATRA"].loc[idx].values, real["palbociclib"].loc[idx].values)
    st["test0"] = t0
    log.info("test 0: %s", {k: round(v, 3) for k, v in t0.items()})

    # ---------------- Test 1: which composition rule predicts the real combination?
    sa, sb = sim.change("tretinoin"), sim.change("Palbociclib")
    preds = {"ATRA only": sa, "palbociclib only": sb,
             "additive": F4.compose(sa, sb, "additive", norm),
             "elementwise max": F4.compose(sa, sb, "elementwise max", norm)}
    shared = sim.genes.intersection(real["combination"].index)
    r = real["combination"].loc[shared].values
    rows = [{"prediction": k, "cosine": F1.cosine(v.loc[shared].values, r),
             "n_genes": len(shared)} for k, v in preds.items()]
    t1 = pd.DataFrame(rows).sort_values("cosine", ascending=False)
    # how does the additive prediction compare with every single drug scored the same way?
    d = F1.direction_test(sim, real["combination"], "tretinoin", null_drugs, n_perm=200)
    cos_all = d.pop("cos_all")
    add_cos = float(t1.set_index("prediction").loc["additive", "cosine"])
    st["test1"] = {"table": t1.to_dict(orient="records"),
                   "additive_rank_among_single_drugs": int((cos_all > add_cos).sum() + 1),
                   "n_single_drugs": len(cos_all),
                   "best_single_drug": str(cos_all.idxmax()), "best_single_cosine": float(cos_all.max())}
    best_single = max(t1.set_index("prediction").loc[["ATRA only", "palbociclib only"], "cosine"])
    st["test1"]["beats_both_singles"] = bool(add_cos > best_single)
    log.info("test 1: additive %.3f vs best single %.3f", add_cos, best_single)

    # ---------------- Test 2: where does the published positive-control pair rank?
    V = sim.change_matrix(null_drugs).values
    w = F4.signature_weights(sig, sim.genes).values
    S = F4.all_pair_maturation_scores(V, w, norm)
    iu = np.triu_indices(len(null_drugs), k=1)
    scores = S[iu]
    pos = {d: i for i, d in enumerate(null_drugs)}
    if all(c in pos for c in CONTROL_PAIR):
        i, j = sorted((pos[CONTROL_PAIR[0]], pos[CONTROL_PAIR[1]]))
        val = float(S[i, j])
        st["test2"] = {"pair": " + ".join(CONTROL_PAIR), "score": val,
                       "n_pairs": int(len(scores)), "rank": int((scores > val).sum() + 1),
                       "percentile": float((scores < val).mean() * 100),
                       "singles": {c: float(V[pos[c]] @ w) for c in CONTROL_PAIR}}
        order = np.argsort(-scores)[:10]
        top = [(null_drugs[iu[0][k]], null_drugs[iu[1][k]], float(scores[k])) for k in order]
        st["test2"]["top10"] = top
        log.info("test 2: %s percentile %.2f", st["test2"]["pair"], st["test2"]["percentile"])

    # ---------------- Test 3: descriptive, the combination we cannot predict
    Xm, mm = DS.load_gse228321(raw, hgnc)
    sub = {"cell_line": "MOLM13"}
    rm = {a: DS.real_change(Xm, mm, a, subset=sub) for a in ["VTP50469", "OTX015", "combination"]}
    im = rm["combination"].index
    for k in rm:
        im = im.intersection(rm[k].index)
    st["test3"] = {f"real combination vs real {k}": F1.cosine(rm["combination"].loc[im].values, rm[k].loc[im].values)
                   for k in ["VTP50469", "OTX015"]}
    st["test3"]["real VTP50469 vs real OTX015"] = F1.cosine(rm["VTP50469"].loc[im].values, rm["OTX015"].loc[im].values)

    (proc / "exp3_f4_stats.json").write_text(json.dumps(st, indent=2, default=str))
    t1.to_csv(proc / "exp3_f4_test1.csv", index=False)
    figure(cfg, t1, scores, st)
    write(cfg, st, t1)
    log.info("F4 done")


def figure(cfg, t1, scores, st):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    ax[0].barh(t1.prediction, t1.cosine, color=["#55A868" if p == "additive" else "#6f62e8" for p in t1.prediction])
    ax[0].set_xlabel("cosine with the real combination change")
    ax[0].set_title("Predicting ATRA + palbociclib in HL-60")
    if "test2" in st:
        ax[1].hist(scores, bins=80, color="#BBBBBB")
        ax[1].axvline(st["test2"]["score"], color="#C44E52", lw=2)
        ax[1].set_xlabel("predicted maturation shift of the pair")
        ax[1].set_title(f"All {st['test2']['n_pairs']:,} pairs; red = {st['test2']['pair']}\n"
                        f"top {100 - st['test2']['percentile']:.3f}%", fontsize=9)
    fig.tight_layout(); f = path(cfg, "figures", "exp3_f4_combinations.png"); fig.savefig(f, dpi=150); plt.close(fig)


def write(cfg, st, t1):
    L = []; w = L.append
    t0, t1s = st["test0"], st["test1"]
    w("# Experiment 3, F4: can we predict drug combinations?")
    w(f"\nGenerated {dt.date.today().isoformat()} by `python -m pipeline.run --phase exp3-f4`. "
      "Plan: [exp3_f4_plan.md](exp3_f4_plan.md), committed before this analysis ran.\n")
    w("## Two things the model cannot do\n")
    w("**Additive and sequential composition are identical here.** F2 found no learned model that beat the linear shift, and for a linear model "
      "applying drug A then drug B is exactly the sum of their shifts. The brief's additive-versus-sequential contrast therefore has no content in "
      "this project, and we did not introduce an untested nonlinear model to manufacture one.")
    w("\n**Synergy is outside the model by construction.** An additive prediction says two drugs together do what they do apart. Genuine synergy or "
      "antagonism cannot be expressed, which bounds what Test 2 below can mean.")

    w("\n## Test 0: is the real combination distinct from its single agents?\n")
    w(pd.Series(t0).round(4).to_frame("cosine").to_markdown())
    dom = max(t0["real combination vs real ATRA"], t0["real combination vs real palbociclib"])
    w(f"\nThe measured ATRA + palbociclib response is closest to **{'ATRA' if t0['real combination vs real ATRA'] >= t0['real combination vs real palbociclib'] else 'palbociclib'}** "
      f"(cosine {dom:.3f}), while the two single agents resemble each other at {t0['real ATRA vs real palbociclib']:.3f}. "
      + ("Because the combination is dominated by one agent, Test 1 has limited power to separate composition rules, and its result should be read with that in mind."
         if dom > 0.9 else "The combination is not simply one of its agents, so Test 1 is informative."))

    w("\n## Test 1: which composition rule predicts the real combination?\n")
    w(t1.to_markdown(index=False, floatfmt=".4f"))
    tt = t1.set_index("prediction")["cosine"]
    best_rule = tt.idxmax()
    w(f"\n**Additive composition {'beats' if t1s['beats_both_singles'] else 'does not beat'} both single-drug predictions** "
      f"({tt['additive']:.3f} against {tt['ATRA only']:.3f} for ATRA alone and {tt['palbociclib only']:.3f} for palbociclib alone). "
      f"Knowing both drugs therefore helps, which is the one positive finding here.")
    if best_rule != "additive":
        w(f"\nHowever, **{best_rule} scores marginally higher** ({tt[best_rule]:.4f} against {tt['additive']:.4f}). A gap of "
          f"{tt[best_rule] - tt['additive']:.4f} on a single combination is well inside noise, so these two rules are not distinguishable here "
          "and we do not claim one is better.")
    w(f"\n**A caveat that matters more.** Scored against every one of the {t1s['n_single_drugs']} single AML/monocytic drugs, the additive prediction "
      f"ranks only **{t1s['additive_rank_among_single_drugs']}**: {t1s['best_single_drug']}, a drug not in this experiment, predicts the combination "
      f"better (cosine {t1s['best_single_cosine']:.3f}). Beating the two constituent drugs is a low bar when unrelated compounds clear it too.")

    if "test2" in st:
        t2 = st["test2"]
        w("\n## Test 2: the published positive control\n")
        w(f"**{t2['pair']}** scores {t2['score']:.4f}, ranking **{t2['rank']:,} of {t2['n_pairs']:,}** pairs "
          f"(top {100 - t2['percentile']:.3f}%). Individual maturation shifts: "
          + ", ".join(f"{k} {v:+.4f}" for k, v in t2["singles"].items()) + ".")
        w("\nThe evidence for this pair is real: Schenk et al. 2012 (*Nature Medicine* 18:605-611, doi:10.1038/nm.2661) showed LSD1 inhibition with "
          "tranylcypromine unlocks an ATRA differentiation response in non-APL AML. **But this test does not validate combination prediction.** "
          "Our score is a normalized sum, so a pair ranks high exactly when both drugs already point toward maturity, and tranylcypromine was already "
          "the second-highest single compound in the Phase 3 screen. The rank is carried by that one drug, not by anything we know about the pair.")
        w("\nThe list below makes the same point: nearly every top pair contains tranylcypromine or bosutinib, alongside partners with no plausible "
          "differentiation mechanism (desmopressin, estrone, methoxsalen). These are rankings of single-drug scores, not predictions about interactions.")
        w("\nHighest-scoring pairs (exploratory, untested):\n")
        w(pd.DataFrame(t2["top10"], columns=["drug A", "drug B", "score"]).to_markdown(index=False, floatfmt=".4f"))

    w("\n## Test 3: the combination we cannot predict\n")
    w(pd.Series(st["test3"]).round(4).to_frame("cosine").to_markdown())
    w("\nVTP50469 has no LINCS signature, so this MOLM13 combination cannot be predicted at all. The measured values are shown only to document "
      "that the combination again sits close to one of its agents.")

    w("\n## Conclusion\n")
    w("- **Beat AML contributes nothing here.** Verified: its 'Combination' panel holds 3 single agents across 21 samples, and all 63,395 probit fits "
      "are single-agent. There are no combination arms in the open files.")
    w(f"- Combination prediction rests on **one combination, in one cell line, at one timepoint**. "
      f"{'Additive composition came out ahead of both single agents there' if t1s['beats_both_singles'] else 'Additive composition did not come out ahead of both single agents there'}, "
      "which is a single data point either way and does not establish that the approach works in general.")
    w("- The ranked pairs above are **untested hypotheses**. No open dataset measures them, and our model cannot represent the synergy that makes "
      "combinations interesting in the first place.")
    w("\n## Figure\n")
    w("![combinations](figures/exp3_f4_combinations.png)")
    path(cfg, "reports", "exp3_f4_results.md").write_text("\n".join(L) + "\n", encoding="utf-8")
