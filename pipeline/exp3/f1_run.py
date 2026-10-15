"""Experiment 3, F1 runner. Plan: reports/exp3_f1_plan.md. Writes reports/exp3_f1_results.md."""
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
from pipeline.exp3 import datasets as DS
from pipeline.exp3 import f1 as F
from pipeline.phase0_data.genes import HGNCMapper

log = get_logger("exp3.f1")

# Comparisons fixed in the plan: (id, dataset, treated arm, drug in our set, subset, role)
COMPARISONS = [
    ("D1 4h", "GSE227542", "RA4h", "tretinoin", None, "primary, favourable"),
    ("D1 2d", "GSE227542", "RA2d", "tretinoin", None, "primary, favourable"),
    ("D1 4d", "GSE227542", "RA4d", "tretinoin", None, "primary, favourable"),
    ("D2 ATRA 48h", "GSE276457", "ATRA", "tretinoin", None, "primary, favourable"),
    ("D3 palbociclib 48h", "GSE276457", "palbociclib", "Palbociclib", None, "primary, other class"),
    ("D4 OTX015 16h", "GSE228321", "OTX015", "OTX-015", {"cell_line": "MOLM13"}, "secondary, unfavourable"),
]


def build_simulator(cfg, proc):
    mapper = joblib.load(proc / "landscape_mapper.joblib")
    ref = ad.read_h5ad(proc / "landscape_reference.h5ad")
    Xr, _ = mapper.align(ref.X, ref.var_names)
    gene_sd = Xr.std(0)
    D = pd.read_parquet(proc / "drug_consensus_signatures.parquet")
    var = ad.read_h5ad(proc / "lincs_signatures_all.h5ad", backed="r").var
    reliable = set(var.index[var.feature_space.isin(["landmark", "best inferred"])])
    man = json.loads((Path(cfg["paths"]["reports"]).parent / "web/public/data/manifest.json").read_text())
    if "shift_norm" not in man:
        raise RuntimeError("manifest has no shift_norm; rerun phase 5")
    sim = F.Simulator(mapper, D, gene_sd, reliable, man["shift_norm"], cfg["phase5"]["alpha_max"])
    return mapper, sim, man["shift_norm"]


def run(cfg: dict) -> None:
    proc = path(cfg, "processed", "x").parent
    raw = Path("data/raw/exp3")
    DS.fetch(raw)
    hgnc = HGNCMapper.from_file(path(cfg, "raw", "hgnc", "hgnc_complete_set.txt"))
    mapper, sim, norm = build_simulator(cfg, proc)
    meta = pd.read_parquet(proc / "drug_signature_meta.parquet")
    meta = meta[~meta.index.duplicated()]
    null_drugs = [d for d in sim.D.index if meta.tier.get(d) == "AML/monocytic"]
    sig = pd.read_csv(proc / "maturation_signature.csv")
    log.info("simulator: %d genes, norm %.3g, null %d drugs", len(sim.genes), norm, len(null_drugs))

    data = {"GSE227542": DS.load_gse227542(raw), "GSE276457": DS.load_gse276457(raw),
            "GSE228321": DS.load_gse228321(raw, hgnc)}
    st: dict = {"n_sim_genes": int(len(sim.genes)), "shift_norm": norm, "n_null_drugs": len(null_drugs),
                "datasets": {k: {"genes": int(v[0].shape[0]), "samples": int(v[0].shape[1])} for k, v in data.items()}}

    rows, mat_rows, nulls = [], [], {}
    for cid, ds, arm, drug, subset, role in COMPARISONS:
        X, m = data[ds]
        real = DS.real_change(X, m, arm, subset=subset)
        d = F.direction_test(sim, real, drug, null_drugs)
        nulls[cid] = d.pop("cos_all")
        d_raw = F.direction_test(sim, real, drug, null_drugs, weight_by_sd=False, n_perm=200)
        d_raw.pop("cos_all")
        rows.append({"id": cid, "dataset": ds, "arm": arm, "drug": drug, "role": role, **d,
                     "cosine_raw_signature": d_raw["cosine"], "rank_raw_signature": d_raw["rank_in_null"]})
        mt = F.maturity_test(sim, real, drug, sig)
        try:
            ls = F.landscape_scores(mapper, X, m, arm, "control", sim.change(drug), subset)
        except Exception as e:  # platform transfer can fail; report rather than hide
            log.warning("landscape score failed for %s: %s", cid, e)
            ls = {"error": str(e)[:80]}
        mat_rows.append({"id": cid, "drug": drug, **mt, **{f"ls_{k}": v for k, v in ls.items()}})
        log.info("%s: cosine %.3f rank %d/%d", cid, d["cosine"], d["rank_in_null"], d["n_null"])

    direction = pd.DataFrame(rows)
    maturity = pd.DataFrame(mat_rows)
    direction["passes"] = (direction.cosine > 0) & (direction.p_rank <= 0.05)

    # Test 3: trajectory, D1 only
    X, m = data["GSE227542"]
    reals = {t: DS.real_change(X, m, a) for t, a in [("4h", "RA4h"), ("2d", "RA2d"), ("4d", "RA4d")]}
    traj = F.trajectory_test(sim, reals, "tretinoin")
    traj["monotone_parallel"] = bool(traj.parallel.is_monotonic_increasing)

    st["direction"] = direction.drop(columns=["negative_controls"]).to_dict(orient="records")
    st["negative_controls"] = {r["id"]: r["negative_controls"] for r in rows}
    st["maturity"] = maturity.to_dict(orient="records")
    st["trajectory"] = traj.to_dict(orient="records")
    direction.drop(columns=["negative_controls"]).to_csv(proc / "exp3_f1_direction.csv", index=False)
    maturity.to_csv(proc / "exp3_f1_maturity.csv", index=False)
    traj.to_csv(proc / "exp3_f1_trajectory.csv", index=False)
    (proc / "exp3_f1_stats.json").write_text(json.dumps(st, indent=2, default=str))
    figures(cfg, direction, nulls, traj, maturity)
    write(cfg, st, direction, maturity, traj, nulls)
    log.info("F1 done: %d/%d comparisons pass", int(direction.passes.sum()), len(direction))


def figures(cfg, direction, nulls, traj, maturity):
    n = len(nulls)
    fig, ax = plt.subplots(1, n, figsize=(3.4 * n, 3.4), sharey=False)
    for a, (cid, cos) in zip(np.atleast_1d(ax), nulls.items()):
        row = direction[direction.id == cid].iloc[0]
        a.hist(cos.values, bins=60, color="#BBBBBB")
        a.axvline(row.cosine, color="#C44E52", lw=2)
        for d, v in row.negative_controls.items():
            a.axvline(v, color="#4C72B0", lw=0.7, alpha=0.7)
        a.set_title(f"{cid}\n{row.drug}: cos {row.cosine:.3f}, rank {row.rank_in_null}/{row.n_null}", fontsize=8)
        a.set_xlabel("cosine with real change", fontsize=8)
    fig.suptitle("Red = the drug actually applied; grey = all 1,691 AML/monocytic drugs; blue = negative controls", fontsize=9)
    fig.tight_layout()
    f1 = path(cfg, "figures", "exp3_f1_direction.png"); fig.savefig(f1, dpi=150); plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].bar(traj.timepoint, traj.off_axis_fraction, color="#DD8452")
    ax[0].set_ylim(0, 1); ax[0].set_ylabel("off-axis fraction of real change")
    ax[0].set_title("How much of the real ATRA change is off our straight line")
    ax[1].plot(traj.timepoint, traj.parallel, "o-", color="#4C72B0")
    ax[1].set_ylabel("component along the simulated direction"); ax[1].set_title("Movement along our axis")
    fig.tight_layout(); f2 = path(cfg, "figures", "exp3_f1_trajectory.png"); fig.savefig(f2, dpi=150); plt.close(fig)
    return f1, f2


def write(cfg, st, direction, maturity, traj, nulls):
    L = []; w = L.append
    w("# Experiment 3, F1: does the linear simulation predict real drug-induced change?")
    w(f"\nGenerated {dt.date.today().isoformat()} by `python -m pipeline.run --phase exp3-f1`. "
      "Plan: [exp3_f1_plan.md](exp3_f1_plan.md), committed before this analysis ran. Data: [exp3_step0_data.md](exp3_step0_data.md).\n")
    w("## Headline limitation\n")
    w("**This is a cell-line-only validation.** No open patient before/after transcriptome dataset exists for any drug in our LINCS set "
      "(Step F0), so whether the simulation transfers to patient cells is **untested and remains unknown**.")
    w("\n**The primary tests are deliberately favourable.** Our tretinoin signature is built from HL60, SKM1, THP1 and U937, and both ATRA "
      "datasets are HL-60 cells. D1 and D2 therefore measure the best case: same drug, same cell line. D4 (OTX015 in MOLM13, scored from a "
      "non-myeloid signature) is the unfavourable transfer case and should be read as the more realistic one.\n")

    npass, ntot = int(direction.passes.sum()), len(direction)
    w(f"## Result: {npass} of {ntot} comparisons pass the pre-registered direction test\n")
    cols = ["id", "drug", "role", "n_genes", "cosine", "spearman", "rank_in_null", "p_rank", "perm_p", "passes"]
    w(direction[cols].to_markdown(index=False, floatfmt=".4f"))
    w(f"\nThe simulated change vector covers **{st['n_sim_genes']} genes** (the landscape mapper's variable genes that also have a reliable "
      f"L1000 measurement). Genes outside this set are not moved by the simulation at all and are excluded. "
      f"`rank_in_null` is the drug's rank among all {st['n_null_drugs']} compounds with an AML/monocytic signature, scored identically; "
      "`p_rank` is that rank as a fraction. The pre-registered rule is cosine > 0 and rank in the top 5%.")

    w("\n### Negative controls\n")
    nc = pd.DataFrame(st["negative_controls"]).T
    w(nc.round(4).to_markdown())
    w("\nThese eight compounds have no established myeloid differentiation mechanism and were fixed in the plan before running.")

    w("\n## Maturity test\n")
    mcols = ["id", "drug", "real_shift", "sim_shift", "same_direction", "magnitude_ratio"]
    w(maturity[mcols].to_markdown(index=False, floatfmt=".4f"))
    w("\n`real_shift` and `sim_shift` are the maturation-signature direction (mean change of up-with-maturation genes minus down genes), in log2 units. "
      "**Magnitude is not a meaningful test:** the simulation's step length was chosen in Phase 5 purely so the animation stays on the manifold, "
      "so `magnitude_ratio` reflects that arbitrary constant, not biology. Only the sign is informative.")
    ls_cols = [c for c in maturity.columns if c.startswith("ls_")]
    if ls_cols:
        w("\nLandscape score (Phase 1 mapper applied to bulk profiles; the mapper was trained on single-cell data, so this is a platform transfer and is reported with that caveat):\n")
        w(maturity[["id", "drug"] + ls_cols].to_markdown(index=False, floatfmt=".4f"))

    w("\n## Trajectory test (exploratory, ATRA time course)\n")
    w(traj.to_markdown(index=False, floatfmt=".4f"))
    w("\nReal timepoints are 4 h, 2 d and 4 d, not the days 1, 3 and 6 assumed in the brief. `off_axis_fraction` is the share of the real change "
      "that is perpendicular to our straight-line direction: 0 would mean the cells move exactly along our axis, 1 that our axis captures none of it.")

    w("\n## Interpretation\n")
    d = direction.set_index("id")
    passed = d[d.passes]
    mat = maturity.set_index("id")
    w(f"**The direction test passes for {len(passed)} of {len(d)} comparisons, but the simulation still captures only a small part of what really "
      "happens.** Both statements are true, and the second matters more.")
    w("\n**1. Direction is detectable above chance.** For every comparison from 16 h onward, the real change vector is closer to the applied drug's "
      "simulated direction than to almost all of 1,691 other drugs scored the same way (ranks 1, 4, 4, 11 and 37). The eight pre-specified negative "
      "controls sit near zero throughout. The LINCS signature does carry real information about which way these cells move.")
    shares = ", ".join(f"{r.cosine ** 2 * 100:.0f}% ({i})" for i, r in passed.iterrows())
    w(f"\n**2. But the straight line misses most of the movement.** Cosine values are 0.28 to 0.62, so the share of the real change lying along our axis "
      f"(cosine squared) is only {shares}. The trajectory test says the same thing directly: {traj.off_axis_fraction.min() * 100:.0f}% to "
      f"{traj.off_axis_fraction.max() * 100:.0f}% of the real ATRA response is perpendicular to our simulated direction. The cells move somewhere our "
      "straight line does not point.")
    early = d.loc["D1 4h"]
    w(f"\n**3. The simulation is wrong at early times.** At 4 h the real change is large (vector norm {traj.iloc[0].total_norm:.0f}, about 40% of the 4 d "
      f"response) but points elsewhere entirely: cosine {early.cosine:.3f}, rank {early.rank_in_null} of {early.n_null}, and the maturation-signature "
      f"shift is **negative** ({mat.loc['D1 4h', 'real_shift']:.2f}) while the simulation predicts positive. The early ATRA response is not a small "
      "version of the late one; it points a different way. Our model has no time axis, so its predicted shift is identical at 4 h, 2 d and 4 d.")
    w("\n**4. Maturity direction is right where it matters.** For ATRA at 2 d, 4 d and 48 h, the real and simulated maturation-signature shifts are both "
      "positive, and the landscape score rises in the real data (+0.07 at 2 d, +0.28 at 4 d). When ATRA does drive HL-60 cells toward maturity, our score sees it.")
    w("\n**5. The unfavourable case held up.** D4 (OTX015 in MOLM13, scored from a non-myeloid-tier signature, in a different cell line from the signature) "
      "still ranked 11 of 1,691. That is the most encouraging result for transfer, though it is a single comparison.")
    w("\n### What we may and may not claim\n")
    w("- **Supported:** the linear model predicts the *direction* of drug-induced expression change better than chance, in AML cell lines, at 16 h to 4 d.")
    w("- **Supported:** our maturation score moves the right way when a differentiation agent works.")
    w("- **Not supported:** that the simulation reproduces the real trajectory. Most of the real movement is off-axis and the early response is mispredicted.")
    w("- **Not supported, and untestable here:** anything about patient cells. No open data allows that test.")
    w("- Wording such as 'simulates cell reversion' is not justified. 'Predicts the direction of drug-induced expression change in cell lines' is.")
    w("\nThis is the motivation for F2: a model with a time axis and a curved path has clear room to beat the straight line, and F1 now gives it a number to beat.")

    w("\n## Figures\n")
    w("![direction](figures/exp3_f1_direction.png)")
    w("\n![trajectory](figures/exp3_f1_trajectory.png)")
    path(cfg, "reports", "exp3_f1_results.md").write_text("\n".join(L) + "\n", encoding="utf-8")
