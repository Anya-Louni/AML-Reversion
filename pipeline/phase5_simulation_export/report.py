"""Phase 5 report."""
from __future__ import annotations

import datetime as dt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from pipeline.common import path


def write_report(cfg, sims: pd.DataFrame, size: int, thr: float, norm: dict) -> None:
    p5 = cfg["phase5"]
    piv = sims.pivot_table(index="patient", columns="drug", values="pct_matured_of_stuck")
    fig, ax = plt.subplots(figsize=(14, 5))
    im = ax.imshow(piv.values, cmap="viridis", vmin=0, vmax=100, aspect="auto")
    ax.set_yticks(range(len(piv)), piv.index, fontsize=8); ax.set_xticks(range(piv.shape[1]), piv.columns, rotation=90, fontsize=7)
    plt.colorbar(im, ax=ax, fraction=0.02, label="% of stuck cells that cross the threshold")
    ax.set_title("Simulated maturation (linear expression-shift model); blank = drug not in that patient's top list")
    fig.tight_layout()
    f = path(cfg, "figures", "phase5_simulation_outcomes.png"); fig.savefig(f, dpi=150); plt.close(fig)

    L = []; w = L.append
    w("# Phase 5: simulation export")
    w(f"\nGenerated {dt.date.today().isoformat()} by `python -m pipeline.run --phase 5`.\n")
    w("> **What the motion means.** Cells move because their log expression is shifted linearly along a drug's LINCS consensus signature: "
      "x' = x + alpha * v, where v is the drug's z-score per gene (landmark and best-inferred genes only) times the gene's SD in the normal reference, rescaled to a common "
      f"length for every drug (see caveats), and alpha rises linearly from 0 to {p5['alpha_max']} over {p5['n_frames']} frames. Each frame is placed with the fixed Phase 1 map. This is a "
      "visualization of a hypothesis. It is not a model of cell fate, dose, time, toxicity or resistance mechanisms.\n")
    w("## Setup\n")
    w(f"- Patients: the {sims.patient.nunique()} van Galen patients with malignant cells at diagnosis; up to {p5['max_cells_per_patient']} cells each "
      "(random subsample with the fixed seed; all cells if fewer).")
    w(f"- Drugs per patient: the top {p5['top_n_drugs']} named compounds in that patient's Phase 3 ranking whose signature comes from AML/monocytic lines, "
      f"plus {', '.join(p5['always_include'])} as a reference positive control.")
    w(f"- Maturity threshold: {thr:.3f} (Phase 2). A cell that starts below it is **matured** if it ends above it, **resistant** if it never crosses it, "
      "and **transient** if it crosses and falls back.")
    w(f"- Export: `web/public/data/`, {size/1e6:.1f} MB in total. Trajectories are uint16-quantized (x, y, score) per frame per cell.")
    w("\n## Outcomes\n")
    t = sims[["patient", "drug", "n", "already_mature", "stuck_start", "matured", "resistant", "transient", "pct_matured_of_stuck",
              "mean_score_start", "mean_score_end"]]
    w(t.to_markdown(index=False, floatfmt=".2f"))
    w(f"\n![simulation outcomes](figures/{f.name})")
    w("\n## Caveats specific to the motion\n")
    by = sims.groupby("drug").agg(pct_matured=("pct_matured_of_stuck", "mean"), offmanifold_end=("offmanifold_end", "mean"),
                                  raw_shift_norm=("raw_shift_norm", "mean")).sort_values("pct_matured", ascending=False)
    rho = by[["pct_matured", "raw_shift_norm"]].corr(method="spearman").iloc[0, 1]
    scan = ", ".join(f"{k}: {v*100:.1f}%" for k, v in norm["scan"].items())
    w(f"- **Common shift length.** Raw LINCS signatures differ hugely in magnitude; with raw signatures, 68% of final-frame cells ended up off the normal "
      f"manifold and the animation mostly showed magnitude. Each drug's shift is therefore rescaled to the same length, so the motion reflects direction (what the "
      f"ranking scores). The length was set by a fixed rule: the largest candidate keeping at most {cfg['phase5']['max_offmanifold_end']*100:.0f}% of final-frame "
      f"cells off-manifold on average. Scan (length: mean off-manifold share) = {scan}; chosen length **{norm['target_norm']}**.")
    w(f"- After rescaling, % matured and raw signature size correlate at Spearman {rho:.2f} across drugs. Large raw signatures no longer win by size; the "
      "negative sign means drugs with the largest raw LINCS responses change many genes that are not on the maturation axis, so once length is equalized they move cells least.")
    w(f"- **Off-manifold endpoints.** At the last frame, on average {sims.offmanifold_end.mean()*100:.0f}% of cells are further from any normal reference cell than the Phase 2 "
      f"off-manifold cut (vs {sims.offmanifold_start.mean()*100:.0f}% at the start). Those cells are drawn at their nearest normal neighbours, so their final position overstates how "
      "'normal' the shifted expression is. The animation shows a direction of change, not a claim that cells become normal monocytes.\n")
    w(by.round(3).to_markdown(floatfmt=".3f"))
    w("\n## Choices and justifications\n")
    rows = [("Shift in units of reference gene SD.", "LINCS z-scores are relative to plate controls; scaling by each gene's spread in the reference keeps shifts comparable across genes with very different expression levels."),
            ("Same alpha schedule for every drug.", "Differences between drugs then come only from their signatures, not from a per-drug tuning."),
            ("Only landmark and best-inferred genes are moved.", "The other inferred genes are least reliable; moving them would add noise."),
            ("Clip expression at zero.", "Log-normalized expression cannot be negative."),
            ("Always include tretinoin.", "Lets every patient be compared with the clinically used differentiation agent."),
            ("Subsample to 300 cells.", "Keeps the static site small and fast while showing each patient's distribution.")]
    w(pd.DataFrame(rows, columns=["Choice", "Justification"]).to_markdown(index=False))
    path(cfg, "reports", "phase5_simulation_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
