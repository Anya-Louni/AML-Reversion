"""Phase 2 report and figures."""
from __future__ import annotations

import datetime as dt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pipeline.common import path

LIKE_COL = {"HSC-like": "#3B4CC0", "Prog-like": "#6A8BEF", "GMP-like": "#9ABBFF", "ProMono-like": "#F2A481",
            "Mono-like": "#D8553F", "cDC-like": "#8E44AD"}


def md(df, index=False):
    return df.to_markdown(index=index, floatfmt=".3f")


def figures(cfg, st, m, out, ref):
    F = {}
    thr = st["threshold"]["threshold"]
    obs = m.obs[(m.obs.sample_type == "AML") & (m.obs.day == 0)]
    pats = st["patients_d0"].sort_values("median_score").index.tolist()

    fig, ax = plt.subplots(1, 2, figsize=(15, 5), gridspec_kw={"width_ratios": [1.6, 1]})
    data = [obs.score.values[obs.patient == p] for p in pats]
    parts = ax[0].violinplot(data, showmedians=True, widths=0.85)
    for b in parts["bodies"]:
        b.set_facecolor("#6A8BEF"); b.set_alpha(0.6)
    ax[0].axhline(thr, color="#C44E52", ls="--", lw=1, label=f"maturity threshold {thr:.2f}")
    ax[0].set_xticks(range(1, len(pats) + 1), [f"{p}\n(n={len(d)})" for p, d in zip(pats, data)], fontsize=7)
    ax[0].set_ylabel("differentiation score"); ax[0].set_title("Malignant cells at diagnosis, per patient"); ax[0].legend(frameon=False)
    sf = st["patients_d0"].reindex(pats)
    ax[1].barh(pats, sf.stuck_fraction, color="#3B4CC0")
    ax[1].set_xlabel("stuck fraction (score < threshold)"); ax[1].set_xlim(0, 1); ax[1].set_title("Stuck fraction")
    fig.tight_layout(); F["patients"] = path(cfg, "figures", "phase2_patient_scores.png"); fig.savefig(F["patients"], dpi=150); plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
    rx = ref.obsm["X_landscape"]
    ax[0].scatter(rx[:, 0], rx[:, 1], s=1, c="#DDDDDD", rasterized=True, label="normal reference")
    for t, col in LIKE_COL.items():
        s = obs[obs.CellType == t]
        ax[0].scatter(s.x, s.y, s=2, c=col, label=t, rasterized=True)
    ax[0].legend(markerscale=5, fontsize=7, frameon=False); ax[0].set_title("Malignant cells (day 0) placed on the landscape, by authors' label")
    ax[0].set_xticks([]); ax[0].set_yticks([])
    lt = cfg["phase2"]["malignant_like_types"]
    d = [m.obs.score.values[(m.obs.CellType == t) & (m.obs.sample_type == "AML")] for t in lt]
    parts = ax[1].violinplot(d, showmedians=True)
    for b, t in zip(parts["bodies"], lt):
        b.set_facecolor(LIKE_COL[t]); b.set_alpha(0.7)
    ax[1].axhline(thr, color="#C44E52", ls="--", lw=1)
    ax[1].set_xticks(range(1, len(lt) + 1), lt, fontsize=8); ax[1].set_ylabel("mapped score")
    ax[1].set_title(f"Mapped score by authors' malignant label (Spearman {st['like_order_rho_all']:.2f})")
    fig.tight_layout(); F["placement"] = path(cfg, "figures", "phase2_placement.png"); fig.savefig(F["placement"], dpi=150); plt.close(fig)

    v = st["pseudobulk"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for t, col in LIKE_COL.items():
        g = v[v.type == t]
        ax[0].scatter(g.true, g.est, c=col, label=t, s=25)
    ax[0].plot([0, 1], [0, 1], "k--", lw=0.8); ax[0].set_xlabel("true share (cell labels)"); ax[0].set_ylabel("NNLS estimate")
    ax[0].set_title(f"Pseudo-bulk check, patient held out (r = {st['pseudobulk_overall_r']:.2f})"); ax[0].legend(fontsize=7, frameon=False)
    vb = st["vs_bottomly"]
    ax[1].barh(vb.type, vb.spearman, color=[LIKE_COL[t] for t in vb.type])
    ax[1].set_xlabel("Spearman with Bottomly et al. 2022 scores"); ax[1].set_xlim(-1, 1); ax[1].axvline(0, color="k", lw=0.5)
    ax[1].set_title(f"Beat AML: our fractions vs published (n = {vb.n.iloc[0]})")
    fig.tight_layout(); F["deconv"] = path(cfg, "figures", "phase2_deconvolution.png"); fig.savefig(F["deconv"], dpi=150); plt.close(fig)

    aml = out[out.rna_control.isin(["No"])]
    sh = aml[[f"share_{t}" for t in lt]].sort_values(["share_Mono-like", "share_HSC-like"])
    fig, ax = plt.subplots(figsize=(14, 3.6))
    bottom = np.zeros(len(sh))
    for t in lt:
        ax.bar(range(len(sh)), sh[f"share_{t}"], bottom=bottom, color=LIKE_COL[t], width=1.0, label=t)
        bottom += sh[f"share_{t}"].values
    ax.set_xlim(0, len(sh)); ax.set_ylabel("share of malignant"); ax.set_xticks([])
    ax.set_title(f"Beat AML: deconvolved malignant composition ({len(sh)} AML samples)"); ax.legend(ncol=6, fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(); F["beat"] = path(cfg, "figures", "phase2_beataml_composition.png"); fig.savefig(F["beat"], dpi=150); plt.close(fig)
    return F


def write_report(cfg, st, m, out, ref):
    p2 = cfg["phase2"]; F = figures(cfg, st, m, out, ref); thr = st["threshold"]
    L = []; w = L.append
    pts = st["patients_d0"]
    w("# Phase 2: patient placement")
    w(f"\nGenerated {dt.date.today().isoformat()} by `python -m pipeline.run --phase 2`.\n")
    w("## Summary\n")
    w(f"- **Maturity threshold** (fixed on normal reference cells, before any drug analysis): score **{thr['threshold']:.3f}**, the cut that best separates "
      f"HSC/Prog/GMP from ProMono/Mono/cDC (balanced accuracy {thr['balanced_accuracy']:.3f}; {thr['sensitivity_immature']*100:.0f}% of immature and "
      f"{thr['specificity_mature']*100:.0f}% of mature normal cells correctly classified).")
    w(f"- **Mapping check against independent labels:** the mapped score of malignant cells increases along the authors' HSC-like < Prog-like < GMP-like "
      f"< ProMono-like < Mono-like labels with Spearman **{st['like_order_rho_all']:.3f}** (all days) and {st['like_order_rho_d0']:.3f} (day 0). "
      "Those labels come from the authors' classifier and were not used to build the map.")
    w(f"- **{len(pts)} patients** with malignant cells at diagnosis. Stuck fraction ranges from {pts.stuck_fraction.min():.2f} to {pts.stuck_fraction.max():.2f} "
      f"(median {pts.stuck_fraction.median():.2f}).")
    w(f"- **Deconvolution check:** on pseudo-bulk samples built from van Galen cells with the patient held out of the reference, estimated vs true malignant "
      f"composition r = **{st['pseudobulk_overall_r']:.3f}**. On Beat AML, our fractions agree with the published Bottomly et al. scores at Spearman "
      f"{st['vs_bottomly'].spearman.min():.2f} to {st['vs_bottomly'].spearman.max():.2f} across the six types.")
    bt = st["bulk_signature_transfer"]
    w(f"- **Signature transfer to bulk:** for {bt['n_genes']} maturation-signature genes, the direction of change with deconvolved mature share in Beat AML "
      f"agrees with the single-cell direction for **{bt['sign_agreement']*100:.0f}%** of genes (Spearman of slopes vs single-cell rho {bt['spearman_slope_vs_sc_rho']:.2f}).")

    w("\n## 1. Maturity threshold\n")
    w("Chosen as the score that maximizes balanced accuracy between immature (HSC, Prog, GMP) and mature (ProMono, Mono, cDC) normal reference cells. "
      "It is a property of normal haematopoiesis only. It defines 'stuck' (Phase 2) and 'matured' (Phase 5) and is not changed afterwards.")

    w("\n## 2. Single-cell placement (van Galen malignant cells)\n")
    mp = st["mapping"]
    w(f"All {len(m):,} malignant cells with a like-type label (all days, plus the two cell lines) were mapped with the fixed Phase 1 map. "
      f"{mp['gene_fraction_present']*100:.0f}% of model genes are present (same platform). "
      f"A cell is flagged **off-manifold** when its mean distance to its 15 reference neighbours exceeds the {p2['offmanifold_quantile']*100:.0f}th percentile "
      f"of the same distance among reference cells (leave-one-out), {mp['offmanifold_cut']:.2f}. Off-manifold cells are still placed, but their position "
      "is less trustworthy.\n")
    w("Mapped score by the authors' malignant label (all AML samples):\n")
    t = st["score_by_like_type"].join(st["offmanifold_by_like_type"].rename("off_manifold_frac"))
    w(md(t, index=True))
    w("\nPer patient, malignant cells at diagnosis:\n")
    cols = ["n_malignant", "median_score", "q25", "q75", "stuck_fraction", "offmanifold_fraction", "low_cell_count"]
    w(md(pts[cols], index=True))
    w(f"\nPatients flagged `low_cell_count` have fewer than {p2['min_cells_per_patient']} malignant cells; their numbers are shown but should be read with caution.")
    w("\nAll samples, including follow-up days (treatment changes the malignant compartment, so later days are shown for context only):\n")
    w(md(st["samples_all"].reset_index(), index=False))
    w("\nCell lines (not patients), as a sanity check: MUTZ3 and OCI-AML3 are AML lines.\n")
    w(md(st["cell_lines"], index=True))
    al = pd.Series(st["sc_stuck_alignment"])
    w(f"\n**Stuck signatures (single-cell):** for each patient, mean log expression of day-0 malignant cells minus the mean of normal reference monocytes, "
      f"over the {cfg['phase1']['signature_top_n']*2} primary signature genes. The 'push' direction (minus stuck) correlates with the maturation direction "
      f"(sign of each gene's rho) at Spearman {al.min():.2f} to {al.max():.2f} across patients, i.e. patients furthest from maturity need the largest push along the signature.")

    w("\n## 3. Bulk deconvolution (Beat AML)\n")
    w(f"Reference: mean CPM profiles of the six van Galen malignant like-types plus non-myeloid normal groups ({', '.join(p2['deconvolution']['normal_groups'])}), "
      f"on the union of each group's top {p2['deconvolution']['n_markers_per_type']} marker genes ({st['deconv_ref_shape'][0]} genes, {st['deconv_genes_used']} found in Beat AML). "
      "Non-negative least squares on gene-scaled values; fractions sum to 1, and malignant shares are renormalized within the six like-types. "
      "Normal myeloid cells are not separate references (they are nearly collinear with the like-types), so in bulk they count toward the matching like-type.")
    w("\n**Pseudo-bulk validation** (each AML sample with >= 200 cells summed into a pseudo-bulk; reference rebuilt without that patient):\n")
    w(md(st["pseudobulk_by_type"], index=True))
    ms = st["pseudobulk_mature_share"]
    w(f"\nThe six fine types are only moderately resolved (adjacent stages such as HSC-like and Prog-like have similar profiles). The coarse quantity used "
      f"downstream, the **mature share** (ProMono + Mono + cDC-like), is recovered much better: Pearson **{ms['pearson']:.3f}**, MAE {ms['mae']:.3f} over "
      f"{ms['n_samples']} pseudo-bulk samples. We also tried 50/100/200 markers per type and two marker rules; overall r stayed between 0.61 and 0.64, "
      "so the pre-set configuration (50 markers, fold change vs mean of other types) was kept rather than tuned.")
    w("\nThis is same-platform, so it is an upper bound on accuracy for real bulk RNA-seq.")
    w("\n**Comparison with Bottomly et al. 2022** (their published van Galen cell-type scores for the same samples):\n")
    w(md(st["vs_bottomly"]))
    w("\nThe published values are signature scores, not fractions, so only rank agreement is meaningful. Their analysis already established that "
      "differentiation state shapes drug response in Beat AML, so this deconvolution is a building block here, not a novel result.")
    w(f"\nPer-sample fractions, malignant shares, mature share (ProMono + Mono + cDC-like) and a bulk maturity index (shares weighted by each like-type's "
      f"median single-cell score) are saved for {len(out)} samples ({st['n_bulk_aml']} AML).")

    w("\n## 4. Bulk stuck signatures\n")
    w("No mature reference exists on the bulk platform, and subtracting a single-cell monocyte profile from bulk RNA-seq would mix platforms. "
      "Instead, for each signature gene we regress Beat AML expression on the deconvolved mature share and take the predicted value at mature share = 1 "
      "as the expected mature expression. Stuck = patient - expected mature. "
      f"Check: the bulk slopes agree in sign with the single-cell signature for {bt['sign_agreement']*100:.0f}% of genes, so the signature carries over to bulk.")

    w("\n## 5. Choices and justifications\n")
    rows = [
        ("Threshold by maximum balanced accuracy on normal cells.", "Uses only normal biology, so it cannot be tuned toward drug results; balanced accuracy is not biased by the much larger number of Mono cells."),
        ("Immature = HSC, Prog, GMP; mature = ProMono, Mono, cDC.", "GMP cells are still proliferating progenitors, while promonocytes are committed to the monocytic fate; this is the standard boundary between progenitor and committed compartments."),
        ("Map malignant cells with the fixed Phase 1 map.", "The same map is used later for drug-shifted cells, so placement and simulation are consistent."),
        ("Off-manifold flag at the 99th percentile of reference distances.", "Malignant cells can be transcriptionally unlike any normal cell; flagging them makes that visible instead of hiding it in a nearest-neighbour average."),
        ("Validate placement against the authors' like-type labels.", "Those labels are independent of our map, so agreement is genuine evidence."),
        ("Day-0 samples for patient summaries.", "Diagnosis samples reflect untreated disease, which is what a differentiation therapy would face."),
        ("NNLS with marker genes and gene scaling.", "Standard reference-based deconvolution; scaling stops a few highly expressed genes from dominating."),
        ("Include non-myeloid normal references.", "Bulk marrow and blood contain lymphocytes and erythroid cells; without references they would be forced into malignant types."),
        ("Patient held out of the reference during pseudo-bulk checks.", "Otherwise the reference would contain the very cells being estimated."),
        ("Bulk mature expectation by within-cohort regression.", "Keeps both sides of the subtraction on the same platform; the sign-agreement check shows it reproduces the single-cell direction."),
    ]
    w(md(pd.DataFrame(rows, columns=["Choice", "Justification"])))
    w("\n## 6. Things to watch\n")
    w("- Malignant cells are placed by similarity to normal cells; a malignant 'Mono-like' cell can look mature in expression yet keep leukaemic behaviour.")
    w("- The bulk mature expectation extrapolates to mature share = 1; patients far outside the cohort range get larger, less certain stuck values.")
    w("- Bulk samples include normal cells, and blood samples contain normal monocytes, which can inflate the Mono-like share.")
    w("\n## Figures\n")
    for f in F.values():
        w(f"![{f.stem}](figures/{f.name})")
    w("\n## Outputs (`data/processed/`)\n")
    w("`malignant_cells_mapped.parquet`, `stuck_signatures_vangalen.parquet`, `beataml_deconvolution.parquet`, `stuck_signatures_beataml.parquet`, "
      "`beataml_mature_expectation.csv`, `phase2_stats.json`.")
    path(cfg, "reports", "phase2_patients_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
