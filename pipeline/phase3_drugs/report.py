"""Phase 3 report and figures."""
from __future__ import annotations

import datetime as dt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pipeline.common import path


def md(df, index=False):
    return df.to_markdown(index=index, floatfmt=".3f")


def figures(cfg, st, G, out):
    F = {}
    bg = G[G.tier == "AML/monocytic"]
    ctl = st["controls"].dropna(subset=["percentile_in_myeloid_screen"])
    fig, ax = plt.subplots(1, 2, figsize=(14, 4.8))
    ax[0].hist(bg.cosine, bins=60, color="#BBBBBB")
    for _, r in ctl.iterrows():
        col = "#C44E52" if r.tier == "tier1" else "#4C72B0"
        ax[0].axvline(r.cosine, color=col, lw=1.2, alpha=0.9)
        ax[0].text(r.cosine, ax[0].get_ylim()[1] * (0.95 if r.tier == "tier1" else 0.75), r.drug, rotation=90, fontsize=6, color=col, va="top")
    ax[0].set_xlabel("generic maturation score (cosine with maturation signature)")
    ax[0].set_ylabel("compounds"); ax[0].set_title(f"{len(bg)} compounds in AML/monocytic lines (red = tier 1, blue = tier 2 controls)")
    ax[1].scatter(G.cosine, G.connectivity, s=4, c="#4C72B0", alpha=0.5, rasterized=True)
    ax[1].set_xlabel("cosine"); ax[1].set_ylabel("KS connectivity")
    ax[1].set_title(f"Two scoring methods (Spearman {st['generic']['cosine_vs_connectivity']:.2f})")
    fig.tight_layout(); F["controls"] = path(cfg, "figures", "phase3_positive_controls.png"); fig.savefig(F["controls"], dpi=150); plt.close(fig)

    cos = out["vangalen"][0]
    top = cos.mean().nlargest(30).index
    M = cos[top]
    fig, ax = plt.subplots(figsize=(14, 5))
    v = np.abs(M.values).max()
    im = ax.imshow(M.values, cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
    ax.set_yticks(range(len(M)), M.index, fontsize=8); ax.set_xticks(range(len(top)), top, rotation=90, fontsize=7)
    plt.colorbar(im, ax=ax, fraction=0.02, label="reversion score (cosine)")
    ax.set_title("Reversion scores, van Galen patients x 30 drugs with the highest mean score")
    fig.tight_layout(); F["heatmap"] = path(cfg, "figures", "phase3_ranking_heatmap.png"); fig.savefig(F["heatmap"], dpi=150); plt.close(fig)

    cosb = out["beataml"][0]
    topb = cosb.mean().sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(14, 4))
    q = cosb[topb.index]
    ax.boxplot([q[c].values for c in q.columns], showfliers=False)
    ax.set_xticks(range(1, len(q.columns) + 1), q.columns, rotation=90, fontsize=5)
    ax.axhline(0, color="k", lw=0.5); ax.set_ylabel("reversion score (cosine)")
    ax.set_title(f"Beat AML: reversion score distribution per drug across {len(cosb)} patients (Beat AML panel drugs)")
    fig.tight_layout(); F["beat"] = path(cfg, "figures", "phase3_beataml_scores.png"); fig.savefig(F["beat"], dpi=150); plt.close(fig)
    return F


def write_report(cfg, st, G, out):
    p3 = cfg["phase3"]; F = figures(cfg, st, G, out); g = st["generic"]
    L = []; w = L.append
    w("# Phase 3: drug ranking")
    w(f"\nGenerated {dt.date.today().isoformat()} by `python -m pipeline.run --phase 3`.\n")
    ct = st["control_tests"]
    w("## Summary\n")
    w(f"- **{st['n_drugs']:,} compounds** scored ({st['n_signatures']:,} LINCS signatures). Line tier used per compound: {st['tier_counts']}.")
    w(f"- **Positive controls, tier 1** (retinoids, vitamin D, PKC activators) vs all other compounds profiled in AML/monocytic lines: "
      f"AUROC **{ct['tier1_cosine']['auroc']:.3f}** (permutation p = {ct['tier1_cosine']['perm_p']:.4f}, {ct['tier1_cosine']['n_controls']} controls, "
      f"{ct['tier1_cosine']['n_background']} compounds) with the cosine score; {ct['tier1_connectivity']['auroc']:.3f} (p = {ct['tier1_connectivity']['perm_p']:.4f}) with connectivity.")
    w(f"- **Tier 2** (weaker evidence): AUROC {ct['tier2_cosine']['auroc']:.3f} (p = {ct['tier2_cosine']['perm_p']:.4f}).")
    w(f"- Robustness of the generic score: cosine vs connectivity Spearman {g['cosine_vs_connectivity']:.2f}; primary genes vs landmark-only {g['primary_vs_landmark']:.2f}; "
      f"primary vs full signature (with ribosomal genes) {g['primary_vs_full_signature']:.2f}. Correlation with signature strength (TAS): {g['cosine_vs_tas']:.2f}.")
    w(f"- **Patient specificity:** rankings are dominated by drug-level effects. The average Spearman between two patients' drug rankings is "
      f"{st['vangalen_between_patient_rank_corr']:.2f} (van Galen) and {st['beataml_between_patient_rank_corr']:.2f} (Beat AML).")

    w("\n## 1. Drug signatures\n")
    w(f"Besides the Phase 0 set, we fetched every compound signature in AML/monocytic lines ({', '.join(p3['wide_cell_lines'])}), up to "
      f"{p3['wide_max_sigs_per_cell']} per compound and line (same quality ranking as Phase 0). Per compound, signatures are averaged within each cell line, "
      "then across lines of the best available tier: AML/monocytic lines first, then CD34/K562, then the core non-myeloid lines. "
      "The tier and lines used are reported for every drug. Values are Level 5 z-scores.")
    w(f"\nScoring genes: the 400 primary maturation-signature genes that are L1000 landmark or best-inferred ({g['n_genes']} genes; {g['n_genes_landmark']} landmark-only for the sensitivity check).")

    w("\n## 2. Scores\n")
    w("- **Generic maturation score** (one number per drug): cosine between the drug signature and the maturation direction (each gene weighted by its "
      "Spearman rho with the differentiation score). Positive = the drug raises maturation genes and lowers stem/progenitor genes.")
    w("- **Patient reversion score:** cosine between the drug signature and the patient's push vector (minus the stuck signature from Phase 2). "
      "Positive = the drug moves expression in the direction this patient needs to reach maturity.")
    w(f"- **Connectivity (second method):** CMap KS connectivity score (Lamb et al. 2006) with the top {p3['connectivity_set_size']} genes the patient "
      "needs raised and lowered as query sets. Rank-based, so it ignores signature magnitude.")
    w(f"\nMethod agreement on patient rankings (median per-patient Spearman, cosine vs connectivity): van Galen {st['vangalen_method_agreement']:.2f}, "
      f"Beat AML {st['beataml_method_agreement']:.2f}.")

    w("\n## 3. Positive controls\n")
    w("Evaluated before any patient-level validation. Background: all compounds whose consensus comes from AML/monocytic lines.\n")
    w(md(st["controls"]))
    w("\n")
    w(md(pd.DataFrame(st["control_tests"]).T, index=True))
    w("\nAUROC is the probability that a random control scores above a random other compound (0.5 = chance). Controls without AML/monocytic-line "
      "signatures (for example enasidenib and ivosidenib, profiled only in non-haematopoietic lines) cannot be placed in this background and are listed as such.")
    w("\n**Top 25 compounds by generic maturation score** (AML/monocytic lines):\n")
    w(md(st["top_generic"][["cosine", "connectivity", "cells", "n_sigs", "mean_tas", "cosine_landmark"]], index=True))
    w("\n**Bottom 10** (push away from maturity):\n")
    w(md(st["bottom_generic"][["cosine", "connectivity", "cells", "n_sigs"]], index=True))

    w("\n## 4. Per-patient rankings\n")
    rv = st["ranked_vg"]
    for p, df in rv.groupby("patient"):
        top = df.nsmallest(8, "rank")[["rank", "drug", "cosine", "connectivity", "tier", "cells", "in_beataml_panel"]]
        w(f"\n**{p}**\n")
        w(md(top))
    w(f"\nFull rankings for all {rv.drug.nunique():,} compounds and {rv.patient.nunique()} patients: `data/processed/drug_rankings_vangalen.parquet`. "
      "Beat AML patient scores (panel drugs only) are in `reversion_cosine_beataml.parquet` and `reversion_connectivity_beataml.parquet`.")

    w("\n## 5. Choices and justifications\n")
    rows = [
        ("Score every compound profiled in AML/monocytic lines, not just the Beat AML panel.", "Positive controls are only meaningful against a large, realistic background, and the web app can show repurposing candidates outside the panel. Validation (Phase 4) still uses only panel drugs."),
        ("Average signatures within each line, then across lines.", "Prevents a line with many replicates from dominating, and cancels some line-specific noise."),
        ("Prefer AML/monocytic lines, then CD34/K562, then non-myeloid.", "Drug effects on expression depend on cell context; the tier is recorded so results can be stratified."),
        ("Restrict to landmark + best-inferred genes.", "The 'inferred' L1000 genes are the least accurate; landmark-only is the sensitivity check."),
        ("Cosine as primary score.", "Uses the magnitude and sign of every gene and is simple to interpret; KS connectivity is the established rank-based alternative, used as the second method."),
        ("Generic score weighted by Spearman rho.", "Genes most tightly linked to maturation count most."),
        ("Positive-control tiers fixed in config before scoring.", "Prevents choosing controls that happen to score well."),
        ("AUROC with label-permutation p-value.", "Rank-based and needs no assumptions about score distributions."),
    ]
    w(md(pd.DataFrame(rows, columns=["Choice", "Justification"])))
    w("\n## 6. Things to watch\n")
    w("- Drugs profiled in several lines get averaged, less noisy signatures, which can make their cosine larger in magnitude than single-line drugs.")
    w("- LINCS exposures are short (mostly 6-24 h); differentiation takes days, so early transcriptional changes are only a proxy.")
    w("- Patient rankings are highly correlated with each other: most of the signal is 'which drugs push toward maturity at all', not 'which drug for which patient'. Phase 4 tests whether the patient-specific part carries information.")
    w("\n## Figures\n")
    for f in F.values():
        w(f"![{f.stem}](figures/{f.name})")
    path(cfg, "reports", "phase3_drug_ranking_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
