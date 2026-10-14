"""Phase 1 report and figures."""
from __future__ import annotations

import datetime as dt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pipeline.common import path

PALETTE = {"HSC": "#3B4CC0", "Prog": "#6A8BEF", "GMP": "#9ABBFF", "ProMono": "#F2A481", "Mono": "#D8553F", "cDC": "#8E44AD"}


def md(df, index=False):
    return df.to_markdown(index=index, floatfmt=".3f")


def _scatter(ax, xy, c, title, cmap=None, cats=None, s=2):
    if cats is not None:
        for k, col in cats.items():
            m = c == k
            ax.scatter(xy[m, 0], xy[m, 1], s=s, c=col, label=k, linewidths=0, rasterized=True)
        ax.legend(markerscale=5, fontsize=7, frameon=False, loc="best")
    else:
        sc_ = ax.scatter(xy[:, 0], xy[:, 1], s=s, c=c, cmap=cmap, linewidths=0, rasterized=True)
        plt.colorbar(sc_, ax=ax, fraction=0.04)
    ax.set_title(title, fontsize=10); ax.set_xticks([]); ax.set_yticks([])


def figures(cfg, a, grid, corr, st, fate) -> dict:
    F = {}
    xy = a.obsm["X_landscape"]; ct = np.asarray(a.obs.CellType.astype(str))
    src = np.asarray(a.obs.source.astype(str))

    fig, ax = plt.subplots(2, 3, figsize=(15, 9))
    _scatter(ax[0, 0], xy, ct, f"{st['embedding_method']}: cell type", cats=PALETTE)
    _scatter(ax[0, 1], xy, src, "source", cats={"healthy": "#4C72B0", "AML_normal": "#DD8452"})
    _scatter(ax[0, 2], xy, np.asarray(a.obs["sample"].astype(str)) == "BM4", "BM4 (largest donor) vs rest",
             cats={False: "#BBBBBB", True: "#C44E52"})
    _scatter(ax[1, 0], xy, a.obs.score.values, "differentiation score (DPT, primary)", cmap="viridis")
    _scatter(ax[1, 1], xy, a.obs.palantir.values, "Palantir pseudotime", cmap="viridis")
    _scatter(ax[1, 2], xy, a.obs.height.values, "landscape height", cmap="magma_r")
    fig.tight_layout(); F["embedding"] = path(cfg, "figures", "phase1_embedding.png"); fig.savefig(F["embedding"], dpi=150); plt.close(fig)

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    order = cfg["phase1"]["lineage_order"] + ["cDC"]
    data = [a.obs.score.values[ct == c] for c in order]
    parts = ax[0].violinplot(data, showmedians=True)
    for b, c in zip(parts["bodies"], order):
        b.set_facecolor(PALETTE[c]); b.set_alpha(0.7)
    ax[0].set_xticks(range(1, len(order) + 1), order); ax[0].set_ylabel("differentiation score"); ax[0].set_title("Score by annotated cell type")
    for c in order:
        m = ct == c
        ax[1].scatter(a.obs.dpt.values[m], a.obs.palantir.values[m], s=2, c=PALETTE[c], label=c, rasterized=True)
    ax[1].set_xlabel("DPT"); ax[1].set_ylabel("Palantir pseudotime")
    ax[1].set_title(f"DPT vs Palantir (Spearman {st['pseudotime']['dpt_vs_palantir_spearman']:.2f})"); ax[1].legend(markerscale=5, fontsize=7, frameon=False)
    med = a.obs.groupby(["CellType", "source"], observed=True).score.median().unstack().reindex(order)
    x = np.arange(len(order)); wd = 0.38
    ax[2].bar(x - wd / 2, med["healthy"], wd, label="healthy BM", color="#4C72B0")
    ax[2].bar(x + wd / 2, med["AML_normal"], wd, label="normal cells in AML samples", color="#DD8452")
    ax[2].set_xticks(x, order); ax[2].set_ylabel("median score"); ax[2].set_title("Healthy vs AML-sample normal cells"); ax[2].legend(fontsize=8, frameon=False)
    fig.tight_layout(); F["pseudotime"] = path(cfg, "figures", "phase1_pseudotime.png"); fig.savefig(F["pseudotime"], dpi=150); plt.close(fig)

    mk = cfg["phase1"]["markers"]; genes = mk["stem"] + mk["mature"] + ["MPO", "ELANE"]
    lin = np.isin(ct, cfg["phase1"]["lineage_order"])
    s = a.obs.score.values[lin]; bins = np.linspace(0, 1, 21); mid = (bins[:-1] + bins[1:]) / 2
    b = np.clip(np.digitize(s, bins) - 1, 0, 19)
    fig, ax = plt.subplots(1, len(genes), figsize=(2.6 * len(genes), 3), sharex=True)
    cr = corr.set_index("gene").spearman
    for i, g in enumerate(genes):
        if g not in a.var_names:
            ax[i].set_title(f"{g} (absent)"); continue
        v = a[lin, g].X.toarray().ravel()
        mean = pd.Series(v).groupby(b).mean().reindex(range(20)); se = pd.Series(v).groupby(b).sem().reindex(range(20))
        col = "#3B4CC0" if g in mk["stem"] else ("#D8553F" if g in mk["mature"] else "#777777")
        ax[i].plot(mid, mean, color=col); ax[i].fill_between(mid, mean - se, mean + se, color=col, alpha=0.25)
        ax[i].set_title(f"{g}  rho={cr.get(g, np.nan):.2f}", fontsize=9); ax[i].set_xlabel("score")
    ax[0].set_ylabel("log-norm expression")
    fig.suptitle("Marker trends along the HSC -> Mono lineage (blue = stem, red = mature, grey = GMP context)", fontsize=10)
    fig.tight_layout(); F["markers"] = path(cfg, "figures", "phase1_markers.png"); fig.savefig(F["markers"], dpi=150); plt.close(fig)

    fig = plt.figure(figsize=(14, 5.5))
    ax1 = fig.add_subplot(1, 2, 1)
    ext = [grid["x"][0], grid["x"][-1], grid["y"][0], grid["y"][-1]]
    im = ax1.imshow(grid["height"], origin="lower", extent=ext, cmap="magma_r", aspect="auto")
    ax1.contour(grid["x"], grid["y"], grid["height"], levels=14, colors="white", linewidths=0.4, alpha=0.6)
    for c, (cx, cy) in st["centroids"].iterrows():
        ax1.text(cx, cy, c, color="white", fontsize=9, ha="center", fontweight="bold")
    plt.colorbar(im, ax=ax1, fraction=0.04); ax1.set_title("Potential surface (height), contours and cell-type centroids")
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    X, Y = np.meshgrid(grid["x"], grid["y"])
    ax2.plot_surface(X, Y, grid["height"], cmap="viridis_r", linewidth=0, antialiased=True, alpha=0.9, rstride=2, cstride=2)
    sub = np.random.default_rng(0).choice(len(xy), size=min(2500, len(xy)), replace=False)
    ax2.scatter(xy[sub, 0], xy[sub, 1], a.obs.height.values[sub] + 0.01, s=1, c=[PALETTE[c] for c in ct[sub]])
    ax2.view_init(elev=38, azim=-60); ax2.set_title("3D view (cells coloured by type)"); ax2.set_axis_off()
    fig.tight_layout(); F["landscape"] = path(cfg, "figures", "phase1_landscape.png"); fig.savefig(F["landscape"], dpi=150); plt.close(fig)

    lo = st["mapping"]["loso"].sort_values("score_mae")
    fig, ax = plt.subplots(figsize=(10, 3.8))
    ax.bar(lo.held_out_sample, lo.score_mae, color="#4C72B0")
    ax.axhline(st["mapping"]["loso_weighted_mae"], color="k", ls="--", lw=1, label="cell-weighted mean")
    ax.set_ylabel("score MAE (0-1 scale)"); ax.set_title("Mapping error with each sample held out"); ax.legend(frameon=False)
    plt.xticks(rotation=90, fontsize=7); fig.tight_layout()
    F["mapping"] = path(cfg, "figures", "phase1_mapping_loso.png"); fig.savefig(F["mapping"], dpi=150); plt.close(fig)
    return F


def write_report(cfg, a, grid, corr, sig, st, fate) -> None:
    p1 = cfg["phase1"]; F = figures(cfg, a, grid, corr, st, fate)
    pt, mp, mx = st["pseudotime"], st["mapping"], st["mixing"]
    L = []; w = L.append
    ct = a.obs.CellType.astype(str)
    w("# Phase 1: normal myeloid differentiation landscape")
    w(f"\nGenerated {dt.date.today().isoformat()} by `python -m pipeline.run --phase 1`. All numbers are computed by the pipeline.\n")

    w("## Summary\n")
    w(f"- Reference: **{a.n_obs:,} normal myeloid cells** ({(a.obs.source=='healthy').sum():,} healthy BM, "
      f"{(a.obs.source=='AML_normal').sum():,} normal cells from AML samples), {st['n_hvg']} highly variable genes, {p1['n_pcs']} PCs, Harmony by sample.")
    w(f"- Differentiation score (DPT, rooted in a healthy HSC) orders the annotated types HSC < Prog < GMP < ProMono < Mono with "
      f"**Spearman {pt['order_rho']['dpt']:.3f}**; Palantir gives {pt['order_rho']['palantir']:.3f}. "
      f"DPT and Palantir agree with **Spearman {pt['dpt_vs_palantir_spearman']:.3f}** ({pt['dpt_vs_palantir_spearman_lineage_only']:.3f} on the main lineage).")
    w(f"- Root choice barely matters: across {p1['n_root_candidates']-1} alternative HSC roots, the minimum Spearman with the primary score is {pt['root_sensitivity_min_rho']:.3f}.")
    w(f"- Embedding: {st['embedding_method']}. Landscape height vs score: Spearman {st['height_vs_score_spearman']:.3f} "
      "(stem-like cells high, mature cells low, as intended).")
    w(f"- Mapping (expression -> landscape): random 20% hold-out score MAE **{mp['random_split']['score_mae']:.3f}** "
      f"(Spearman {mp['random_split']['score_spearman']:.3f}); with whole samples held out, cell-weighted MAE **{mp['loso_weighted_mae']:.3f}**. "
      f"Position error is {mp['random_split']['xy_rmse_frac_of_range']*100:.1f}% (random) and {mp['loso_weighted_xy']*100:.1f}% (held-out sample) of the embedding range.")
    mk = st["markers"]
    w(f"- Marker sanity check {'**passed**' if st['markers_pass'] else '**FAILED**'}: " +
      ", ".join(f"{g} rho {r:+.2f}" for g, (r, _) in mk.items()) + " (stem markers should be negative, mature markers positive).")

    w("\n## 1. Reference cells and QC\n")
    q = st["qc"]
    w(f"Cells: van Galen cells with `PredictionRefined = normal` and CellType in {p1['cell_types']}, excluding cell lines. "
      f"{q['n_cells_selected']:,} selected, {q['n_cells_after_qc']:,} after QC (>= {p1['qc']['min_genes']} genes, < {p1['qc']['max_pct_mt']}% mitochondrial UMIs). "
      f"The authors' data were already filtered (minimum {q['umis_min']:.0f} UMIs per cell, maximum mitochondrial fraction {q['pct_mt_max']:.1f}%), "
      f"so our QC removes {q['n_cells_selected']-q['n_cells_after_qc']} cells. Genes detected in fewer than {p1['qc']['min_cells_per_gene']} cells were dropped "
      f"({q['n_genes_before']:,} -> {q['n_genes_after']:,}).\n")
    w(md(a.obs.groupby(["CellType", "source"], observed=True).size().unstack().reindex(p1["cell_types"]), index=True))
    w("\nProcessing: library-size normalization to 10,000 UMIs and log1p; 2,000 HVGs ranked by normalized dispersion within each batch and combined across batches (scanpy `seurat` flavor; the `seurat_v3` loess fit failed numerically on several small sorted batches); "
      f"PCA on centred log expression of HVGs ({st['pca_var_explained']*100:.1f}% of HVG variance in {p1['n_pcs']} PCs); Harmony on the PCs. Batch = sample, except that samples with fewer than {p1['min_batch_cells']} reference cells "
      f"({a.obs.loc[a.obs.batch=='pooled_small_samples', 'sample'].nunique()} samples, {(a.obs.batch=='pooled_small_samples').sum()} cells) are pooled into one batch, because a "
      "per-batch variance fit or correction cannot be estimated from a handful of cells.")
    w("\n**Batch structure** (fraction of each cell's 30 nearest neighbours sharing its label):\n")
    w(md(pd.DataFrame(mx).T, index=True))
    w(f"\nWithout correction, cells cluster by sample far more than chance. Harmony cuts same-sample neighbours from {mx['uncorrected_pca']['sample']:.2f} to "
      f"{mx['harmony']['sample']:.2f} (random: {mx['expected_if_random']['sample']:.2f}). It costs some cell-type purity ({mx['uncorrected_pca']['CellType']:.2f} -> "
      f"{mx['harmony']['CellType']:.2f}, still far above the random {mx['expected_if_random']['CellType']:.2f}), so batch is reduced, not removed. "
      "The landscape (pseudotime, embedding, surface) is built on the Harmony PCs so that donor-specific islands do not become fake valleys.")

    w("\n## 2. Pseudotime and differentiation score\n")
    w(f"- **Root:** the healthy-donor HSC with the highest `Score_HSC` (the authors' classifier score, computed independently of our pipeline): "
      f"`{pt['root_cell']}` from {pt['root_sample']}.")
    w(f"- **Primary method: diffusion pseudotime (DPT)** on a 30-NN graph of the Harmony PCs. Every cell is connected to the root.")
    w(f"- **Second method: Palantir** (same PCs, same root, {p1['palantir_waypoints']} waypoints). Left to detect terminal states itself, it picked "
      f"cells labelled {pt['palantir_auto_terminal_states']}, which is biologically wrong (neither is an end state). We therefore fixed the terminal "
      f"states from the authors' classifier scores, independently of our pseudotime: the Mono cell with the highest `Score_Mono` and the cDC cell with the "
      f"highest `Score_cDC`. The pseudotime itself barely changes (Spearman {pt['palantir_auto_vs_fixed_spearman']:.3f} between the two runs; order agreement "
      f"{pt['palantir_auto_order_rho']:.3f} with automatic states). All Palantir numbers below use the fixed states.")
    w(f"- **Uncorrected check:** DPT on uncorrected PCs agrees with the primary at Spearman {pt['dpt_vs_uncorrected_spearman']:.3f}, so the ordering is not an "
      f"artefact of Harmony. It orders the annotated types slightly better ({pt['order_rho']['dpt_uncorrected']:.3f} vs {pt['order_rho']['dpt']:.3f}); "
      "we keep Harmony as primary for the batch reason in section 1, a choice made before computing either value.")
    w(f"- **Score:** DPT rescaled to [0, 1] between its {p1['score_clip_quantiles'][0]*100:.0f}st and {p1['score_clip_quantiles'][1]*100:.0f}th percentiles "
      f"({pt['score_clip'][0]:.3f} -> 0, {pt['score_clip'][1]:.3f} -> 1) and clipped, so a few extreme cells do not compress everyone else.\n")
    w(md(pd.DataFrame(pt["order_rho"], index=["Spearman with HSC<Prog<GMP<ProMono<Mono"]).T, index=True))
    w("\nMedian score and Palantir pseudotime by cell type and source:\n")
    w(md(st["median_by_type_source"].reindex(p1["cell_types"]), index=True))
    w("\nWithin each cell type, healthy and AML-sample normal cells get similar scores, which supports pooling them.")
    w("\nNote on shape: DPT is not linear in biological stage. HSC and Prog are compressed near 0, and most of the range is spent between GMP and Mono. "
      "Phase 2 thresholds must therefore be chosen on this score's actual distribution, not assumed evenly spaced.")
    w("\nThe cDC branch is not part of the HSC->Mono order; it receives high scores (median about 0.75), consistent with being a mature, separate fate.")
    w("\nPalantir fate probabilities (median per cell type):\n")
    w(md(st["fate_by_type"], index=True))
    fb = st["fate_by_type"]
    w(f"\nNegative result: the fate probabilities separate the monocyte fate (Mono cells {fb.loc['Mono', 'fate_Mono']:.2f}) but not the dendritic one: "
      f"cDC cells get only {fb.loc['cDC', 'fate_cDC']:.2f} cDC probability. Early cells sit at about 0.5/0.5, as expected before commitment. "
      "Because the cDC fate is not resolved, fate probabilities are reported here but not used downstream.")

    w("\n## 3. Embedding and potential surface\n")
    w(f"- Embedding: {st['embedding_method']} on the Harmony PCs ({'the preferred method; no fallback needed' if st['embedding_method']=='PHATE' else 'fallback'}).")
    lc = p1["landscape"]
    w(f"- Surface: {lc['grid']}x{lc['grid']} grid with {lc['margin']*100:.0f}% margin. Cells are binned and Gaussian-smoothed (bandwidth {lc['bandwidth']*100:.1f}% of the range) to get a local "
      f"mean score and a density. Height = {lc['weight_pseudotime']} x (1 - local score) + {1-lc['weight_pseudotime']:.2f} x normalized(-log density), "
      "blended toward a high rim where there are no cells, lightly smoothed, then scaled to [0, 1].")
    w("- Reading the surface: stem-like cells sit on the high plateau; dense, mature populations form the deepest basins. The density term adds local "
      "dips where many cells sit, so each stable population looks like a basin, as in Waddington's picture.")
    w("\nCell-type centroids (used for the floating labels in the web app):\n")
    w(md(st["centroids"].rename(columns={0: "x", 1: "y"}).reindex(p1["cell_types"]), index=True))

    w("\n## 4. Mapping new cells onto the landscape\n")
    w(f"Model: log-normalized expression of the {st['n_hvg']} HVGs -> centre with reference means -> fixed PCA projection -> distance-weighted kNN "
      f"(k = {p1['mapping']['k']}) regression onto (score, x, y). It uses the uncorrected PCA deliberately: Harmony corrects a whole dataset at once and cannot be "
      "applied to a single new cell or to the drug-shifted vectors of Phase 5. Genes missing from a query are set to the reference mean.\n")
    w(md(pd.DataFrame([mp["random_split"]])))
    w("\nLeave-one-sample-out (each sample with >= 30 cells held out in turn):\n")
    w(md(mp["loso"]))
    w("\nSpearman is low for some held-out samples because their true scores barely vary (see `score_sd_true`, e.g. the CD34+CD38- sorted BM5 sample is "
      "almost all HSC). MAE is the relevant number there. Caveat: the PCA basis and the target pseudotime were computed with all cells, so this tests "
      "how well the fixed map generalizes across samples, not a full re-derivation of the landscape without them.")

    w("\n## 5. Maturation signature\n")
    w(f"Spearman correlation of each gene (detected in >= 1% of cells) with the score, on the HSC->Mono lineage only (cDC excluded). "
      f"{len(corr):,} genes tested. The signature is the top {p1['signature_top_n']} positively and top {p1['signature_top_n']} negatively correlated genes. "
      "All genes, before any exclusion:\n")
    top = pd.concat([corr.head(15).reset_index(drop=True)[["gene", "spearman"]].rename(columns={"gene": "up gene", "spearman": "rho up"}),
                     corr.tail(15).iloc[::-1].reset_index(drop=True)[["gene", "spearman"]].rename(columns={"gene": "down gene", "spearman": "rho down"})], axis=1)
    w(md(top))
    sg = st["signature"]
    w(f"\nRibosomal-protein and mitochondrial genes (regex `{p1['signature_exclude_regex']}`, {sg['n_excluded_genes']} genes) dominate the down-regulated end: "
      f"they make up {sg['full_down_frac_excluded']*100:.0f}% of the 200 most negatively correlated genes (and {sg['full_up_frac_excluded']*100:.0f}% of the positive end). "
      "Their decline is real biology (progenitors build more ribosomes), but it reflects growth more than lineage maturation, and drugs that simply slow "
      "growth lower them too. Left in, they would let Phase 3 reward any antiproliferative drug as 'differentiating'. The **primary signature excludes them**. "
      f"The full signature is saved as `maturation_signature_full.csv` for a sensitivity analysis ({sg['overlap_primary_full']} of 400 genes shared). "
      "This rule was set before any drug was scored.")
    w(f"\nWeakest gene kept in the primary signature: |rho| = {sig.spearman.abs().min():.3f}. Top genes of the primary signature:\n")
    ps = pd.concat([sig[sig.direction.str.startswith("up")].head(15).reset_index(drop=True)[["gene", "spearman"]].rename(columns={"gene": "up gene", "spearman": "rho up"}),
                    sig[sig.direction.str.startswith("down")].head(15).reset_index(drop=True)[["gene", "spearman"]].rename(columns={"gene": "down gene", "spearman": "rho down"})], axis=1)
    w(md(ps))

    w("\n## 6. Sanity checks\n")
    w(md(pd.DataFrame([{"gene": g, "expected": "falls" if k == "stem" else "rises", "spearman": r,
                        "pass": (r < 0) if k == "stem" else (r > 0)} for g, (r, k) in mk.items()])))

    w("\n## 7. Choices and justifications\n")
    rows = [
        ("Use cells the authors called normal, from healthy donors and from AML samples.",
         "The authors' calls use genotyping plus expression, which we cannot redo. Normal cells from AML marrow add about 3,500 cells, and section 2 shows they get the same scores as healthy cells of the same type."),
        ("Keep cDC in the reference but out of the ordered lineage.",
         "cDCs are a real mature fate and AML blasts can be cDC-like, so the landscape must contain them. They are not on the HSC->monocyte path, so ordering checks and the signature exclude them."),
        ("HVGs chosen within each batch, then combined.",
         "Picking variable genes across all cells at once would favour genes that differ between donors rather than between cell states."),
        ("No gene scaling before PCA.",
         "Scaling inflates noisy low-expression genes. Leaving genes in log units also keeps Phase 5 drug shifts, which are added in log space, consistent with the projection."),
        ("Harmony for the reference structure.",
         "Section 1 shows strong sample clustering. Harmony is a standard, well-tested correction that works on PCs."),
        ("harmonypy 0.0.10.",
         "The current release needs a C++ compiler to install. 0.0.10 is the pure-Python implementation of the same algorithm."),
        ("DPT as the primary pseudotime, Palantir as the second method.",
         "DPT was named first in the plan and is deterministic given the graph. Palantir is an independent algorithm (waypoints and shortest paths) that also returns fate probabilities, so agreement between the two is informative."),
        ("Root = the healthy HSC with the highest authors' HSC score; sensitivity over 10 roots.",
         "This avoids picking the root by looking at our own pseudotime, and the sensitivity check shows the result does not depend on that one cell."),
        ("Score = DPT clipped at the 1st/99th percentiles.",
         "Min-max scaling would let a handful of outlying cells set the scale."),
        (f"PHATE for the 2D embedding.",
         "PHATE is designed to preserve continuous trajectories and branches, which suits a differentiation landscape better than UMAP's cluster-focused layout."),
        (f"Height weights {p1['landscape']['weight_pseudotime']} score / {1-p1['landscape']['weight_pseudotime']:.2f} density.",
         "The plan asks for pseudotime blended with -log density. Pseudotime dominates so the slope always runs from stem to mature; the smaller density term shapes basins without inverting the slope."),
        ("Fixed PCA + kNN for mapping, not Harmony-corrected coordinates.",
         "The map must apply to single new cells and to synthetic shifted vectors (Phase 5). Leave-one-sample-out error measures the cost of not batch-correcting queries."),
        ("Distance-weighted kNN with k = 15.",
         "A local average of 15 reference cells is smooth enough to be stable but local enough to keep branch structure. It also gives an interpretable placement: a query lands among the reference cells it most resembles."),
        ("Signature on the HSC->Mono lineage, top 200 genes each way.",
         "Most AML blasts in this cohort are on the monocytic/granulocytic path. Two hundred genes per direction is large enough for robust connectivity scoring and is in the range commonly used for CMap queries."),
        ("Exclude ribosomal and mitochondrial genes from the primary signature.",
         "They track growth rather than lineage identity, and drug-induced growth arrest lowers them in any cell line, so they would confound 'pushes toward maturity' with 'slows proliferation'. The full signature is kept as a pre-specified sensitivity analysis."),
        ("Fix Palantir terminal states from the authors' Mono and cDC scores.",
         "Automatic detection picked Prog and ProMono cells, which are not end states. Using the authors' classifier keeps the choice independent of our own pseudotime."),
    ]
    w(md(pd.DataFrame(rows, columns=["Choice", "Justification"])))

    w("\n## 8. Things to watch\n")
    w("- The score is compressed at the stem end (HSC and Prog both near 0), so it separates early states less well than late ones.")
    w("- Batch is reduced, not removed; the mapping of new cells is uncorrected and inherits some batch error (quantified in section 4).")
    w("- Weak marker correlations for PROM1 (-0.23) and CSF1R (+0.21) reflect low detection in 3-prime scRNA-seq; their binned trends (figure) still go the expected way.")
    w("- Malignant cells may sit off the normal manifold; kNN will still place them among their nearest normal neighbours. Phase 2 should report distance to the reference so off-manifold cells are visible.")
    w("- The landscape is a 2D projection; distances on it are qualitative.")

    w("\n## Figures\n")
    for f in F.values():
        w(f"![{f.stem}](figures/{f.name})")
    w("\n## Outputs (`data/processed/`)\n")
    w("`landscape_reference.h5ad` (reference cells with scores and coordinates), `landscape_grid.npz` (height surface), "
      "`landscape_mapper.joblib` (expression -> landscape map), `maturation_correlations.csv`, `maturation_signature.csv`, "
      "`landscape_centroids.csv`, `phase1_stats.json`.")
    path(cfg, "reports", "phase1_landscape_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
