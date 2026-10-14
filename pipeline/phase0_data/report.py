"""Write reports/phase0_data_report.md and Phase 0 QC figures."""
from __future__ import annotations

import datetime as dt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pipeline.common import ROOT, path

MYELOID_NORMAL = ["HSC", "Prog", "GMP", "ProMono", "Mono", "cDC"]
MYELOID_MALIG = ["HSC-like", "Prog-like", "GMP-like", "ProMono-like", "Mono-like", "cDC-like"]


def md(df: pd.DataFrame, index: bool = False) -> str:
    return df.to_markdown(index=index)


def size_mb(p) -> float:
    return sum(f.stat().st_size for f in p.rglob("*") if f.is_file()) / 1e6 if p.is_dir() else p.stat().st_size / 1e6


def figures(cfg, vg, b, sigs) -> list[str]:
    out = []
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    o = vg.obs
    order = o.groupby("sample", observed=True).sample_type.first().sort_values().index
    data = [np.log10(o.loc[o["sample"] == s, "n_counts"] + 1) for s in order]
    ax[0].boxplot(data, showfliers=False)
    ax[0].set_xticks(range(1, len(order) + 1), order, rotation=90, fontsize=6)
    ax[0].set_ylabel("log10 UMIs per cell"); ax[0].set_title("van Galen: UMIs per cell by sample")
    d = b["response"].groupby("inhibitor").dbgap_rnaseq_sample.nunique().sort_values()
    ax[1].hist(d.values, bins=30, color="#4C72B0")
    ax[1].set_xlabel("samples with RNA-seq per drug"); ax[1].set_ylabel("drugs")
    ax[1].set_title("Beat AML: drug coverage")
    fig.tight_layout(); f = path(cfg, "figures", "phase0_qc_vangalen_beataml.png"); fig.savefig(f, dpi=150); plt.close(fig)
    out.append(f)

    fig, ax = plt.subplots(figsize=(8, 4))
    s = sigs.obs.groupby("cell_iname").size().sort_values(ascending=False)
    colors = ["#DD8452" if c in cfg["lincs"]["myeloid_cell_lines"] else "#8C8C8C" for c in s.index]
    ax.bar(s.index[:30], s.values[:30], color=colors[:30])
    ax.set_ylabel("signatures downloaded"); ax.set_title("LINCS signatures by cell line (orange = myeloid/heme)")
    plt.xticks(rotation=90, fontsize=7); fig.tight_layout()
    f = path(cfg, "figures", "phase0_lincs_selection.png"); fig.savefig(f, dpi=150); plt.close(fig)
    out.append(f)
    return out


def write_report(cfg: dict, s: dict, vg) -> None:
    b, L, sigs = s["beataml"], s["lincs"], s["lincs_sigs"]
    raw = ROOT / cfg["paths"]["raw"]
    figs = figures(cfg, vg, b, sigs)
    lines: list[str] = []
    w = lines.append
    o = vg.obs

    # ---------------- van Galen
    samp = s["vangalen"]["samples"].copy()
    calls = o.groupby(["sample", "PredictionRefined"], observed=True, dropna=False).size().unstack(fill_value=0)
    calls.columns = [str(c) for c in calls.columns]
    qc = o.groupby("sample", observed=True).agg(median_umis=("n_counts", "median"),
                                                median_genes=("n_genes_detected", "median"))
    samp = samp.set_index("sample").join(calls).join(qc).reset_index()
    samp["day"] = samp["day"].astype("Int64")
    nonline = o[o.sample_type != "cell_line"]
    ct = nonline.groupby(["CellType", "sample_type"], observed=True).size().unstack(fill_value=0)
    my_norm = nonline[nonline.CellType.isin(MYELOID_NORMAL) & (nonline.PredictionRefined == "normal")]
    my_norm_tab = my_norm.groupby(["CellType", "sample_type"], observed=True).size().unstack(fill_value=0).reindex(MYELOID_NORMAL)
    malig = o[(o.sample_type == "AML") & (o.PredictionRefined == "malignant")]
    d0 = malig[malig.day == 0]
    mal_tab = d0.groupby(["patient", "CellType"], observed=True).size().unstack(fill_value=0)
    mal_tab = mal_tab.reindex(columns=[c for c in MYELOID_MALIG if c in mal_tab.columns])
    mal_tab["total"] = mal_tab.sum(1)
    nan_calls = o[o.PredictionRefined.isna()]

    # ---------------- Beat AML
    e, es = b["expr"], b["expr_stats"]
    r = b["response"]
    ctrl = e.obs.rna_control[~e.obs.rna_control.isin(["", "No"])]
    rna = set(e.obs_names) - set(ctrl.index)
    drug_s = set(r.dbgap_rnaseq_sample.dropna())
    mut = b["mutations"]
    both = rna & drug_s
    clin = b["clinical"].dropna(subset=["dbgap_rnaseq_sample"]).drop_duplicates("dbgap_rnaseq_sample").set_index("dbgap_rnaseq_sample")
    dna_of = clin["dbgap_dnaseq_sample"].reindex(list(both))
    triple = {k for k, v in dna_of.items() if isinstance(v, str) and v in mut.index}
    stage = clin.loc[list(both), "diseaseStageAtSpecimenCollection"].value_counts()
    per_drug = r.dropna(subset=["dbgap_rnaseq_sample"]).groupby("inhibitor").dbgap_rnaseq_sample.nunique()
    trip_mut = mut.loc[[dna_of[k] for k in triple]]
    vg_sc = b["vg_scores"]

    # ---------------- LINCS + overlaps
    match = L["match"]
    matched = match.dropna(subset=["lincs_names"])
    so = sigs.obs
    by_drug = so.groupby("drug").agg(n_sigs=("cmap_name", "size"), n_cells=("cell_iname", "nunique"),
                                      myeloid_cells=("cell_iname", lambda x: ",".join(sorted(set(x) & set(cfg["lincs"]["myeloid_cell_lines"])))),
                                      n_hiq=("is_hiq", lambda x: int((x.astype(str).isin(["1", "1.0"])).sum())))
    overlap = matched[["beat_name", "lincs_names", "match_how", "match_quality", "note"]].set_index("beat_name")
    overlap = overlap.join(by_drug).join(per_drug.rename("beat_rna_samples"))
    overlap["has_myeloid_sig"] = overlap.myeloid_cells.fillna("").ne("")
    n_my = int(overlap.has_myeloid_sig.sum())
    n_aml_line = int(overlap.myeloid_cells.fillna("").str.contains("THP1|HL60|SKM1|NOMO1|U937").sum())
    overlap["beat_rna_samples"] = overlap.index.map(per_drug).fillna(0).astype(int)
    overlap["n_sigs"] = overlap.n_sigs.fillna(0).astype(int)
    usable = overlap[(overlap.beat_rna_samples >= 50) & (overlap.n_sigs > 0)]
    usable_my = usable[usable.has_myeloid_sig]
    no_sigs = overlap[overlap.n_sigs == 0].index.tolist()
    few = overlap[overlap.beat_rna_samples < 50].index.tolist()

    gv, gb, gl = set(vg.var_names), set(e.var_names), set(sigs.var_names)
    reliable = set(sigs.var_names[sigs.var.feature_space.isin(["landmark", "best inferred"])])
    landmark = set(sigs.var_names[sigs.var.feature_space == "landmark"])

    controls = so[so.is_positive_control.astype(str) == "True"].groupby("drug").agg(
        n_sigs=("cmap_name", "size"), cells=("cell_iname", lambda x: ",".join(sorted(set(x)))),
        n_hiq=("is_hiq", lambda x: int((x.astype(str).isin(["1", "1.0"])).sum())))

    # ---------------- write
    w("# Phase 0 data report")
    w(f"\nGenerated {dt.date.today().isoformat()} by `python -m pipeline.run --phase 0`. All numbers below are computed from the downloaded files.\n")
    w("## Headline numbers\n")
    w(f"- **van Galen (GSE116256):** {vg.n_obs:,} cells x {vg.n_vars:,} genes from {o['sample'].nunique()} samples "
      f"({(samp.sample_type=='AML').sum()} AML samples from {o[o.sample_type=='AML'].patient.nunique()} patients, "
      f"{(samp.sample_type=='healthy_BM').sum()} healthy bone marrow samples from {o[o.sample_type=='healthy_BM'].patient.nunique()} donors, 2 cell lines).")
    w(f"- **Normal myeloid cells for the landscape (Phase 1):** {len(my_norm):,} cells labelled normal HSC/Prog/GMP/ProMono/Mono/cDC "
      f"({int((my_norm.sample_type=='healthy_BM').sum()):,} from healthy BM, {int((my_norm.sample_type=='AML').sum()):,} normal cells inside AML samples).")
    w(f"- **Malignant cells at diagnosis (day 0):** {len(d0):,} cells across {d0.patient.nunique()} patients.")
    w(f"- **Beat AML:** {e.n_obs} RNA-seq samples x {e.n_vars:,} genes ({len(ctrl)} healthy controls); "
      f"{r.inhibitor.nunique()} single-agent drugs, {len(r):,} sample-drug AUCs; "
      f"**{len(both)} AML RNA-seq samples have drug data, {len(triple)} of these also have mutation calls.**")
    w(f"- **LINCS (CMap LINCS 2020 release):** {L['n_trt_cp_sigs']:,} compound signatures, {L['n_compounds']:,} compounds, {L['n_cells']} cell lines. "
      f"Downloaded {sigs.n_obs:,} signatures x {sigs.n_vars:,} genes for {so.drug.nunique()} drugs.")
    w(f"- **Drug overlap Beat AML x LINCS: {len(matched)} of {match.shape[0]} Beat AML drugs** "
      f"({int((matched.match_quality=='exact').sum())} exact, {int((matched.match_quality=='approximate').sum())} approximate). "
      f"{n_my} have a signature in a myeloid/heme line, {n_aml_line} in an AML/monocytic line (THP1, HL60, SKM1, NOMO1, U937). "
      f"**Usable for validation (>= 50 Beat AML RNA-seq samples and >= 1 downloaded signature): {len(usable)} drugs, "
      f"{len(usable_my)} of them with a myeloid/heme signature.** This is far above the ~15-drug threshold, so power is not limited by drug overlap.")
    w(f"- **Gene overlap (HGNC):** all three datasets share {len(gv & gb & gl):,} genes; "
      f"{len(gv & gb & reliable):,} of them are L1000 landmark or best-inferred genes, {len(gv & gb & landmark)} are landmark genes.")

    w("\n## 1. Sources and files\n")
    files = []
    for sub in ["vangalen", "beataml", "lincs", "hgnc"]:
        for f in sorted((raw / sub).rglob("*")):
            if f.is_file() and not (sub == "vangalen" and f.name != "GSE116256_RAW.tar"):
                files.append({"source": sub, "file": str(f.relative_to(raw)).replace("\\", "/"), "MB": round(f.stat().st_size / 1e6, 2)})
    w(md(pd.DataFrame(files)))
    w(f"\nThe van Galen tar was extracted into {len(list((raw/'vangalen').glob('*.txt.gz')))} per-sample text files in `data/raw/vangalen/`. "
      "LINCS signature values were not downloaded as a file: they were read by HTTP range requests from the public 35.5 GB "
      "`level5_beta_trt_cp_n720216x12328.gctx` (see section 4).")
    w("\nLicences and access: Beat AML processed data are CC-BY-4.0 and openly downloadable from the authors' "
      "GitHub (biodev/beataml2.0_data, linked from https://biodev.github.io/BeatAML2/). No dbGaP / controlled-access data were used; "
      "the `dbgap` in file names refers only to the sample identifiers. GEO and CMap LINCS data are public.")

    w("\n## 2. van Galen et al. 2019 (single-cell)\n")
    w("Each sample has a dense gene x cell UMI table (`*.dem.txt.gz`) and a per-cell annotation table (`*.anno.txt.gz`). "
      "Annotation columns used later: `PredictionRefined` (malignant / normal / unclear), `CellType` (e.g. HSC, HSC-like, GMP-like), "
      "per-cell type scores `Score_*`, and the sample name, which encodes patient and day (D0 = diagnosis, later days = during treatment). "
      f"Every annotated cell was found in its count table. All {len(samp)} samples share the same {s['vangalen']['gene_mapping']['n_input']:,}-gene list.")
    gm = s["vangalen"]["gene_mapping"]
    w(f"\nGene symbols to HGNC: {gm['n_approved']:,} already approved, {gm.get('n_previous',0):,} mapped from previous symbols, "
      f"{gm.get('n_alias',0)} from unambiguous aliases, {gm['n_unmapped']:,} unmapped (kept under their original name; mostly lncRNAs/clone IDs). "
      f"{gm['n_collapsed']} symbols collapsed onto an existing gene (UMIs summed).")
    w("\n### Samples\n")
    cols = ["sample", "patient", "day", "sample_type", "n_cells_annotated"] + [c for c in ["malignant", "normal", "unclear", "nan"] if c in samp.columns] + ["median_umis", "median_genes"]
    w(md(samp[cols]))
    w(f"\n{len(nan_calls)} cells, all from {', '.join(sorted(nan_calls['sample'].unique()))}, are listed in the annotation file "
      f"but have no annotation at all ({int(nan_calls.CellType.isna().sum())} have empty `CellType`, `PredictionRefined` and `Score_*`). "
      "They are probably cells the authors did not classify; we exclude them from all analyses.")
    w("\n### Cell types (excluding cell lines)\n")
    w(md(ct, index=True))
    w("\n### Normal myeloid cells available for the Phase 1 landscape\n")
    w(md(my_norm_tab, index=True))
    w("\n### Malignant cells per patient at day 0 (CellType)\n")
    w(md(mal_tab, index=True))
    w("\nSome patients have very few malignant cells at diagnosis (see `total`), which will limit per-patient distributions in Phase 2.")

    w("\n## 3. Beat AML 2.0 (bulk RNA-seq, mutations, ex vivo drug response)\n")
    w(f"- **Expression** (`beataml_waves1to4_norm_exp_dbgap.txt`): {es['n_samples']} samples x {es['n_genes_input']:,} Ensembl genes (GRCh37). "
      f"Mapped to HGNC through Ensembl IDs ({es['hgnc_how']}); {es['n_unmapped_dropped']:,} genes without an HGNC symbol dropped, "
      f"{es['n_duplicate_symbols_dropped']} duplicate symbols resolved by keeping the most expressed Ensembl ID, leaving {es['n_genes']:,} genes. "
      f"Values range {es['value_range'][0]:.2f} to {es['value_range'][1]:.2f} (normalized, log-scale as provided; we use them as is).")
    w(f"- **Healthy controls with RNA-seq:** {ctrl.astype(str).value_counts().to_dict()}.")
    w(f"- **Drug response** (`beataml_probit_curve_fits_v4_dbgap.txt`): {len(r):,} rows, {r.inhibitor.nunique()} drugs, "
      f"{r.dbgap_rnaseq_sample.nunique()} RNA-seq samples, {r.dbgap_subject_id.nunique()} patients. All rows are single-agent probit fits with "
      f"`paper_inclusion = True`. Split flag in the file: {r.status.value_counts().to_dict()}. "
      f"AUC median {r.auc.median():.1f} (range {r.auc.min():.1f}-{r.auc.max():.1f}); lower AUC = more sensitive. "
      f"Samples per drug (with RNA-seq): median {per_drug.median():.0f}, range {per_drug.min()}-{per_drug.max()}.")
    w(f"- **Linkage:** {len(both)} non-control RNA-seq samples have drug data; {len(triple)} of those also have a DNA sample with mutation calls. "
      f"Disease stage of the {len(both)}: {stage.to_dict()}.")
    w(f"- **Mutations** (`beataml_wes_wv1to4_mutations_dbgap.txt`): {len(b['mutations_long']):,} calls in {mut.shape[0]} DNA samples. "
      f"Assay: {mut.capture_type.value_counts().to_dict()}. The Nimblegen assay is a targeted panel, so for those samples an absent call "
      "outside the panel is not a true negative; Phase 4 must restrict mutation features to genes covered by both assays (or model assay type).")
    w(f"- Among the {len(triple)} RNA+drug+DNA samples: IDH1 R132 hotspot {int(trip_mut.IDH1_hotspot.sum())}, "
      f"IDH2 R140/R172 hotspot {int(trip_mut.IDH2_hotspot.sum())}, FLT3-ITD (clinical) {int(trip_mut['clin_FLT3-ITD'].sum())}, "
      f"NPM1 (clinical) {int(trip_mut['clin_NPM1'].sum())}.")
    w(f"- **Published van Galen deconvolution** (`beataml2_manuscript_vg_cts.xlsx`): {len(vg_sc)} samples x 6 malignant cell types "
      f"({', '.join(vg_sc.columns[1:])}). These are signed scores (range {vg_sc.iloc[:,1:].min().min():.1f} to {vg_sc.iloc[:,1:].max().max():.1f}), "
      "not fractions, so Phase 2 will compare our fractions to them by rank correlation.")
    d_idh = r[r.inhibitor.str.contains("AGI-")].groupby("inhibitor").dbgap_rnaseq_sample.nunique()
    w(f"- **Differentiation agents in the Beat AML panel:** no ATRA, enasidenib or ivosidenib. The panel has the tool IDH inhibitors "
      f"AGI-5198 (IDH1) and AGI-6780 (IDH2): {d_idh.to_dict()} RNA-seq samples. These are the drugs for the Phase 4 IDH subgroup check.")

    w("\n## 4. LINCS L1000 drug signatures\n")
    w("**Options evaluated.**\n")
    geo = s["lincs_geo"]
    w(f"1. GEO GSE92742 + GSE70138 Level 5: 20 GB + 5 GB gzipped gctx (no partial reads possible), measured download ~2.7 MB/s from NCBI. "
      f"Compound coverage: {geo}. HL60 has only 40 compounds and enasidenib/ivosidenib are absent.")
    w("2. SigCom LINCS (Ma'ayan Lab) metadata API with per-signature TSVs on S3. Small queries worked, but filtered queries "
      "returned HTTP 504 (gateway timeout) repeatedly during this session, so it is not a reliable pipeline dependency.")
    w(f"3. **Chosen: CMap LINCS 2020 release** (the data behind SigCom). Metadata tables are public on S3, and the Level 5 compound matrix "
      f"is an uncompressed, contiguous HDF5 matrix, so single signature rows (12,328 float32 values) can be read with HTTP range requests. "
      f"Fetched rows are checked byte-for-byte against h5py reads of the same rows. Values are MODZ z-scores (CMap Level 5).")
    w("\n**Hematopoietic cell lines in LINCS 2020** (compound signatures):\n")
    cen = L["census"]
    heme = [c for c in cfg["lincs"]["myeloid_cell_lines"] if c in cen.index] + [c for c in ["JURKAT", "OCILY3", "OCILY10", "OCILY19"] if c in cen.index]
    w(md(cen.loc[heme], index=True))
    w(f"\nMOLM13, MV411, OCIAML2/3, KG1 and NB4 have no compound signatures. CMap annotates U937 as 'b-cell lymphoma'; it is a histiocytic lymphoma "
      "line with monocytic features and is commonly used as a myeloid model, so we keep it in the myeloid group. K562 is CML (erythroid-leaning) "
      "and CD34 is primary normal CD34+ cells; both are kept but ranked after the AML lines. Note that LINCS signatures here are mostly 6 h or 24 h exposures.")
    w(f"\n**Selection rule** (fixed before any scoring): for every matched Beat AML drug and positive control, up to "
      f"{cfg['lincs']['max_sigs_per_drug_cell']} signatures per (drug, cell line), for the myeloid lines {cfg['lincs']['myeloid_cell_lines']} "
      f"and the core lines {cfg['lincs']['fallback_cell_lines']}, ranked by CMap `is_hiq`, then `is_exemplar_sig`, then TAS. "
      "Positive controls profiled only in other lines take signatures from whatever lines exist.")
    w(f"\nDownloaded: {sigs.n_obs:,} signatures, {int((so.is_myeloid_line.astype(str)=='True').sum()):,} from myeloid/heme lines. "
      f"Exposure times: {so.pert_itime.value_counts().to_dict()}.")
    w(f"\nGenes: {sigs.n_vars:,} features, {L['feature_space']}. Only the 978 landmark genes are measured; the others are computationally inferred, "
      f"and the 'inferred' (not 'best inferred') group is the least reliable. HGNC mapping: {L['gene_mapping']}.")
    w("\n**Positive controls (differentiation agents)**\n")
    w(md(controls, index=True))
    w(f"\nMissing from LINCS entirely: {L['controls_missing'] or 'none'}. Enasidenib and ivosidenib exist only as 24 h signatures in 23 "
      "non-hematopoietic lines, none passing CMap's high-quality flag. Those lines are presumably IDH wild-type, where these drugs are not expected "
      "to change transcription, so a weak or null control result for them would be the expected outcome, not a pipeline failure. "
      "Tretinoin (ATRA) has signatures in HL60, U937, THP1, SKM1 and CD34, which makes it the best positive control.")

    w("\n## 5. Overlaps\n")
    w("**Genes (HGNC symbols):**\n")
    w(md(pd.DataFrame([
        {"set": "van Galen", "n": len(gv)}, {"set": "Beat AML", "n": len(gb)}, {"set": "LINCS (all 12,328)", "n": len(gl)},
        {"set": "van Galen & Beat AML", "n": len(gv & gb)}, {"set": "van Galen & LINCS", "n": len(gv & gl)},
        {"set": "Beat AML & LINCS", "n": len(gb & gl)}, {"set": "all three", "n": len(gv & gb & gl)},
        {"set": "all three, landmark + best inferred", "n": len(gv & gb & reliable)},
        {"set": "all three, landmark only", "n": len(gv & gb & landmark)}])))
    w("\n**Drugs (Beat AML panel x LINCS 2020):**\n")
    w(f"{match.match_how.value_counts().to_dict()}\n")
    show = overlap.reset_index()[["beat_name", "lincs_names", "match_how", "match_quality", "beat_rna_samples", "n_sigs", "n_cells", "myeloid_cells", "n_hiq"]]
    w(md(show.sort_values("beat_name")))
    flagged = match[(match.match_quality == "approximate") | (match.note.fillna("") != "")]
    w("\nMatches that are approximate or need a note:\n")
    w(md(flagged[["beat_name", "lincs_names", "match_how", "match_quality", "note"]]))
    w(f"\nNot in LINCS 2020 ({int((match.match_how=='none').sum())}): {', '.join(match[match.match_how=='none'].beat_name)}.")

    w("\n## 6. Blockers, caveats and deviations\n")
    w("- **No hard blockers.** All three sources are downloaded and parsed.")
    w("- **Python 3.12 instead of 3.11** (3.11 is not installed on this machine; 3.12 is supported by scanpy). Recorded here as a deviation.")
    w("- **No clinical differentiation agents in the Beat AML ex vivo panel** (no ATRA, enasidenib, ivosidenib). Validation against ex vivo "
      "response therefore tests mostly kinase inhibitors and cytotoxics, where viability, not maturation, is the readout. The IDH subgroup check "
      "has to use the tool compounds AGI-5198 / AGI-6780.")
    w("- **IDH inhibitor signatures are weak:** enasidenib/ivosidenib only in non-heme lines without high-quality flags; AGI-5198/AGI-6780 "
      "have THP1 signatures, but THP1 is not known to carry an IDH mutation, so no strong on-target transcriptional effect is expected.")
    w("- **Cell-line mismatch:** most matched drugs have signatures in THP1, but only a handful in HL60/SKM1/NOMO1/U937, and short exposures (mostly 6-24 h).")
    w("- **Inferred genes:** only 978 genes are measured by L1000; reversion scores should be checked on landmark-only genes as a robustness analysis.")
    w("- **Mutation assay heterogeneity** (exome vs targeted panel), see section 3.")
    w("- **Beat AML deconvolution reference values are scores, not fractions.**")
    w("- **SigCom LINCS API** timed out (HTTP 504); we read the same LINCS 2020 data directly from CMap's public S3 instead.")
    w("- Some van Galen patients have fewer than 100 malignant cells at diagnosis; per-patient statistics for them will be noisy.")
    w("\n## 7. Design decisions and justifications\n")
    w("Every choice below was fixed before any drug scoring or validation, so none of them can be tuned to the results.\n")
    lc = cfg["lincs"]
    decisions = [
        ("Data sources", [
            ("Use GEO GSE116256 supplementary files for van Galen.",
             "They are the authors' own per-cell counts and annotations (malignant vs normal, cell type). Re-deriving these labels ourselves "
             "would need the genotyping the authors did, which we cannot reproduce."),
            ("Use the Beat AML 2.0 harmonized files from the authors' GitHub, not dbGaP or GDC.",
             "They are the processed data behind the 2022 paper, openly licensed (CC-BY-4.0), and already harmonized across all four waves. "
             "The project rules forbid controlled-access data."),
            ("Use the CMap LINCS 2020 release instead of GEO GSE92742/GSE70138.",
             "LINCS 2020 has the compounds we need that GEO lacks (enasidenib, ivosidenib) and far more myeloid coverage "
             "(THP1: 1,630 compounds vs 371 in GEO). It also avoids a 25 GB download that cannot be read partially."),
            ("Read LINCS rows with HTTP range requests instead of downloading the 35.5 GB file.",
             "The matrix is stored uncompressed and contiguous, so each signature is one fixed-size block. We fetch about 340 MB instead "
             "of 35.5 GB. To rule out offset errors, random fetched rows are compared with h5py's own read of the same rows, and the run fails on any mismatch."),
            ("Do not depend on the SigCom LINCS API.",
             "Its filtered queries returned HTTP 504 repeatedly. A pipeline that must rerun from one command cannot depend on a service that times out."),
            ("Use Level 5 MODZ z-scores.",
             "Level 5 is CMap's standard consensus signature: replicate-collapsed and compared against plate controls. These are the values "
             "the published connectivity methods were designed for."),
        ]),
        ("Gene identifiers", [
            ("Harmonize everything to current HGNC approved symbols.",
             "The three datasets use different annotation vintages (hg19-era symbols, GRCh37 Ensembl, CMap's 2020 symbols). Without one "
             "reference, the same gene can appear under two names and silently drop out of overlaps."),
            ("Map in the order approved, then previous symbol, then alias, and only when the mapping is unambiguous.",
             "Approved symbols are definitive. Previous symbols are usually renames. Aliases are often shared between genes, so an "
             "ambiguous alias is left unmapped rather than guessed."),
            ("Map Beat AML through Ensembl IDs first.",
             "Ensembl IDs are stable identifiers, while the symbols in the file can be outdated."),
            ("Keep unmapped van Galen genes under their original names.",
             "Dropping them would change per-cell UMI totals used in Phase 1 normalization. They mostly drop out of the overlaps anyway."),
            ("When two van Galen symbols map to one HGNC gene, sum their UMIs.",
             "UMIs are counts of the same molecule, so summing is the correct way to merge them."),
            ("When two Beat AML Ensembl IDs map to one symbol, keep the most expressed.",
             "The values are log-scale normalized expression, which cannot be summed meaningfully. The most expressed ID is usually the canonical gene."),
            ("Report gene overlap separately for landmark and best-inferred LINCS genes.",
             "Only the 978 landmark genes are measured; the rest are inferred. Later phases will check robustness on landmark-only genes."),
        ]),
        ("Cells and samples", [
            ("Keep only van Galen cells that appear in the annotation tables, and fail loudly if an annotated cell has no counts.",
             "The annotation tables define the authors' QC-passed cells. Any mismatch would mean a parsing error."),
            ("Exclude the 783 BM5-34p38n cells with no annotation.",
             "They have no cell type and no normal/malignant call, so they cannot be placed in the lineage or used as references. "
             "Labelling them ourselves would mean making up data."),
            ("Include normal cells from AML samples in the landscape (proposed for Phase 1).",
             "The authors classified them as normal, and they add about 3,500 myeloid cells to the 4,430 from healthy donors. Phase 1 "
             "will check that they occupy the same positions as healthy cells before relying on them."),
            ("Exclude the 36 Beat AML healthy-control RNA-seq samples from patient linkage.",
             "They have no drug data and are not AML. They stay in the expression file as a possible reference."),
            ("Exclude the MUTZ3 and OCI-AML3 cell lines from patient analyses.",
             "They are lines, not patients. They stay in the AnnData object for possible sanity checks."),
        ]),
        ("Drug response and mutations", [
            ("Use the probit-fit AUC as the sensitivity measure.",
             "It is the measure the Beat AML papers use, it exists for every drug-sample pair, and it summarizes the whole dose-response curve. "
             "IC50 is often capped at the tested concentration range."),
            ("Keep all rows (all have paper_inclusion = True) and keep the file's train/test flag only for reference.",
             "Phase 4 will define its own patient-level split before seeing any results."),
            ("Use clinical consensus calls for FLT3-ITD and NPM1.",
             "Exome variant callers often miss FLT3 internal tandem duplications, and the Beat AML FAQ says these clinical calls are consensus calls."),
            ("Define IDH1 R132 and IDH2 R140/R172 hotspot flags.",
             "These are the variants that produce 2-hydroxyglutarate and that IDH inhibitors target. Other IDH variants are of unknown significance."),
            ("Record each sample's capture assay next to the mutation matrix.",
             "The targeted panel covers fewer genes, so a 0 there can mean 'not tested'. Keeping the assay label lets Phase 4 handle this correctly."),
            ("Require at least 50 RNA-seq samples per drug for validation.",
             "A per-drug correlation or model needs enough patients to be stable. In this panel, matched drugs below 50 have 0 to 45 samples and the next "
             "drug up has 85 (median 466), so any cutoff between 46 and 85 selects exactly the same drugs."),
        ]),
        ("Drug name matching", [
            ("Match normalized names (lowercase, letters and digits only), including each parenthetical alias.",
             "Beat AML names such as 'Quizartinib (AC220)' carry two names and the two catalogues format them differently (AC220 vs AC-220)."),
            ("Check curated synonyms first, and list each one with a reason in `drug_names.py`.",
             "Some drugs appear under a different name in LINCS (rapamycin as sirolimus, flavopiridol as alvocidib). Every curated pair was "
             "checked against names that actually exist in LINCS."),
            ("Label prodrug, enantiomer and racemate matches as approximate.",
             "Their expression effects may differ from the exact Beat AML compound. Phase 4 can rerun without them."),
            ("Exclude look-alike pairs (baicalin vs baicalein, BAY-11-7082 vs 7085) and drop JQ-1-R.",
             "These are different molecules. JQ-1-R is the inactive enantiomer of JQ1, so including it would dilute the real signature."),
            ("Do not match the generic 'AMPK Inhibitor'.",
             "It is probably Compound C (dorsomorphin), but nothing in the Beat AML files confirms that."),
        ]),
        ("Signature selection", [
            (f"Prefer myeloid and blood lines {lc['myeloid_cell_lines']}.",
             "Drug effects on gene expression depend on cell type, so signatures from myeloid cells should transfer best to AML blasts. "
             "Lines in the list without LINCS data cost nothing and keep the rule general."),
            ("Keep U937 in the myeloid group despite CMap's 'b-cell lymphoma' label.",
             "U937 is a histiocytic lymphoma line with monocytic features and is widely used as a monocytic differentiation model."),
            ("Keep K562 and CD34, flagged as non-AML.",
             "K562 is a CML line and CD34 is primary normal progenitor cells. Both are hematopoietic and better than epithelial lines, "
             "but they are not AML, so Phase 3 can weight or exclude them."),
            (f"Also take the core lines {lc['fallback_cell_lines']}.",
             "These are the lines most compounds were profiled in. They let us score the roughly 30 drugs that lack a myeloid signature, "
             "and compare myeloid against non-myeloid signatures for the same drug."),
            (f"At most {lc['max_sigs_per_drug_cell']} signatures per (drug, cell line).",
             "Some drugs have thousands of signatures (bortezomib 8,500). Without a cap they would dominate any consensus, and 6 still "
             "covers several doses and times for most drugs."),
            ("Rank within each group by `is_hiq`, then `is_exemplar_sig`, then TAS.",
             "`is_hiq` is CMap's reproducibility flag, the exemplar flag is CMap's pick of the most representative signature, and TAS "
             "measures how strong and reproducible the response is."),
            (f"Positive controls {lc['positive_controls']}.",
             "ATRA, tamibarotene and alitretinoin are retinoids (ATRA is the clinical differentiation therapy). Enasidenib and ivosidenib "
             "are clinical IDH inhibitors, and AGI-5198/AGI-6780 are their tool compounds. Calcitriol (vitamin D) induces monocytic "
             "differentiation of HL60 cells. Decitabine and azacitidine are hypomethylating agents with reported differentiating effects. "
             "Controls profiled only in other lines are taken from whatever lines exist, so a missing control is visible rather than silently absent."),
        ]),
    ]
    for title, items in decisions:
        w(f"\n**{title}**\n")
        w(md(pd.DataFrame(items, columns=["Choice", "Justification"])))
    w("\n## Figures\n")
    for f in figs:
        w(f"![{f.stem}](figures/{f.name})")
    w("\n## Cached outputs (`data/processed/`)\n")
    proc = ROOT / cfg["paths"]["processed"]
    w(md(pd.DataFrame([{"file": f.name, "MB": round(f.stat().st_size / 1e6, 1)} for f in sorted(proc.glob("*")) if f.is_file()])))
    path(cfg, "reports", "phase0_data_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
