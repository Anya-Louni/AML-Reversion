"""Parse Beat AML 2.0 open processed data (Bottomly et al. 2022, CC-BY-4.0)."""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from pipeline.common import get_logger
from pipeline.phase0_data.genes import HGNCMapper

log = get_logger("beataml")

IDH_HOTSPOTS = {"IDH1": ("p.R132",), "IDH2": ("p.R140", "p.R172")}


def load_expression(raw: Path, hgnc: HGNCMapper) -> tuple[ad.AnnData, dict]:
    df = pd.read_csv(raw / "beataml_waves1to4_norm_exp_dbgap.txt", sep="\t")
    meta_cols = ["stable_id", "display_label", "description", "biotype"]
    samples = [c for c in df.columns if c not in meta_cols]
    df = df.copy()  # de-fragment the wide frame before adding columns
    sym = hgnc.map_ensembl(df.stable_id).values
    how = np.where(pd.notna(sym), "ensembl", None)
    # Fall back to the provided symbol when the Ensembl ID is not in HGNC.
    for i in np.where(pd.isna(sym))[0]:
        s, h = hgnc.map_symbol(df.display_label.iat[i])
        sym[i], how[i] = (s, f"symbol_{h}") if s else (None, "unmapped")
    df["hgnc"], df["hgnc_how"] = sym, how
    n_unmapped = int(pd.isna(sym).sum())
    df = df.dropna(subset=["hgnc"])
    # Several Ensembl IDs can map to one symbol; keep the one with the highest mean expression.
    df["mean"] = df[samples].mean(axis=1)
    n_before = len(df)
    df = df.sort_values("mean", ascending=False).drop_duplicates("hgnc")
    X = df[samples].T.values.astype(np.float32)
    var = df.set_index("hgnc")[["stable_id", "display_label", "biotype", "hgnc_how"]]
    a = ad.AnnData(X=X, obs=pd.DataFrame(index=pd.Index(samples, name="dbgap_rnaseq_sample")), var=var)
    stats = {"n_genes_input": n_before + n_unmapped, "n_unmapped_dropped": n_unmapped,
             "n_duplicate_symbols_dropped": n_before - len(df), "n_genes": a.n_vars, "n_samples": a.n_obs,
             "value_range": (float(X.min()), float(X.max())), "hgnc_how": var.hgnc_how.value_counts().to_dict()}
    return a, stats


def load_drug_response(raw: Path) -> pd.DataFrame:
    p = pd.read_csv(raw / "beataml_probit_curve_fits_v4_dbgap.txt", sep="\t")
    keep = ["dbgap_subject_id", "dbgap_dnaseq_sample", "dbgap_rnaseq_sample", "inhibitor", "status",
            "min_conc", "max_conc", "ic50", "auc", "converged", "curve_type", "all_gt_50", "all_lt_50"]
    return p[keep]


def load_mutations(raw: Path, clinical: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (long mutation calls, per-DNA-sample binary gene matrix with assay type)."""
    m = pd.read_csv(raw / "beataml_wes_wv1to4_mutations_dbgap.txt", sep="\t")
    mat = (m.assign(v=1).pivot_table(index="dbgap_sample_id", columns="symbol", values="v",
                                     aggfunc="max", fill_value=0).astype(np.int8))
    for gene, prefixes in IDH_HOTSPOTS.items():
        hits = m[(m.symbol == gene) & m.hgvsp_short.fillna("").str.startswith(prefixes)].dbgap_sample_id
        mat[f"{gene}_hotspot"] = mat.index.isin(hits).astype(np.int8)
    cap = m.groupby("dbgap_sample_id").capture_type.agg(lambda s: "|".join(sorted(set(s))))
    mat.insert(0, "capture_type", cap.reindex(mat.index))
    # Clinical consensus calls (FLT3-ITD is not reliably called from exome data).
    c = clinical.dropna(subset=["dbgap_dnaseq_sample"]).drop_duplicates("dbgap_dnaseq_sample")
    c = c.set_index("dbgap_dnaseq_sample")
    for col in ["FLT3-ITD", "NPM1"]:
        mat[f"clin_{col}"] = c[col].reindex(mat.index).map({"positive": 1, "negative": 0})
    return m, mat


def load_all(raw: Path, hgnc: HGNCMapper) -> dict:
    clinical = pd.read_excel(raw / "beataml_wv1to4_clinical.xlsx", sheet_name="summary")
    mapping = pd.read_excel(raw / "beataml_waves1to4_sample_mapping.xlsx")
    expr, expr_stats = load_expression(raw, hgnc)
    resp = load_drug_response(raw)
    muts_long, muts = load_mutations(raw, clinical)
    vg = pd.read_excel(raw / "beataml2_manuscript_vg_cts.xlsx")
    families = pd.read_excel(raw / "beataml_drug_families.xlsx", sheet_name="drug_family")

    ctrl = mapping.dropna(subset=["dbgap_rnaseq_sample"]).set_index("dbgap_rnaseq_sample")["rna_control"]
    expr.obs["rna_control"] = ctrl.reindex(expr.obs_names).astype(str).replace("nan", "")
    clin_r = clinical.dropna(subset=["dbgap_rnaseq_sample"]).drop_duplicates("dbgap_rnaseq_sample")
    clin_r = clin_r.set_index("dbgap_rnaseq_sample")
    for col in ["dbgap_subject_id", "dbgap_dnaseq_sample", "diseaseStageAtSpecimenCollection",
                "specimenType", "ELN2017", "FLT3-ITD", "NPM1", "cohort"]:
        expr.obs[col] = clin_r[col].reindex(expr.obs_names).astype(str).replace("nan", "")
    return {"clinical": clinical, "mapping": mapping, "expr": expr, "expr_stats": expr_stats,
            "response": resp, "mutations_long": muts_long, "mutations": muts, "vg_scores": vg,
            "drug_families": families}
