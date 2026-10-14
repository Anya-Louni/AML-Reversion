"""Phase 0: data acquisition, parsing, harmonization and audit."""
from __future__ import annotations

import json

import anndata as ad
import pandas as pd

from pipeline.common import get_logger, path
from pipeline.phase0_data import beataml, lincs, vangalen
from pipeline.phase0_data.download import fetch_all
from pipeline.phase0_data.genes import HGNCMapper
from pipeline.phase0_data.report import write_report

log = get_logger("phase0")


def run(cfg: dict) -> None:
    fetch_all(cfg)
    hgnc = HGNCMapper.from_file(path(cfg, "raw", "hgnc", "hgnc_complete_set.txt"))
    proc = path(cfg, "processed", "x").parent
    summ: dict = {}

    vg_dir = path(cfg, "raw", "vangalen", "x").parent
    vangalen.extract(vg_dir / "GSE116256_RAW.tar", vg_dir)
    vg, summ["vangalen"] = vangalen.build(vg_dir, hgnc)
    vg.write_h5ad(proc / "vangalen_raw.h5ad", compression="gzip")

    b = beataml.load_all(path(cfg, "raw", "beataml", "x").parent, hgnc)
    b["expr"].write_h5ad(proc / "beataml_expr.h5ad", compression="gzip")
    b["response"].to_parquet(proc / "beataml_drug_response.parquet")
    b["mutations"].to_parquet(proc / "beataml_mutations_matrix.parquet")
    b["mutations_long"].to_parquet(proc / "beataml_mutations_long.parquet")
    b["clinical"].astype(str).to_parquet(proc / "beataml_clinical.parquet")
    b["vg_scores"].to_parquet(proc / "beataml_vg_celltype_scores.parquet")
    summ["beataml"] = b

    sig_path = proc / "lincs_signatures.h5ad"
    lraw = path(cfg, "raw", "lincs", "lincs2020", "x").parent
    beat_drugs = b["response"].inhibitor.unique()
    if sig_path.exists():
        log.info("using cached %s", sig_path.name)
        sigs = ad.read_h5ad(sig_path)
        summ["lincs"] = json.loads((proc / "lincs_summary.json").read_text())
        summ["lincs"]["match"] = pd.read_csv(proc / "drug_matching.csv")
        summ["lincs"]["census"] = pd.read_csv(proc / "lincs_cell_census.csv", index_col=0)
    else:
        sigs, ls = lincs.build(lraw, proc, beat_drugs, cfg, hgnc)
        sigs.write_h5ad(sig_path, compression="gzip")
        ls["match"].to_csv(proc / "drug_matching.csv", index=False)
        ls["census"].to_csv(proc / "lincs_cell_census.csv")
        (proc / "lincs_summary.json").write_text(json.dumps(
            {k: v for k, v in ls.items() if k not in ("match", "census")}, indent=2, default=str))
        summ["lincs"] = ls
    summ["lincs_sigs"] = sigs

    geo = {}
    for url in cfg["sources"]["lincs"]["geo_meta"]:
        s = pd.read_csv(path(cfg, "raw", "lincs", url.rsplit("/", 1)[1]), sep="\t", low_memory=False,
                        usecols=["pert_iname", "pert_type", "cell_id"])
        s = s[s.pert_type == "trt_cp"]
        geo[url.rsplit("/", 1)[1].split("_")[0]] = {"n_sigs": len(s), "n_compounds": s.pert_iname.nunique(),
                                                     "n_cells": s.cell_id.nunique()}
    summ["lincs_geo"] = geo

    write_report(cfg, summ, vg)
    log.info("phase 0 done")
