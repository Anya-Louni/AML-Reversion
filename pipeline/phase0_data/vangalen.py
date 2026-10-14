"""Build an AnnData object from van Galen et al. 2019 (GEO GSE116256).

The GEO archive holds, per sample, a dense gene x cell UMI table (`*.dem.txt.gz`)
and a per-cell annotation table (`*.anno.txt.gz`) with columns such as
PredictionRefined (malignant/normal) and CellType. Nanopore genotyping tables
(`*.nanopore.txt.gz`) are not used.
"""
from __future__ import annotations

import re
import tarfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from pipeline.common import get_logger
from pipeline.phase0_data.genes import HGNCMapper, summarize_mapping

log = get_logger("vangalen")

FNAME = re.compile(r"^(GSM\d+)_(.+)\.(dem|anno)\.txt\.gz$")


def sample_info(sample: str) -> dict:
    """Parse a sample name such as 'AML328-D29', 'BM5-34p' or 'OCI-AML3'."""
    m = re.match(r"^(AML\w+?)-D(\d+)$", sample)
    if m:
        return {"sample": sample, "patient": m.group(1), "day": int(m.group(2)), "sample_type": "AML"}
    if sample.startswith("BM"):
        return {"sample": sample, "patient": sample.split("-")[0], "day": np.nan, "sample_type": "healthy_BM"}
    return {"sample": sample, "patient": sample, "day": np.nan, "sample_type": "cell_line"}


def extract(tar_path: Path, out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path) as t:
        members = [m for m in t.getmembers() if FNAME.match(Path(m.name).name)]
        for m in members:
            if not (out_dir / Path(m.name).name).exists():
                t.extract(m, out_dir, filter="data")
    return sorted(out_dir.glob("*.txt.gz"))


def read_sample(dem: Path, anno: Path) -> ad.AnnData:
    counts = pd.read_csv(dem, sep="\t", index_col=0)
    meta = pd.read_csv(anno, sep="\t", index_col=0)
    missing = meta.index.difference(counts.columns)
    if len(missing):
        raise ValueError(f"{anno.name}: {len(missing)} annotated cells missing from counts")
    counts = counts[meta.index]  # keep annotated cells only, in annotation order
    X = sp.csr_matrix(counts.values.T.astype(np.float32))
    a = ad.AnnData(X=X, obs=meta, var=pd.DataFrame(index=counts.index.astype(str)))
    a.uns["n_cells_in_counts"] = counts.shape[1]
    return a


def build(raw_dir: Path, hgnc: HGNCMapper) -> tuple[ad.AnnData, dict]:
    files = {}
    for f in raw_dir.glob("*.txt.gz"):
        m = FNAME.match(f.name)
        if m:
            files.setdefault(m.group(2), {})[m.group(3)] = f
    adatas, per_sample = [], []
    for sample in sorted(files):
        pair = files[sample]
        if set(pair) != {"dem", "anno"}:
            log.warning("skipping %s: has %s", sample, sorted(pair))
            continue
        a = read_sample(pair["dem"], pair["anno"])
        info = sample_info(sample)
        for k, v in info.items():
            a.obs[k] = v
        per_sample.append({**info, "n_cells_counts": a.uns.pop("n_cells_in_counts"),
                           "n_cells_annotated": a.n_obs, "n_genes": a.n_vars})
        adatas.append(a)
        log.info("%s: %d cells", sample, a.n_obs)
    genes0 = adatas[0].var_names
    same = all(a.var_names.equals(genes0) for a in adatas)
    adata = ad.concat(adatas, join="outer" if not same else "inner", merge="same")
    adata.obs_names_make_unique()

    # Harmonize gene symbols; collapse symbols that map to the same HGNC gene by summing UMIs.
    m = hgnc.map_symbols(adata.var_names)
    target = m["hgnc"].fillna(m["original"]).values
    stats = summarize_mapping(m)
    if len(set(target)) < len(target):
        codes, uniq = pd.factorize(target)
        agg = sp.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)),
                            shape=(len(codes), len(uniq)), dtype=np.float32)
        X = (adata.X @ agg).tocsr()
        var = pd.DataFrame(index=pd.Index(uniq, name="gene"))
        var["original_symbols"] = pd.Series(adata.var_names).groupby(codes).agg("|".join).values
        adata = ad.AnnData(X=X, obs=adata.obs, var=var)
    else:
        adata.var_names = target
        adata.var["original_symbols"] = m["original"].values
    adata.var["hgnc_approved"] = adata.var_names.isin(hgnc.approved)

    for c in adata.obs.columns:
        if adata.obs[c].dtype == object:
            adata.obs[c] = adata.obs[c].astype(str).replace("nan", np.nan)
    adata.obs["n_counts"] = np.asarray(adata.X.sum(1)).ravel()
    adata.obs["n_genes_detected"] = np.asarray((adata.X > 0).sum(1)).ravel()
    summary = {"gene_mapping": stats, "samples": pd.DataFrame(per_sample),
               "identical_gene_lists_across_samples": same}
    return adata, summary
