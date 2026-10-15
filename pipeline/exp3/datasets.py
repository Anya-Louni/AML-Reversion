"""Loaders for the Experiment 3 treated-cell datasets (see reports/exp3_step0_data.md).

Every loader returns raw counts as a genes x samples DataFrame indexed by HGNC symbol, plus a
sample table giving the arm of each column. Nothing is normalized here.
"""
from __future__ import annotations

import gzip
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline.common import download, get_logger
from pipeline.phase0_data.genes import HGNCMapper

log = get_logger("exp3.datasets")

GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series"
SOURCES = {
    "GSE227542": f"{GEO}/GSE227nnn/GSE227542/suppl/GSE227542_RAW.tar",
    "GSE276457": f"{GEO}/GSE276nnn/GSE276457/suppl/GSE276457_RNASeq_raw_count_data.txt.gz",
    "GSE228321": f"{GEO}/GSE228nnn/GSE228321/suppl/GSE228321_Bulk_RNASeq_RawCount.txt.gz",
}


def fetch(raw: Path) -> None:
    for name, url in SOURCES.items():
        download(url, raw / url.rsplit("/", 1)[1])


def _collapse(df: pd.DataFrame) -> pd.DataFrame:
    """Sum counts of rows mapping to the same symbol."""
    return df.groupby(level=0).sum()


def load_gse227542(raw: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """HL-60 ATRA time course: 0 h, 4 h, 2 d, 4 d, three replicates each."""
    tar = raw / "GSE227542_RAW.tar"
    cols, meta = {}, []
    with tarfile.open(tar) as t:
        for m in t.getmembers():
            g = re.match(r"(GSM\d+)_(Ctrl|RA4h|RA2d|RA4d)_(\d)\.count\.txt\.gz$", Path(m.name).name)
            if not g:
                continue
            with gzip.open(t.extractfile(m), "rt") as f:
                s = pd.read_csv(f, sep="\t", header=None, index_col=0).iloc[:, 0]
            name = f"{g.group(2)}_{g.group(3)}"
            cols[name] = s
            meta.append({"sample": name, "gsm": g.group(1),
                         "arm": "control" if g.group(2) == "Ctrl" else g.group(2),
                         "time": {"Ctrl": "0h", "RA4h": "4h", "RA2d": "2d", "RA4d": "4d"}[g.group(2)]})
    X = pd.DataFrame(cols)
    X.index = X.index.astype(str)
    X = X[~X.index.str.startswith(("__", "N_"))]
    return _collapse(X), pd.DataFrame(meta).sort_values("sample").reset_index(drop=True)


def load_gse276457(raw: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """HL-60 48 h: untreated, ATRA, palbociclib (CDK4i), combination; three replicates each."""
    df = pd.read_csv(raw / "GSE276457_RNASeq_raw_count_data.txt.gz", sep="\t")
    X = df.set_index("Geneid").drop(columns=[c for c in ("Chr", "Start", "End", "Strand", "Length") if c in df])
    X.index = X.index.astype(str)
    arm = {"cont": "control", "ATRA": "ATRA", "CDK4i": "palbociclib", "combi": "combination"}
    meta = [{"sample": c, "arm": next(a for k, a in arm.items() if c.startswith(k))} for c in X.columns]
    return _collapse(X), pd.DataFrame(meta)


def load_gse228321(raw: Path, hgnc: HGNCMapper) -> tuple[pd.DataFrame, pd.DataFrame]:
    """MOLM13, MV4-11 and OCI-AML3 at 16 h: control, VTP50469, OTX015, VTP+OTX (MOLM13 only)."""
    X = pd.read_csv(raw / "GSE228321_Bulk_RNASeq_RawCount.txt.gz", sep="\t", index_col=0)
    sym = hgnc.map_ensembl([i.split(".")[0] for i in X.index])
    X.index = pd.Index(sym.values, name="gene")
    X = _collapse(X[X.index.notna()])
    meta = []
    for c in X.columns:
        line = "MOLM13" if "MOLM" in c else "MV4-11" if "MV411" in c else "OCI-AML3"
        if "VTP_OTX" in c:
            a = "combination"
        elif "OTX" in c:
            a = "OTX015"
        elif "VTP" in c:
            a = "VTP50469"
        else:
            a = "control"
        meta.append({"sample": c, "cell_line": line, "arm": a})
    return X, pd.DataFrame(meta)


def log2cpm(X: pd.DataFrame) -> pd.DataFrame:
    """Counts -> log2(CPM + 1), per sample."""
    return np.log2(X / X.sum(axis=0) * 1e6 + 1)


def real_change(X: pd.DataFrame, meta: pd.DataFrame, treated: str, control: str = "control",
                subset: dict | None = None) -> pd.Series:
    """Mean log2 fold change of a treated arm against its control, over genes detected in the contrast."""
    m = meta
    if subset:
        for k, v in subset.items():
            m = m[m[k] == v]
    tr = m.loc[m.arm == treated, "sample"]
    ct = m.loc[m.arm == control, "sample"]
    if len(tr) == 0 or len(ct) == 0:
        raise ValueError(f"no samples for {treated} vs {control} with {subset}")
    cols = list(tr) + list(ct)
    keep = X[cols].sum(axis=1) > 0
    L = log2cpm(X.loc[keep, cols])
    return L[tr].mean(axis=1) - L[ct].mean(axis=1)
