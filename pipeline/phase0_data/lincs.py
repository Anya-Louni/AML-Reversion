"""LINCS L1000 (CMap LINCS 2020 release) signatures for the drugs we need.

The Level 5 compound matrix (level5_beta_trt_cp_n720216x12328.gctx, 35.5 GB) is
an uncompressed HDF5 file on public S3 with a contiguous (signature x gene)
float32 matrix. We read only the selected signature rows with HTTP range
requests instead of downloading the whole file.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import anndata as ad
import fsspec
import h5py
import numpy as np
import pandas as pd
import requests

from pipeline.common import get_logger
from pipeline.phase0_data.drug_names import match_drugs
from pipeline.phase0_data.genes import HGNCMapper, summarize_mapping

log = get_logger("lincs")

LEVEL5_URL = ("https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level5/"
              "level5_beta_trt_cp_n720216x12328.gctx")
SIG_COLS = ["sig_id", "pert_id", "cmap_name", "cell_iname", "pert_idose", "pert_dose", "pert_dose_unit",
            "pert_itime", "pert_time", "nsample", "tas", "cc_q75", "is_hiq", "qc_pass", "is_exemplar_sig",
            "pct_self_rank_q25", "project_code", "pert_type"]


def load_siginfo(raw: Path) -> pd.DataFrame:
    s = pd.read_csv(raw / "siginfo_beta.txt", sep="\t", usecols=SIG_COLS, low_memory=False)
    return s[s.pert_type == "trt_cp"].drop(columns="pert_type")


def compound_aliases(raw: Path) -> dict[str, str]:
    c = pd.read_csv(raw / "compoundinfo_beta.txt", sep="\t")
    return {a: r.cmap_name for r in c.dropna(subset=["compound_aliases"]).itertuples()
            for a in str(r.compound_aliases).split("|")}


def select_signatures(sig: pd.DataFrame, drugs: dict[str, str], cells: list[str], cap: int) -> pd.DataFrame:
    """Pick up to `cap` signatures per (LINCS drug, cell line) for the given cells.

    `drugs` maps LINCS cmap_name -> label used downstream. Ranking inside each group:
    high-quality flag, exemplar flag, then transcriptional activity score (TAS).
    """
    s = sig[sig.cmap_name.isin(drugs) & sig.cell_iname.isin(cells)].copy()
    s["drug"] = s.cmap_name.map(drugs)
    s = s.sort_values(["is_hiq", "is_exemplar_sig", "tas"], ascending=False)
    return s.groupby(["cmap_name", "cell_iname"], sort=False).head(cap).reset_index(drop=True)


class RemoteLevel5:
    """Random access to rows of the remote Level 5 matrix."""

    def __init__(self, url: str = LEVEL5_URL, cache_dir: Path | None = None):
        self.url = url
        f = fsspec.open(url, "rb", block_size=4 * 2**20).open()
        self.h5 = h5py.File(f, "r")
        dset = self.h5["0/DATA/0/matrix"]
        if dset.chunks is not None or dset.compression is not None:
            raise RuntimeError("expected a contiguous, uncompressed matrix")
        self.dtype = dset.dtype
        self.shape = dset.shape
        self.offset = dset.id.get_offset()
        self.row_bytes = self.shape[1] * self.dtype.itemsize
        cache = cache_dir / "lincs2020_level5_ids.npz" if cache_dir else None
        if cache is not None and cache.exists():
            z = np.load(cache, allow_pickle=True)
            col, row = z["col"], z["row"]
        else:
            col = self.h5["0/META/COL/id"][:].astype(str)
            row = self.h5["0/META/ROW/id"][:].astype(str)
            if cache is not None:
                np.savez(cache, col=col, row=row)
        self.sig_ids = pd.Index(col)
        self.gene_ids = np.array(row)

    def read_h5py(self, i: int) -> np.ndarray:
        return self.h5["0/DATA/0/matrix"][i, :]

    def _fetch(self, i: int, session: requests.Session) -> np.ndarray:
        start = self.offset + i * self.row_bytes
        hdr = {"Range": f"bytes={start}-{start + self.row_bytes - 1}"}
        for attempt in range(5):
            try:
                r = session.get(self.url, headers=hdr, timeout=60)
                r.raise_for_status()
                if len(r.content) != self.row_bytes:
                    raise IOError("short read")
                return np.frombuffer(r.content, dtype=self.dtype).astype(np.float32)
            except Exception as e:
                log.warning("row %d attempt %d: %s", i, attempt + 1, e)
        raise RuntimeError(f"failed to fetch row {i}")

    def fetch_rows(self, idx: np.ndarray, workers: int = 16) -> np.ndarray:
        out = np.empty((len(idx), self.shape[1]), dtype=np.float32)
        sessions = [requests.Session() for _ in range(workers)]

        def job(k):
            out[k] = self._fetch(int(idx[k]), sessions[k % workers])
            if k and k % 500 == 0:
                log.info("fetched %d / %d rows", k, len(idx))

        with ThreadPoolExecutor(workers) as ex:
            list(ex.map(job, range(len(idx))))
        return out


def build(raw: Path, processed: Path, beat_drugs, cfg: dict, hgnc: HGNCMapper) -> tuple[ad.AnnData, dict]:
    lc = cfg["lincs"]
    sig = load_siginfo(raw)
    match = match_drugs(beat_drugs, sig.cmap_name.dropna().unique(), compound_aliases(raw))

    drugs: dict[str, str] = {}
    for r in match.dropna(subset=["lincs_names"]).itertuples():
        for n in r.lincs_names.split("|"):
            drugs[n] = r.beat_name
    controls_present = [c for c in lc["positive_controls"] if c in set(sig.cmap_name)]
    for c in controls_present:
        drugs.setdefault(c, c)

    cells = list(dict.fromkeys(lc["myeloid_cell_lines"] + lc["fallback_cell_lines"]))
    sel = select_signatures(sig, drugs, cells, lc["max_sigs_per_drug_cell"])
    # Positive controls profiled only outside our preferred lines: take what exists.
    extra = [c for c in controls_present if c not in set(sel.cmap_name)]
    if extra:
        any_cells = sig[sig.cmap_name.isin(extra)].cell_iname.unique().tolist()
        sel = pd.concat([sel, select_signatures(sig, {c: c for c in extra}, any_cells,
                                                lc["max_sigs_per_drug_cell"])], ignore_index=True)
    sel["is_myeloid_line"] = sel.cell_iname.isin(lc["myeloid_cell_lines"])
    sel["is_positive_control"] = sel.drug.isin(lc["positive_controls"]) | sel.cmap_name.isin(lc["positive_controls"])
    log.info("selected %d signatures for %d drugs", len(sel), sel.drug.nunique())

    remote = RemoteLevel5(cache_dir=processed)
    missing = ~sel.sig_id.isin(remote.sig_ids)
    if missing.any():
        log.warning("%d selected sig_ids not in Level 5 matrix; dropped", missing.sum())
        sel = sel[~missing].reset_index(drop=True)
    idx = remote.sig_ids.get_indexer(sel.sig_id)
    X = remote.fetch_rows(idx)
    # Consistency check: range-request bytes must equal h5py's own read for a few rows.
    for k in np.random.default_rng(0).choice(len(idx), size=min(3, len(idx)), replace=False):
        if not np.allclose(X[k], remote.read_h5py(int(idx[k])), equal_nan=True):
            raise RuntimeError("range-request row does not match h5py read")

    gi = pd.read_csv(raw / "geneinfo_beta.txt", sep="\t", dtype={"gene_id": str})
    gi = gi.set_index("gene_id").loc[remote.gene_ids]
    gm = hgnc.map_symbols(gi.gene_symbol)
    var = pd.DataFrame({"entrez_id": gi.index, "lincs_symbol": gi.gene_symbol.values,
                        "feature_space": gi.feature_space.values,
                        "hgnc_how": gm["how"].values},
                       index=pd.Index(gm["hgnc"].fillna(gm["original"]).values, name="gene"))
    obs = sel.set_index("sig_id")
    for c in obs.columns:
        if obs[c].dtype == object:
            obs[c] = obs[c].astype(str)
    a = ad.AnnData(X=X, obs=obs, var=var)
    a.var_names_make_unique()

    census = sig.groupby("cell_iname").agg(n_sigs=("sig_id", "size"), n_compounds=("cmap_name", "nunique"),
                                           n_hiq=("is_hiq", "sum"))
    summary = {"match": match, "census": census.sort_values("n_compounds", ascending=False),
               "n_trt_cp_sigs": len(sig), "n_compounds": sig.cmap_name.nunique(), "n_cells": sig.cell_iname.nunique(),
               "controls_present": controls_present,
               "controls_missing": [c for c in lc["positive_controls"] if c not in controls_present],
               "gene_mapping": summarize_mapping(gm), "feature_space": var.feature_space.value_counts().to_dict()}
    return a, summary
