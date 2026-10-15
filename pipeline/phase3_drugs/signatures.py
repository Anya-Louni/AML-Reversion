"""Drug signatures: wide myeloid fetch and per-drug consensus."""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from pipeline.common import get_logger
from pipeline.phase0_data import lincs

log = get_logger("phase3.signatures")


def fetch_wide(cfg: dict, raw: Path, proc: Path) -> ad.AnnData:
    """All compound signatures in AML/monocytic lines (capped), merged with the Phase 0 set."""
    out = proc / "lincs_signatures_all.h5ad"
    if out.exists():
        return ad.read_h5ad(out)
    p3 = cfg["phase3"]
    base = ad.read_h5ad(proc / "lincs_signatures.h5ad")
    sig = lincs.load_siginfo(raw)
    names = sig[sig.cell_iname.isin(p3["wide_cell_lines"])].cmap_name.dropna().unique()
    match = pd.read_csv(proc / "drug_matching.csv").dropna(subset=["lincs_names"])
    to_beat = {n: r.beat_name for r in match.itertuples() for n in r.lincs_names.split("|")}
    drugs = {n: to_beat.get(n, n) for n in names}
    sel = lincs.select_signatures(sig, drugs, p3["wide_cell_lines"], p3["wide_max_sigs_per_cell"])
    sel = sel[~sel.sig_id.isin(base.obs_names)].reset_index(drop=True)
    log.info("wide screen: %d compounds, %d new signatures to fetch", len(names), len(sel))
    remote = lincs.RemoteLevel5(cache_dir=proc)
    idx = remote.sig_ids.get_indexer(sel.sig_id)
    sel, idx = sel[idx >= 0].reset_index(drop=True), idx[idx >= 0]
    part_dir = proc / "lincs_wide_parts"; part_dir.mkdir(exist_ok=True)
    X = remote.fetch_rows(idx, cache=part_dir / "rows.npy")
    # Offset check without reopening the remote HDF5: refetch rows already verified against h5py in Phase 0.
    chk = base.obs_names[np.random.default_rng(1).choice(base.n_obs, size=3, replace=False)]
    got = remote.fetch_rows(remote.sig_ids.get_indexer(chk))
    if not np.allclose(got, np.asarray(base[chk].X), equal_nan=True):
        raise RuntimeError("range-request rows do not match Phase 0 verified rows")
    new = ad.AnnData(X=X, obs=sel.set_index("sig_id").astype({c: str for c in sel.columns if sel[c].dtype == object and c != "sig_id"}),
                     var=base.var.copy())
    new.obs["is_myeloid_line"] = True
    new.obs["is_positive_control"] = False
    base.obs["source_fetch"] = "phase0"; new.obs["source_fetch"] = "wide"
    for c in base.obs.columns:
        base.obs[c] = base.obs[c].astype(str)
    for c in new.obs.columns:
        new.obs[c] = new.obs[c].astype(str)
    allsig = ad.concat([base, new], merge="same")
    allsig.var = base.var.copy()
    allsig.write_h5ad(out, compression="gzip")
    return allsig


def consensus(sigs: ad.AnnData, p3: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Per (drug, cell) mean signature, then one drug-level signature from the best available line tier.

    Returns (drug x gene consensus, drug metadata, (drug, cell) x gene table).
    """
    obs = sigs.obs[["drug", "cell_iname", "tas", "is_hiq"]].copy()
    obs["tas"] = pd.to_numeric(obs.tas, errors="coerce")
    X = pd.DataFrame(np.asarray(sigs.X), index=sigs.obs_names, columns=sigs.var_names)
    key = obs.drug.astype(str) + "||" + obs.cell_iname.astype(str)
    dc = X.groupby(key.values).mean()
    dc.index = pd.MultiIndex.from_tuples([tuple(i.split("||")) for i in dc.index], names=["drug", "cell"])
    n_sig = key.value_counts()
    tiers = [("AML/monocytic", p3["primary_lines"]), ("CD34/K562", p3["secondary_lines"]), ("non-myeloid", p3["fallback_lines"])]
    rows, meta = {}, []
    for d, g in dc.groupby(level="drug"):
        cells = g.index.get_level_values("cell")
        for tier, lines in tiers:
            use = [c for c in cells if c in lines]
            if use:
                rows[d] = g.loc[[(d, c) for c in use]].mean()
                meta.append({"drug": d, "tier": tier, "cells": ",".join(sorted(use)),
                             "n_sigs": int(sum(n_sig.get(f"{d}||{c}", 0) for c in use)),
                             "mean_tas": float(obs[(obs.drug == d) & obs.cell_iname.isin(use)].tas.mean()),
                             "all_cells": ",".join(sorted(set(cells)))})
                break
    return pd.DataFrame(rows).T, pd.DataFrame(meta).set_index("drug"), dc
