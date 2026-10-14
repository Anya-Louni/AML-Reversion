"""NNLS deconvolution of bulk RNA-seq into van Galen cell types.

Reference profiles are mean CPM of van Galen cell groups: the six malignant
like-types plus non-myeloid normal groups (T/NK, B, erythroid, pDC). Normal
myeloid cells are not separate references, because they are nearly collinear
with the malignant like-types; in bulk they are absorbed into the matching
like-type. This mirrors the Bottomly et al. 2022 use of van Galen types, but with
NNLS fractions instead of signature scores.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.optimize import nnls


def label_cells(obs: pd.DataFrame, like_types, normal_groups: dict) -> pd.Series:
    lab = pd.Series(pd.NA, index=obs.index, dtype=object)
    mal = (obs.PredictionRefined == "malignant") & obs.CellType.isin(like_types) & (obs.sample_type == "AML")
    lab[mal.values] = obs.CellType[mal].astype(str).values
    norm = obs.PredictionRefined == "normal"
    for g, types in normal_groups.items():
        m = (norm & obs.CellType.isin(types)).values
        lab[m] = g
    return lab


def cpm(X) -> sp.csr_matrix:
    X = sp.csr_matrix(X, dtype=np.float64)
    tot = np.asarray(X.sum(1)).ravel()
    return sp.diags(1e6 / np.maximum(tot, 1)) @ X


def build_reference(X_cpm, genes, labels: pd.Series, n_markers: int, min_mean_cpm: float) -> pd.DataFrame:
    """Mean CPM per group, restricted to the union of each group's top marker genes."""
    groups = [g for g in pd.unique(labels.dropna())]
    prof = {}
    for g in groups:
        m = (labels == g).values
        prof[g] = np.asarray(X_cpm[m].mean(0)).ravel()
    P = pd.DataFrame(prof, index=genes)
    markers = set()
    for g in groups:
        others = P.drop(columns=g).mean(1)
        lfc = np.log2((P[g] + 1) / (others + 1))
        ok = P[g] >= min_mean_cpm
        markers |= set(lfc[ok].nlargest(n_markers).index)
    return P.loc[sorted(markers)]


def nnls_fractions(bulk: pd.DataFrame, ref: pd.DataFrame) -> pd.DataFrame:
    """bulk: genes x samples in linear (CPM-like) units. Returns samples x types fractions summing to 1."""
    genes = ref.index.intersection(bulk.index)
    R = ref.loc[genes].values
    scale = R.mean(1, keepdims=True) + 1e-9      # equalize gene weight (otherwise a few huge genes dominate)
    Rs = R / scale
    B = bulk.loc[genes].values / scale
    out = np.array([nnls(Rs, B[:, j])[0] for j in range(B.shape[1])])
    out = out / np.maximum(out.sum(1, keepdims=True), 1e-12)
    return pd.DataFrame(out, index=bulk.columns, columns=ref.columns)


def malignant_shares(frac: pd.DataFrame, like_types) -> pd.DataFrame:
    m = frac[[c for c in like_types if c in frac.columns]]
    return m.div(m.sum(1).replace(0, np.nan), axis=0)
