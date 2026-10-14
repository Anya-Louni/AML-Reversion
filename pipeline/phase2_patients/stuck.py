"""Per-patient 'stuck' signatures restricted to maturation-signature genes.

stuck = patient expression - expression expected for mature myeloid cells.
push  = -stuck, the direction a drug would have to move the patient toward maturity.

Single-cell patients (van Galen): mean log-normalized expression of malignant
cells minus the mean of normal reference monocytes (same platform).

Bulk patients (Beat AML): no mature reference exists on the bulk platform, so the
mature expectation is estimated within the cohort. For every gene, regress
expression on the deconvolved mature share of the malignant compartment and
predict at share = 1. This keeps both sides of the subtraction on the same platform.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def sc_stuck(X_lognorm_patient: np.ndarray, mono_mean: np.ndarray, genes) -> pd.Series:
    return pd.Series(np.asarray(X_lognorm_patient).mean(0) - mono_mean, index=genes)


def bulk_mature_expectation(expr: pd.DataFrame, mature_share: pd.Series) -> pd.DataFrame:
    """expr: samples x genes (log2 CPM). Returns per-gene intercept, slope, expected mature value, r."""
    m = mature_share.loc[expr.index].values
    X = expr.values
    mc = m - m.mean()
    slope = (mc[:, None] * (X - X.mean(0))).sum(0) / (mc ** 2).sum()
    inter = X.mean(0) - slope * m.mean()
    r = (mc[:, None] * (X - X.mean(0))).sum(0) / np.sqrt((mc ** 2).sum() * ((X - X.mean(0)) ** 2).sum(0))
    return pd.DataFrame({"intercept": inter, "slope": slope, "expected_mature": inter + slope, "r": r}, index=expr.columns)


def bulk_stuck(expr: pd.DataFrame, expected: pd.Series) -> pd.DataFrame:
    return expr - expected.loc[expr.columns].values
