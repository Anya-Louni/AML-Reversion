"""Probit dose-response fits on Beat AML raw viability, reproducing the published AUC definition.

Published AUCs range 0 to 286.3 = 100 * (log10(10) - log10(0.0137)): the area under the fitted
viability curve (percent) over log10 concentration. We fit viability/100 ~ Phi(a + b * log10(conc)) by
minimizing binomial deviance (viability clipped to [0, 1]) and integrate the fitted curve.
"""
from __future__ import annotations

import numpy as np
from scipy.special import ndtr
from scipy.optimize import minimize

_GRID = 200


def _deviance(params, x, y):
    a, b = params
    p = np.clip(ndtr(a + b * x), 1e-9, 1 - 1e-9)
    return -2 * np.sum(y * np.log(p) + (1 - y) * np.log(1 - p))


def fit_probit(conc: np.ndarray, viability: np.ndarray) -> tuple[float, float]:
    x = np.log10(np.asarray(conc, float))
    y = np.clip(np.asarray(viability, float) / 100.0, 0, 1)
    res = minimize(_deviance, x0=np.array([0.5, -0.5]), args=(x, y), method="L-BFGS-B",
                   bounds=[(-20, 20), (-20, 20)])
    return float(res.x[0]), float(res.x[1])


def auc(a: float, b: float, cmin: float, cmax: float) -> float:
    xs = np.linspace(np.log10(cmin), np.log10(cmax), _GRID)
    return float(np.trapezoid(100 * ndtr(a + b * xs), xs))


def curve_auc(conc, viability) -> float:
    conc = np.asarray(conc, float)
    a, b = fit_probit(conc, viability)
    return auc(a, b, conc.min(), conc.max())
