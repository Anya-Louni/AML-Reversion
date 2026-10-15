"""Linear expression-shift simulation: x' = x + alpha * s * sd.

This is a visualization of a hypothesis, not a model of cell fate dynamics. Each malignant cell's
log expression is moved along the drug's LINCS consensus signature (z-scores, converted to the
reference's per-gene SD units), with alpha rising linearly over the frames, and every frame is
placed on the landscape with the fixed Phase 1 map.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from pipeline.phase1_landscape.mapping import LandscapeMapper


def mutation_tags(obs: pd.DataFrame) -> dict[str, list[str]]:
    """Mutations detected in the authors' single-cell genotyping (MutTranscripts), per patient.

    Entries look like 'NPM1.W288fs/21' (gene.variant/read id); a tag is listed if any cell carries it.
    """
    out = {}
    for p, s in obs.dropna(subset=["MutTranscripts"]).groupby("patient", observed=True).MutTranscripts:
        tags = set()
        for v in s.astype(str):
            for tok in re.split(r"[\s,;]+", v):
                m = re.match(r"^([A-Za-z0-9-]+)\.([^/]+)", tok)
                if m:
                    tags.add(f"{m.group(1)} {m.group(2)}")
        out[str(p)] = sorted(tags)
    return out


def drug_shift(mapper: LandscapeMapper, drug_sig: pd.Series, gene_sd: np.ndarray) -> np.ndarray:
    """Per model gene shift for alpha = 1; genes without a reliable LINCS value are not moved."""
    s = drug_sig.reindex(mapper.genes).fillna(0.0).values
    return s * gene_sd


def simulate(mapper: LandscapeMapper, X0: np.ndarray, shift: np.ndarray, alphas: np.ndarray) -> np.ndarray:
    """X0: cells x model genes (log-normalized). Returns frames x cells x (score, x, y)."""
    T, N = len(alphas), X0.shape[0]
    frames = np.clip(X0[None, :, :] + alphas[:, None, None] * shift[None, None, :], 0, None)
    P = mapper.knn.predict(mapper.embed(frames.reshape(T * N, -1)))
    return P.reshape(T, N, 3)


def offmanifold_fraction(mapper: LandscapeMapper, X: np.ndarray, cut: float) -> float:
    """Share of cells whose mean distance to their reference neighbours exceeds the Phase 2 cut."""
    d, _ = mapper.knn.kneighbors(mapper.embed(X))
    return float(np.mean(d.mean(1) > cut))


def outcome(traj: np.ndarray, threshold: float) -> dict:
    s0, s_end, s_max = traj[0, :, 0], traj[-1, :, 0], traj[:, :, 0].max(0)
    stuck0 = s0 < threshold
    matured = stuck0 & (s_end >= threshold)
    resistant = stuck0 & (s_max < threshold)
    return {"n": int(len(s0)), "already_mature": int((~stuck0).sum()), "stuck_start": int(stuck0.sum()),
            "matured": int(matured.sum()), "resistant": int(resistant.sum()),
            "transient": int((stuck0 & ~matured & ~resistant).sum()),
            "pct_matured_of_stuck": float(matured.sum() / max(stuck0.sum(), 1) * 100),
            "mean_score_start": float(s0.mean()), "mean_score_end": float(s_end.mean()),
            "status": np.where(~stuck0, 0, np.where(matured, 1, np.where(resistant, 2, 3))).astype(np.uint8)}
