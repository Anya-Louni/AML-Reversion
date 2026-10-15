"""Experiment 3, F1: compare the linear simulation's change vector with real measured change.

The simulated change is exactly what Phase 5 adds to a cell: normalize(s_g * sd_g) * target_norm,
over the landscape mapper's genes that also have a reliable LINCS measurement. Genes outside that
set are not moved by the simulation and are excluded from the comparison (reported as a count).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# Negative controls fixed in reports/exp3_f1_plan.md before the analysis was run.
NEGATIVE_CONTROLS = ["ciprofloxacin", "famotidine", "loratadine", "paracetamol",
                     "lisinopril", "furosemide", "glipizide", "lidocaine"]


class Simulator:
    """The Phase 5 linear shift, exposed as a gene-space change vector."""

    def __init__(self, mapper, D: pd.DataFrame, gene_sd: np.ndarray, reliable: set, target_norm: float, alpha_max: float):
        genes = pd.Index(mapper.genes)
        self.keep = genes.isin(list(reliable)) & genes.isin(D.columns)
        self.genes = genes[self.keep]
        self.sd = pd.Series(gene_sd, index=genes)[self.keep]
        self.D = D[self.genes]
        self.target_norm = target_norm
        self.alpha_max = alpha_max

    def change(self, drug: str, weight_by_sd: bool = True) -> pd.Series:
        """Total change the simulation applies at the final frame (alpha = alpha_max)."""
        s = self.D.loc[drug]
        v = s * self.sd if weight_by_sd else s
        n = np.linalg.norm(v.values)
        return v / max(n, 1e-9) * self.target_norm * self.alpha_max

    def change_matrix(self, drugs: list[str], weight_by_sd: bool = True) -> pd.DataFrame:
        V = self.D.loc[drugs]
        if weight_by_sd:
            V = V * self.sd.values
        n = np.linalg.norm(V.values, axis=1, keepdims=True)
        return pd.DataFrame(V.values / np.maximum(n, 1e-9) * self.target_norm * self.alpha_max,
                            index=drugs, columns=self.genes)


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / max(na * nb, 1e-12))


def direction_test(sim: Simulator, real: pd.Series, drug: str, null_drugs: list[str],
                   weight_by_sd: bool = True, n_perm: int = 1000, seed: int = 0) -> dict:
    """Cosine and Spearman of simulated vs real change, with the rank of `drug` in the null."""
    shared = sim.genes.intersection(real.index)
    r = real.loc[shared].values
    # Score the null set plus the tested drug (which may sit outside the null, e.g. a non-myeloid-tier drug).
    scored = list(dict.fromkeys(list(null_drugs) + [drug]))
    V = sim.change_matrix(scored, weight_by_sd)[shared]
    cos_all = V.values @ r / np.maximum(np.linalg.norm(V.values, axis=1) * np.linalg.norm(r), 1e-12)
    cos_all = pd.Series(cos_all, index=scored)
    null = cos_all[null_drugs]
    c = float(cos_all[drug])
    rank = int((null > c).sum() + 1)
    p = sim.change(drug, weight_by_sd).loc[shared].values
    rho = float(spearmanr(p, r)[0])
    rng = np.random.default_rng(seed)
    perm = np.array([cosine(p[rng.permutation(len(p))], r) for _ in range(n_perm)])
    return {"n_genes": int(len(shared)), "cosine": c, "spearman": rho,
            "rank_in_null": rank, "n_null": len(null_drugs), "percentile": float((null < c).mean() * 100),
            "p_rank": rank / len(null_drugs), "null_mean": float(null.mean()), "null_sd": float(null.std()),
            "perm_p": float((np.sum(perm >= c) + 1) / (n_perm + 1)),
            "negative_controls": {d: float(cos_all[d]) for d in NEGATIVE_CONTROLS if d in cos_all.index},
            "cos_all": null}


def signature_shift(change: pd.Series, sig: pd.DataFrame) -> float:
    """Maturation direction of a change vector: mean change of up genes minus mean of down genes."""
    up = sig.loc[sig.direction.str.startswith("up"), "gene"]
    dn = sig.loc[sig.direction.str.startswith("down"), "gene"]
    u = change.reindex(up).dropna()
    d = change.reindex(dn).dropna()
    if len(u) == 0 or len(d) == 0:
        return float("nan")
    return float(u.mean() - d.mean())


def maturity_test(sim: Simulator, real: pd.Series, drug: str, sig: pd.DataFrame, ln2: float = np.log(2)) -> dict:
    """Direction and magnitude of the maturation shift, real vs simulated, on shared genes."""
    shared = sim.genes.intersection(real.index)
    simc = sim.change(drug).loc[shared] / ln2          # natural-log units -> log2, to match log2FC
    realc = real.loc[shared]
    rs, ss = signature_shift(realc, sig), signature_shift(simc, sig)
    return {"real_shift": rs, "sim_shift": ss, "same_direction": bool(np.sign(rs) == np.sign(ss)),
            "magnitude_ratio": float(ss / rs) if rs != 0 else float("nan"),
            "n_shared_genes": int(len(shared))}


def landscape_scores(mapper, counts: pd.DataFrame, meta: pd.DataFrame, arm: str, control: str,
                     sim_change: pd.Series | None, subset: dict | None = None) -> dict:
    """Phase 1 landscape score of bulk profiles (platform transfer, reported with that caveat)."""
    from pipeline.phase1_landscape.mapping import lognorm
    m = meta
    if subset:
        for k, v in subset.items():
            m = m[m[k] == v]
    tr, ct = list(m.loc[m.arm == arm, "sample"]), list(m.loc[m.arm == control, "sample"])
    X = counts[ct + tr]
    L = pd.DataFrame(np.asarray(lognorm(X.T.values).todense()), index=ct + tr, columns=counts.index)
    pred = mapper.predict(L.values, L.columns)
    pred.index = L.index
    out = {"control": float(pred.loc[ct, "score"].mean()), "treated": float(pred.loc[tr, "score"].mean())}
    out["real_delta"] = out["treated"] - out["control"]
    if sim_change is not None:
        Lsim = L.loc[ct].copy()
        common = Lsim.columns.intersection(sim_change.index)
        Lsim[common] = Lsim[common].values + sim_change[common].values
        ps = mapper.predict(Lsim.clip(lower=0).values, Lsim.columns)
        out["simulated"] = float(ps["score"].mean())
        out["sim_delta"] = out["simulated"] - out["control"]
    return out


def trajectory_test(sim: Simulator, reals: dict[str, pd.Series], drug: str) -> pd.DataFrame:
    """Project each real timepoint onto the simulated direction; report the off-axis residual."""
    rows = []
    for t, real in reals.items():
        shared = sim.genes.intersection(real.index)
        u = sim.change(drug).loc[shared].values
        u = u / max(np.linalg.norm(u), 1e-12)
        r = real.loc[shared].values
        a = float(r @ u)
        resid = r - a * u
        rows.append({"timepoint": t, "n_genes": len(shared), "parallel": a,
                     "residual_norm": float(np.linalg.norm(resid)), "total_norm": float(np.linalg.norm(r)),
                     "off_axis_fraction": float(np.linalg.norm(resid) / max(np.linalg.norm(r), 1e-12)),
                     "cosine": cosine(r, u)})
    return pd.DataFrame(rows)
