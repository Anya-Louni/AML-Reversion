"""Experiment 2 data assembly, models and metrics (plan: reports/exp2_plan.md)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.linear_model import Ridge

ALPHAS = [0.01, 0.1, 1, 10, 100, 1000]
MIN_DRUGS = 20


# ------------------------------------------------------------------ target
def drug_means(Y: pd.DataFrame, dev: pd.Index) -> pd.Series:
    return Y.loc[dev].mean()


def residuals(Y: pd.DataFrame, means: pd.Series) -> pd.DataFrame:
    return Y - means


# ------------------------------------------------------------------ features
def mutation_clinical_features(samples: pd.Index, dev: pd.Index, clin: pd.DataFrame, mut: pd.DataFrame) -> pd.DataFrame:
    """Mutation genes chosen on development samples only (>= 5% mutated among dev samples with DNA)."""
    c = clin.reindex(samples)
    dna = c.dbgap_dnaseq_sample.where(c.dbgap_dnaseq_sample.isin(mut.index))
    genes = [g for g in mut.columns if g != "capture_type" and not g.startswith("clin_") and not g.endswith("_hotspot")]
    M = pd.DataFrame(0.0, index=samples, columns=genes)
    has = dna.notna()
    M.loc[has[has].index] = mut.loc[dna[has].values, genes].values
    dev_has = has.loc[dev]
    keep = [g for g in genes if M.loc[dev_has[dev_has].index, g].mean() >= 0.05]
    X = M[keep].copy()
    X["missing_dna"] = (~has).astype(float).values
    cap = mut.capture_type.reindex(dna.values).fillna("").values
    X["targeted_panel"] = pd.Series(cap, index=samples).str.contains("Nimblegen").astype(float)
    for col in ["FLT3-ITD", "NPM1"]:
        X[f"clin_{col}"] = c[col].map({"positive": 1.0, "negative": 0.0}).fillna(0.0).values
    for t in ["Bone Marrow Aspirate", "Peripheral Blood", "Leukapheresis"]:
        X[f"specimen_{t.split()[0].lower()}"] = (c.specimenType == t).astype(float).values
    age = pd.to_numeric(c.ageAtSpecimenAcquisition, errors="coerce")
    X["age"] = age.fillna(age.loc[dev].median()).values
    return X


# ------------------------------------------------------------------ fitting
class PerDrugRidge:
    """One ridge model per drug: residual(p, d) ~ features(p) [+ drug-specific feature], no intercept,
    features standardized with development statistics."""

    def __init__(self, alpha: float):
        self.alpha = alpha
        self.models: dict = {}
        self.stats: dict = {}

    @staticmethod
    def _design(F: pd.DataFrame, extra: pd.DataFrame | None, d: str) -> pd.DataFrame:
        X = F.copy()
        if extra is not None:
            X["__drug_feature"] = extra[d].reindex(F.index).values
        return X

    def fit(self, R: pd.DataFrame, F: pd.DataFrame, rows: pd.Index, extra: pd.DataFrame | None = None) -> "PerDrugRidge":
        for d in R.columns:
            X = self._design(F.loc[rows], extra, d)
            y = R.loc[rows, d]
            ok = y.notna() & X.notna().all(1)
            if ok.sum() < 10:
                continue
            mu, sd = X[ok].mean(), X[ok].std().replace(0, 1)
            m = Ridge(alpha=self.alpha, fit_intercept=False).fit(((X[ok] - mu) / sd).values, y[ok].values)
            self.models[d], self.stats[d] = m, (mu, sd)
        return self

    def predict(self, F: pd.DataFrame, rows: pd.Index, extra: pd.DataFrame | None = None) -> pd.DataFrame:
        P = pd.DataFrame(np.nan, index=rows, columns=list(self.models))
        for d, m in self.models.items():
            X = self._design(F.loc[rows], extra, d)
            mu, sd = self.stats[d]
            ok = X.notna().all(1)
            P.loc[ok[ok].index, d] = m.predict(((X[ok] - mu) / sd).values)
        return P

    def pure_type_sensitivity(self, F: pd.DataFrame, means: pd.Series) -> pd.DataFrame:
        """sensitivity(k, d) = development drug mean + prediction for a sample made only of cell type k."""
        out = {}
        for d, m in self.models.items():
            mu, sd = self.stats[d]
            E = pd.DataFrame(np.eye(F.shape[1]), columns=F.columns, index=F.columns)
            out[d] = means[d] + m.predict(((E - mu) / sd).values)
        return pd.DataFrame(out, index=F.columns).T


def patient_folds(subjects: pd.Series, k: int = 5, seed: int = 0) -> pd.Series:
    u = np.array(sorted(subjects.unique()))
    np.random.default_rng(seed).shuffle(u)
    return subjects.map({s: i % k for i, s in enumerate(u)})


# ------------------------------------------------------------------ metrics
def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    if a.std() == 0 or b.std() == 0:
        return 0.0
    return float(np.corrcoef(rankdata(a), rankdata(b))[0, 1])


def per_sample_spearman(P: pd.DataFrame, R: pd.DataFrame) -> pd.Series:
    out = {}
    cols = P.columns.intersection(R.columns)
    for s in P.index:
        p, r = P.loc[s, cols], R.loc[s, cols]
        ok = p.notna() & r.notna()
        if ok.sum() >= MIN_DRUGS:
            out[s] = _spearman(p[ok].values, r[ok].values)
    return pd.Series(out, dtype=float)


def per_patient(values: pd.Series, subject: pd.Series) -> pd.Series:
    return values.groupby(subject.loc[values.index]).mean()


def top_k_hits(P: pd.DataFrame, R: pd.DataFrame, k: int = 5, top_frac: float = 0.2) -> pd.Series:
    out = {}
    cols = P.columns.intersection(R.columns)
    for s in P.index:
        p, r = P.loc[s, cols], R.loc[s, cols]
        ok = p.notna() & r.notna()
        if ok.sum() < MIN_DRUGS:
            continue
        p, r = p[ok], r[ok]
        top_obs = set(r.nlargest(max(1, int(round(top_frac * len(r))))).index)
        out[s] = len(set(p.nlargest(k).index) & top_obs) / k
    return pd.Series(out, dtype=float)


def pooled_pearson(P: pd.DataFrame, R: pd.DataFrame) -> float:
    cols = P.columns.intersection(R.columns)
    a, b = P[cols].values.ravel(), R.loc[P.index, cols].values.ravel()
    ok = ~np.isnan(a) & ~np.isnan(b)
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def bootstrap_ci(per_pat: pd.Series, n: int = 10000, seed: int = 0) -> tuple[float, float, float]:
    v = per_pat.dropna().values
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(v), size=(n, len(v)))
    m = v[idx].mean(1)
    return float(v.mean()), float(np.quantile(m, 0.025)), float(np.quantile(m, 0.975))


def sign_flip_test(diff: pd.Series, n: int = 10000, seed: int = 0, alternative: str = "greater") -> float:
    v = diff.dropna().values
    rng = np.random.default_rng(seed)
    signs = rng.choice([-1, 1], size=(n, len(v)))
    null = (signs * v).mean(1)
    obs = v.mean()
    if alternative == "greater":
        return float((np.sum(null >= obs) + 1) / (n + 1))
    return float((np.sum(np.abs(null) >= abs(obs)) + 1) / (n + 1))


def label_permutation_test(P: pd.DataFrame, R: pd.DataFrame, subject: pd.Series, n: int = 2000, seed: int = 0) -> float:
    """Within-sample shuffle of drug labels; statistic = mean over patients of per-sample Spearman."""
    rng = np.random.default_rng(seed)
    cols = P.columns.intersection(R.columns)
    rows, obs_rows = [], []
    for s in P.index:
        p, r = P.loc[s, cols], R.loc[s, cols]
        ok = (p.notna() & r.notna()).values
        if ok.sum() >= MIN_DRUGS:
            a, b = rankdata(p.values[ok]), rankdata(r.values[ok])
            a = (a - a.mean()) / (a.std() or 1); b = (b - b.mean()) / (b.std() or 1)
            rows.append((s, a, b))
    subj = subject.loc[[s for s, _, _ in rows]].values
    obs = pd.Series([(a * b).mean() for _, a, b in rows]).groupby(subj).mean().mean()
    null = np.empty(n)
    for i in range(n):
        vals = [(a[rng.permutation(len(a))] * b).mean() for _, a, b in rows]
        null[i] = pd.Series(vals).groupby(subj).mean().mean()
    return float((np.sum(null >= obs) + 1) / (n + 1))
