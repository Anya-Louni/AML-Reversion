"""Phase 4: validation against Beat AML ex vivo drug response (plan: reports/phase4_plan.md)."""
from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests

from pipeline.common import get_logger, path
from pipeline.phase3_drugs import scoring as S3
from pipeline.phase4_validation import models as M
from pipeline.phase4_validation.report import write_report

log = get_logger("phase4")
MIN_DRUGS = 20
N_PERM = 10000
N_BOOT = 2000


def load(cfg, proc):
    r = pd.read_parquet(proc / "beataml_drug_response.parquet").dropna(subset=["dbgap_rnaseq_sample"])
    r["sens"] = -r.auc
    Y = r.pivot_table(index="dbgap_rnaseq_sample", columns="inhibitor", values="sens", aggfunc="mean")
    R = pd.read_parquet(proc / "reversion_cosine_beataml.parquet")
    samples = Y.index.intersection(R.index)
    counts = Y.loc[samples].notna().sum()
    drugs = [d for d in R.columns if d in counts.index and counts[d] >= 50]
    subj = r.drop_duplicates("dbgap_rnaseq_sample").set_index("dbgap_rnaseq_sample").dbgap_subject_id.loc[samples].astype(str)
    return Y.loc[samples, drugs], R.loc[samples, drugs], subj, r


def features(proc, samples):
    dec = pd.read_parquet(proc / "beataml_deconvolution.parquet").loc[samples]
    frac = dec[[c for c in dec.columns if c.startswith("share_")] + ["bulk_maturity_index"]].fillna(0)
    clin = pd.read_parquet(proc / "beataml_clinical.parquet")
    clin = clin[clin.dbgap_rnaseq_sample.isin(samples)].drop_duplicates("dbgap_rnaseq_sample").set_index("dbgap_rnaseq_sample")
    mut = pd.read_parquet(proc / "beataml_mutations_matrix.parquet")
    dna = clin.dbgap_dnaseq_sample.reindex(samples).replace({"nan": np.nan, "None": np.nan})
    has = dna.isin(mut.index)
    genes = [c for c in mut.columns if c not in ("capture_type",) and not c.startswith("clin_")]
    sub = mut.loc[dna[has].values, genes]
    sub.index = dna[has].index
    keep = [g for g in genes if sub[g].mean() >= 0.05 and not g.endswith("_hotspot")]
    X = sub[keep].reindex(samples).fillna(0)
    X["missing_dna"] = (~has).astype(int).values
    cap = mut.capture_type.reindex(dna.values).values
    X["targeted_panel"] = pd.Series(cap, index=samples).fillna("").str.contains("Nimblegen").astype(int)
    for c in ["FLT3-ITD", "NPM1"]:
        X[f"clin_{c}"] = clin[c].reindex(samples).map({"positive": 1, "negative": 0}).fillna(0).values
    hot = mut[["IDH1_hotspot", "IDH2_hotspot"]].reindex(dna.values)
    hot.index = samples
    return frac, X.astype(float), hot, dec


def run(cfg: dict) -> None:
    seed = cfg["seed"]
    proc = path(cfg, "processed", "x").parent
    Y, R, subj, resp = load(cfg, proc)
    frac, MUT, hot, dec = features(proc, Y.index)
    folds = M.patient_folds(subj, 5, seed)
    st = {"n_samples": len(Y), "n_patients": int(subj.nunique()), "n_drugs": Y.shape[1],
          "n_pairs": int(Y.notna().sum().sum()), "mutation_features": list(MUT.columns),
          "drugs_per_sample_median": float(Y.notna().sum(1).median())}
    log.info("%d samples, %d patients, %d drugs", len(Y), subj.nunique(), Y.shape[1])

    specs = {"B0 random": ("random", None, None), "B1 global drug ranking": ("global", None, None),
             "B2 mutation-only": ("mutation", None, MUT), "B3 differentiation fractions": ("fractions", None, frac),
             "M1 reversion score (raw)": ("raw", R, None), "M1c reversion score (calibrated)": ("calibrated", R, None),
             "M2 combined": ("combined", R, pd.concat([frac, MUT], axis=1)),
             # Exploratory ablation, added after seeing the pre-registered results (not used for claims):
             "X1 fractions + mutations (ablation)": ("fractions", None, pd.concat([frac, MUT], axis=1))}
    res, preds = {}, {}
    for name, (kind, Rr, F) in specs.items():
        if kind == "random":
            draws = [M.metrics(*M.fit_predict("random", Y, folds, seed=seed + 1000 * i), Y, folds, MIN_DRUGS) for i in range(100)]
            met = {k: pd.concat([d[k] for d in draws], axis=1).mean(1) for k in ("per_sample_raw", "per_sample_centred", "per_drug")}
        else:
            P, Ptr = M.fit_predict(kind, Y, folds, R=Rr, F=F, seed=seed)
            met = M.metrics(P, Ptr, Y, folds, MIN_DRUGS)
            preds[name] = (P, met["Pc"])
        res[name] = met
        log.info("%s done", name)

    rows = []
    for name, met in res.items():
        a = M.boot_ci(met["per_sample_raw"], subj, N_BOOT, seed)
        b = M.boot_ci(met["per_sample_centred"], subj, N_BOOT, seed)
        c = M.boot_ci_simple(met["per_drug"], N_BOOT, seed)
        pd_ = met["per_drug"].dropna()
        rows.append({"model": name, "raw_mean": a[0], "raw_lo": a[1], "raw_hi": a[2],
                     "centred_mean": b[0], "centred_lo": b[1], "centred_hi": b[2],
                     "per_drug_mean": c[0], "per_drug_lo": c[1], "per_drug_hi": c[2],
                     "n_samples_scored": int(met["per_sample_raw"].notna().sum())})
    table = pd.DataFrame(rows).set_index("model")
    st["table"] = table

    # ---------- primary tests (M1 raw)
    stat1, p1, null1 = M.perm_test(R, Y, MIN_DRUGS, N_PERM, seed)
    Rc = M.centre(R, Y, folds, {f: R[folds != f] for f in np.unique(folds)})
    Yc = M.centre(Y, Y, folds, None)
    stat2, p2, null2 = M.perm_test(Rc, Yc, MIN_DRUGS, N_PERM, seed + 1)
    st["P1"] = {"statistic": stat1, "p": p1, "significant": p1 < 0.025}
    st["P2"] = {"statistic": stat2, "p": p2, "significant": p2 < 0.025}
    st["nulls"] = (null1, null2)

    # ---------- paired comparisons vs B1 and between models
    def paired(a, b, metric):
        x = pd.concat([res[a][metric], res[b][metric]], axis=1).dropna()
        d = x.iloc[:, 0] - x.iloc[:, 1]
        if (d == 0).all():
            return {"median_diff": 0.0, "p": 1.0}
        return {"mean_diff": float(d.mean()), "p": float(wilcoxon(d, alternative="greater").pvalue)}
    comps = {}
    for m in ["M1 reversion score (raw)", "M1c reversion score (calibrated)", "M2 combined", "B2 mutation-only", "B3 differentiation fractions"]:
        comps[f"{m} > B1 (raw)"] = paired(m, "B1 global drug ranking", "per_sample_raw")
    for b in ["B2 mutation-only", "B3 differentiation fractions"]:
        comps[f"M2 combined > {b} (centred)"] = paired("M2 combined", b, "per_sample_centred")
    comps["M1 raw > B0 (centred)"] = paired("M1 reversion score (raw)", "B0 random", "per_sample_centred")
    st["exploratory"] = pd.DataFrame({
        "M2 combined > X1 ablation (raw)": paired("M2 combined", "X1 fractions + mutations (ablation)", "per_sample_raw"),
        "M2 combined > X1 ablation (centred)": paired("M2 combined", "X1 fractions + mutations (ablation)", "per_sample_centred"),
        "M1c calibrated > B0 (centred)": paired("M1c reversion score (calibrated)", "B0 random", "per_sample_centred"),
    }).T
    st["comparisons"] = pd.DataFrame(comps).T

    # Claim rules from the plan. (The plan's second clause, M2 > B1 and M2 > B2/B3, implies M2 > B1, so it is covered.)
    m1_beats_b1 = comps["M1 reversion score (raw) > B1 (raw)"]["p"] < 0.05
    m2_beats_b1 = comps["M2 combined > B1 (raw)"]["p"] < 0.05
    st["claims"] = {"information_about_response": bool(st["P1"]["significant"]),
                    "patient_specific_value": bool(st["P2"]["significant"] and (m1_beats_b1 or m2_beats_b1)),
                    "m1_beats_b1": bool(m1_beats_b1), "m2_beats_b1": bool(m2_beats_b1)}

    # ---------- per-drug detail for M1
    pdrug = []
    for d in Y.columns:
        ok = Y[d].notna()
        rho, p = spearmanr(R.loc[ok, d], Y.loc[ok, d])
        pdrug.append({"drug": d, "n": int(ok.sum()), "spearman": rho, "p": p})
    pdrug = pd.DataFrame(pdrug)
    pdrug["q"] = multipletests(pdrug.p.fillna(1), method="fdr_bh")[1]
    meta = pd.read_parquet(proc / "drug_signature_meta.parquet")
    pdrug = pdrug.join(meta[["tier", "cells"]], on="drug").sort_values("spearman", ascending=False)
    st["per_drug_m1"] = pdrug

    # ---------- secondary: other score variants
    sig_all = ad.read_h5ad(proc / "lincs_signatures_all.h5ad", backed="r")
    fs = sig_all.var.feature_space
    D = pd.read_parquet(proc / "drug_consensus_signatures.parquet")
    stuck = pd.read_parquet(proc / "stuck_signatures_beataml.parquet").loc[Y.index]
    sig_full = pd.read_csv(proc / "maturation_signature_full.csv")
    variants = {}
    lm = set(sig_all.var_names[fs.eq("landmark").values])
    variants["landmark genes"] = S3.score_patients(D.loc[Y.columns], -stuck[[g for g in stuck.columns if g in lm]], 100)[0]
    variants["connectivity score"] = pd.read_parquet(proc / "reversion_connectivity_beataml.parquet").loc[Y.index, Y.columns]
    # full signature (with ribosomal genes): rebuild bulk stuck for those genes
    be = ad.read_h5ad(proc / "beataml_expr.h5ad")
    E = pd.DataFrame(be[Y.index].X, index=Y.index, columns=be.var_names)
    from pipeline.phase2_patients.stuck import bulk_mature_expectation, bulk_stuck
    genes_full = [g for g in sig_full.gene if g in E.columns and g in D.columns
                  and fs.reindex([g]).isin(["landmark", "best inferred"]).iloc[0]]
    fit = bulk_mature_expectation(E[genes_full], dec.loc[Y.index, "mature_share"])
    variants["full signature (with ribosomal genes)"] = S3.score_patients(D.loc[Y.columns], -bulk_stuck(E[genes_full], fit.expected_mature), 100)[0]
    my = [d for d in Y.columns if meta.tier.get(d) == "AML/monocytic"]
    variants["AML/monocytic-line drugs only"] = R[my]
    sec = []
    for name, V in variants.items():
        Yv = Y[V.columns]
        stat_a, p_a, _ = M.perm_test(V, Yv, MIN_DRUGS, 2000, seed)
        Vc = M.centre(V, Yv, folds, {f: V[folds != f] for f in np.unique(folds)})
        stat_b, p_b, _ = M.perm_test(Vc, M.centre(Yv, Yv, folds, None), MIN_DRUGS, 2000, seed + 1)
        sec.append({"variant": name, "n_drugs": V.shape[1], "raw_mean_rho": stat_a, "raw_p": p_a, "centred_mean_rho": stat_b, "centred_p": p_b})
    st["secondary"] = pd.DataFrame(sec)

    # ---------- IDH subgroup check
    idh = []
    for drug, gene in [("AGI-5198", "IDH1_hotspot"), ("AGI-6780", "IDH2_hotspot")]:
        if drug not in Y.columns:
            continue
        ok = Y[drug].notna() & hot[gene].notna()
        mutm, wt = ok & (hot[gene] == 1), ok & (hot[gene] == 0)
        u = mannwhitneyu(Y.loc[mutm, drug], Y.loc[wt, drug], alternative="greater")
        u2 = mannwhitneyu(R.loc[mutm, drug], R.loc[wt, drug], alternative="two-sided")
        idh.append({"drug": drug, "mutation": gene.replace("_hotspot", " hotspot"), "n_mut": int(mutm.sum()), "n_wt": int(wt.sum()),
                    "median_sens_mut": float(Y.loc[mutm, drug].median()), "median_sens_wt": float(Y.loc[wt, drug].median()),
                    "p_mut_more_sensitive": float(u.pvalue),
                    "median_reversion_mut": float(R.loc[mutm, drug].median()), "median_reversion_wt": float(R.loc[wt, drug].median()),
                    "p_reversion_differs": float(u2.pvalue)})
    st["idh"] = pd.DataFrame(idh)
    st["idh_data"] = (Y, hot)

    # ---------- venetoclax sanity check
    if "Venetoclax" in Y.columns:
        ok = Y.Venetoclax.notna()
        rho, p = spearmanr(dec.loc[Y.index[ok], "mature_share"], Y.loc[ok, "Venetoclax"])
        st["venetoclax"] = {"n": int(ok.sum()), "spearman": float(rho), "p": float(p)}

    # ---------- save
    table.to_csv(proc / "phase4_model_table.csv")
    pdrug.to_csv(proc / "phase4_per_drug_m1.csv", index=False)
    st["idh"].to_csv(proc / "phase4_idh.csv", index=False)
    res["M1 reversion score (raw)"]["per_sample_raw"].to_frame("rho").join(subj).to_csv(proc / "phase4_per_sample_m1.csv")
    ser = {k: v for k, v in st.items() if isinstance(v, (dict, int, float))}
    ser["table"] = table.reset_index().to_dict(orient="records")
    (proc / "phase4_stats.json").write_text(json.dumps(ser, indent=2, default=str))
    write_report(cfg, st, res, Y, R, dec)
    log.info("phase 4 done")
