"""Experiment 2, Step 1: ceiling analysis. Writes reports/exp2_step1_ceiling.md."""
from __future__ import annotations

import datetime as dt
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pipeline.common import get_logger, path
from pipeline.exp2 import ceiling as C
from pipeline.phase4_validation.run import load as load_phase4

log = get_logger("exp2.step1")
# Verdict rule, fixed before the analysis was run (also stated in the report):
PROCEED_CEILING, BORDERLINE_CEILING, MIN_SIGNAL_SHARE = 0.5, 0.3, 0.10


def md(df, index=False):
    return df.to_markdown(index=index, floatfmt=".3f")


def run(cfg: dict) -> None:
    proc = path(cfg, "processed", "x").parent
    raw_dir = path(cfg, "raw", "beataml", "x").parent
    Ysens, _, subj, resp = load_phase4(cfg, proc)
    A = -Ysens                          # AUC, same response as Phase 4 (sign does not affect variances)
    st: dict = {"n_samples": A.shape[0], "n_drugs": A.shape[1], "n_obs": int(A.notna().sum().sum())}

    # ---- 1. additive model on the Phase 4 matrix
    fit = C.additive_fit(A)
    st["shares"] = C.variance_shares(A, fit)
    resid = fit["resid"]
    var_resid = float(np.nanvar(resid.values))

    # ---- 2. measurement noise from technical replicates
    cache = proc / "exp2_replicate_aucs.parquet"
    raw = pd.read_csv(raw_dir / "beataml_wv1to4_raw_inhibitor_v4_dbgap.txt", sep="\t")
    if cache.exists():
        reps = pd.read_parquet(cache)
    else:
        log.info("fitting replicate curves")
        reps = C.replicate_aucs(raw)
        reps.to_parquet(cache)
    pairs = reps.groupby(["dbgap_rnaseq_sample", "inhibitor"]).auc.agg(["var", "count", "first", "last"])
    pairs = pairs[pairs["count"] >= 2]
    st["n_replicate_pairs"] = len(pairs)
    st["n_replicate_samples"] = int(pairs.index.get_level_values(0).nunique())
    in_matrix = pairs[pairs.index.get_level_values(1).isin(A.columns)]
    sigma2_drug = in_matrix.groupby(level=1)["var"].mean()
    sigma2_all = float(in_matrix["var"].mean())
    st["sigma2_single_curve"] = sigma2_all
    st["drugs_with_replicates"] = int(sigma2_drug.size)
    ncur = C.n_curves_per_pair(raw)
    N = pd.DataFrame(np.nan, index=A.index, columns=A.columns)
    nc = ncur[ncur.index.get_level_values(0).isin(A.index) & ncur.index.get_level_values(1).isin(A.columns)]
    for (s, d), v in nc.items():
        N.at[s, d] = v
    N = N.where(A.notna()).fillna(1)
    sig = pd.Series(sigma2_all, index=A.columns)
    sig.loc[sigma2_drug.index] = sigma2_drug[sigma2_drug.index.isin(sig.index)]
    noise = (pd.DataFrame(np.tile(sig.values, (len(A), 1)), index=A.index, columns=A.columns) / N).where(A.notna())
    mean_noise = float(np.nanmean(noise.values))
    R = 1 - mean_noise / var_resid
    st["reliability_replicates"] = {"var_residual": var_resid, "mean_noise_var": mean_noise, "R": R,
                                   "ceiling_corr": float(np.sqrt(max(R, 0))),
                                   "signal_share_of_total": float((var_resid - mean_noise) / st["shares"]["total_var"])}
    Rp = 1 - noise.mean(1) / resid.var(1)
    ceil_p = np.sqrt(Rp.clip(lower=0))
    st["per_sample_ceiling"] = {"mean": float(ceil_p.mean()), "median": float(ceil_p.median()),
                                "q10": float(ceil_p.quantile(0.1)), "q90": float(ceil_p.quantile(0.9))}

    # ---- 3. split-replicate check: two independent curves of the same sample and drug
    mu, row, col = fit["mu"], fit["row"], fit["col"]
    pr = in_matrix.reset_index()
    pr = pr[pr.dbgap_rnaseq_sample.isin(A.index)]
    pr["resA"] = pr["first"] - mu - pr.dbgap_rnaseq_sample.map(row) - pr.inhibitor.map(col)
    pr["resB"] = pr["last"] - mu - pr.dbgap_rnaseq_sample.map(row) - pr.inhibitor.map(col)
    RA = pr.pivot(index="dbgap_rnaseq_sample", columns="inhibitor", values="resA")
    RB = pr.pivot(index="dbgap_rnaseq_sample", columns="inhibitor", values="resB")
    split = C.per_row_spearman(RA, RB, min_n=10)
    st["split_replicate"] = {"n_samples": int(split.size), "mean_spearman": float(split.mean()) if len(split) else None,
                             "median_spearman": float(split.median()) if len(split) else None,
                             "drugs_per_sample_median": float(RA.notna().sum(1).loc[split.index].median()) if len(split) else None}
    st["pearson_replicate_auc"] = float(np.corrcoef(pr["first"], pr["last"])[0, 1])
    # Bias check: the effects above were fitted on published AUCs that include these same curves, which removes
    # some of the shared signal. Refit the effects with the replicated pairs masked out.
    A2 = A.copy()
    for s_, d_ in zip(pr.dbgap_rnaseq_sample, pr.inhibitor):
        A2.at[s_, d_] = np.nan
    f2 = C.additive_fit(A2)
    base2 = f2["mu"] + pr.dbgap_rnaseq_sample.map(f2["row"]) + pr.inhibitor.map(f2["col"])
    split2 = C.per_row_spearman(pr.assign(r=pr["first"] - base2).pivot(index="dbgap_rnaseq_sample", columns="inhibitor", values="r"),
                                pr.assign(r=pr["last"] - base2).pivot(index="dbgap_rnaseq_sample", columns="inhibitor", values="r"), 10)
    st["split_replicate_debiased"] = {"mean_spearman": float(split2.mean()), "ceiling": float(np.sqrt(max(split2.mean(), 0)))}
    st["n_samples_resid_below_noise"] = int((resid.var(axis=1) < st["sigma2_single_curve"]).sum())

    # ---- 4. test-retest across samples of the same patient (pessimistic: includes biological change)
    allr = resp.dropna(subset=["dbgap_rnaseq_sample"])
    allr = allr[allr.inhibitor.isin(A.columns)]
    Afull = allr.pivot_table(index="dbgap_rnaseq_sample", columns="inhibitor", values="auc", aggfunc="mean")
    ffull = C.additive_fit(Afull)
    sub_of = allr.drop_duplicates("dbgap_rnaseq_sample").set_index("dbgap_rnaseq_sample").dbgap_subject_id
    tr = []
    for s, grp in sub_of.groupby(sub_of):
        samples = list(grp.index)
        for i in range(len(samples)):
            for j in range(i + 1, len(samples)):
                v = C.per_row_spearman(ffull["resid"].loc[[samples[i]]], ffull["resid"].loc[[samples[j]]].set_axis([samples[i]]), 20)
                if len(v):
                    tr.append({"subject": s, "a": samples[i], "b": samples[j], "spearman": float(v.iloc[0])})
    tr = pd.DataFrame(tr)
    st["test_retest"] = {"n_pairs": len(tr), "n_subjects": int(tr.subject.nunique()) if len(tr) else 0,
                         "mean_spearman": float(tr.spearman.mean()) if len(tr) else None,
                         "median_spearman": float(tr.spearman.median()) if len(tr) else None}
    rtt = st["test_retest"]["mean_spearman"]
    st["ceiling_test_retest"] = float(np.sqrt(max(rtt, 0))) if rtt is not None else None

    # ---- 5. verdict (rule fixed above)
    c_hi = st["reliability_replicates"]["ceiling_corr"]
    c_split = float(np.sqrt(max(st["split_replicate"]["mean_spearman"] or 0, 0)))
    st["ceiling_split_replicate"] = c_split
    c_lo = st["ceiling_test_retest"]
    share = st["reliability_replicates"]["signal_share_of_total"]
    central = c_split if st["split_replicate"]["n_samples"] else c_hi
    if central >= PROCEED_CEILING and share >= MIN_SIGNAL_SHARE:
        verdict = "enough signal to proceed"
    elif central >= BORDERLINE_CEILING:
        verdict = "borderline"
    else:
        verdict = "not enough"
    st["verdict"] = verdict
    st["central_ceiling"] = central

    # ---- figures and report
    F = figures(cfg, st, pr, ceil_p, tr, split)
    write(cfg, st, F, sigma2_drug)
    (proc / "exp2_step1_stats.json").write_text(json.dumps(st, indent=2, default=str))
    log.info("verdict: %s", verdict)


def figures(cfg, st, pr, ceil_p, tr, split):
    F = {}
    fig, ax = plt.subplots(1, 4, figsize=(19, 4.2))
    sh = st["shares"]
    ax[0].bar(["drug", "patient", "drug x patient\n(residual)"], [sh["drug"], sh["patient"], sh["residual"]], color=["#555555", "#8172B2", "#DD8452"])
    ax[0].set_ylabel("share of AUC variance"); ax[0].set_title("Additive decomposition")
    ax[1].scatter(pr["first"], pr["last"], s=3, alpha=0.3, c="#4C72B0", rasterized=True)
    ax[1].plot([0, 290], [0, 290], "k--", lw=0.7); ax[1].set_xlabel("AUC, replicate curve 1"); ax[1].set_ylabel("AUC, replicate curve 2")
    ax[1].set_title(f"Technical replicates (r = {st['pearson_replicate_auc']:.3f}, n = {len(pr):,})")
    ax[2].hist(split.values, bins=30, color="#DD8452", alpha=0.8, label="split-replicate")
    if len(tr):
        ax[2].hist(tr.spearman, bins=30, color="#8172B2", alpha=0.6, label="different samples, same patient")
    ax[2].set_xlabel("per-sample Spearman of residuals"); ax[2].set_title("Reproducibility of drug x patient residuals"); ax[2].legend(fontsize=8, frameon=False)
    ax[3].hist(ceil_p.values, bins=30, color="#55A868")
    ax[3].axvline(st["central_ceiling"], color="k", ls="--", lw=1, label="central ceiling")
    ax[3].set_xlabel("sqrt(reliability)"); ax[3].set_title("Per-sample ceiling (replicate noise model)"); ax[3].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    F["main"] = path(cfg, "figures", "exp2_step1_ceiling.png"); fig.savefig(F["main"], dpi=150); plt.close(fig)
    return F


def write(cfg, st, F, sigma2_drug):
    sh, rr, sp, tt = st["shares"], st["reliability_replicates"], st["split_replicate"], st["test_retest"]
    L = []; w = L.append
    w("# Experiment 2, Step 1: ceiling analysis")
    w(f"\nGenerated {dt.date.today().isoformat()} by `python -m pipeline.run --phase exp2-1`. Experiment 2 is a follow-up designed after the "
      "pre-registered Phase 4 analysis found no patient-specific value. Phase 4 results are unchanged.\n")
    w(f"## Verdict: **{st['verdict']}**\n")
    w(f"Central ceiling estimate (best achievable per-sample correlation between a perfect predictor and observed drug x patient residuals): "
      f"**{st['central_ceiling']:.2f}** (bias-corrected {st['split_replicate_debiased']['ceiling']:.2f}). Range across noise estimates: {st['ceiling_test_retest']:.2f} (test-retest, pessimistic) to {rr['ceiling_corr']:.2f} "
      f"(replicate variance model, optimistic). Residual signal is {rr['signal_share_of_total']*100:.1f}% of total AUC variance.")
    w(f"\nRule, fixed before running: proceed if the central ceiling is >= {PROCEED_CEILING} and the residual signal is >= {MIN_SIGNAL_SHARE*100:.0f}% of total variance; "
      f"borderline if the ceiling is {BORDERLINE_CEILING} to {PROCEED_CEILING}; otherwise not enough. The central estimate is the split-replicate one, "
      "because it measures exactly the quantity the metric uses (rank agreement of residuals within a sample) with real repeated measurements.")

    w("\n## 1. Data\n")
    w(f"Same matrix as Phase 4: {st['n_samples']} samples x {st['n_drugs']} drugs, {st['n_obs']:,} observed AUCs (probit-fit AUC, lower = more sensitive). "
      "Why only the residual matters: the Experiment 2 metric is a Spearman correlation across drugs within each sample, which is unchanged by adding a constant "
      "to a sample. So a sample's overall sensitivity (patient effect) cannot help; the usable signal is the drug x patient interaction.")

    w("\n## 2. Variance decomposition\n")
    w("AUC = grand mean + drug effect + patient effect + residual, fitted by alternating least squares on observed entries.\n")
    w(md(pd.DataFrame([{"component": "drug", "share": sh["drug"]}, {"component": "patient (sample)", "share": sh["patient"]},
                       {"component": "residual (drug x patient + noise)", "share": sh["residual"]},
                       {"component": "cross terms (unbalanced design)", "share": sh["cross_terms"]}])))

    w("\n## 3. Measurement noise\n")
    w(f"The open raw dose-response file has {st['n_replicate_pairs']:,} sample-drug pairs measured as two or more separate replicate curves "
      f"({st['n_replicate_samples']} samples). Each replicate curve was refitted with the published probit-AUC definition; our pooled refits reproduce the "
      "published AUCs (Pearson 0.9995 on a 1,500-curve check, median absolute difference 0.0).")
    w(f"\n- Two replicate curves of the same sample and drug agree at Pearson **{st['pearson_replicate_auc']:.3f}** in AUC.")
    w(f"- Noise variance of one curve: {st['sigma2_single_curve']:.1f} AUC units squared (estimated per drug for {st['drugs_with_replicates']} drugs, "
      "pooled for the rest), divided by the number of curves behind each published AUC.")
    w(f"- Residual variance: {rr['var_residual']:.1f}; mean noise variance in the matrix: {rr['mean_noise_var']:.1f}; reliability R = {rr['R']:.3f}; ceiling sqrt(R) = {rr['ceiling_corr']:.3f}.")
    w(f"- Per-sample ceiling under this model: mean {st['per_sample_ceiling']['mean']:.3f} (10th to 90th percentile {st['per_sample_ceiling']['q10']:.3f} to "
      f"{st['per_sample_ceiling']['q90']:.3f}).")
    w("\nCaveat: replicate curves are run on the same day from the same sample preparation, so they miss day-to-day and handling variation. This estimate is optimistic.")

    w("\n## 4. Direct reproducibility of residuals\n")
    w(f"- **Split-replicate:** for {sp['n_samples']} samples with at least 10 drugs measured twice (median {sp['drugs_per_sample_median']} drugs), the residuals from "
      f"curve 1 and curve 2 agree at mean per-sample Spearman **{sp['mean_spearman']:.3f}** (median {sp['median_spearman']:.3f}). A perfect predictor of the true residual "
      f"would correlate with one noisy measurement at about sqrt of this, **{st['ceiling_split_replicate']:.2f}**.")
    w(f"- **Test-retest, different samples of the same patient:** {tt['n_pairs']} sample pairs from {tt['n_subjects']} patients, mean Spearman {tt['mean_spearman']:.3f} "
      f"(median {tt['median_spearman']:.3f}); implied ceiling {st['ceiling_test_retest']:.2f}. These samples were taken at different times (for example diagnosis and relapse), "
      "so real biological change lowers this number; it is a pessimistic bound.")
    sd = st["split_replicate_debiased"]
    w(f"- **Bias check:** the split-replicate residuals above subtract drug and sample effects fitted on AUCs that already include these curves, "
      f"which removes part of the shared signal. Refitting the effects with the replicated pairs masked raises the agreement to {sd['mean_spearman']:.3f}, "
      f"ceiling **{sd['ceiling']:.2f}**. The verdict uses the pre-specified, more conservative {st['ceiling_split_replicate']:.2f}; Step 4 reports "
      "results against both.")
    w(f"- {st['n_samples_resid_below_noise']} samples have a total residual variance below single-curve noise (the spike at 0 in the last panel): their responses "
      "are nearly uniform after removing drug and sample effects, so they carry almost no rankable patient-specific signal.")
    w("\nThe split-replicate drugs are the ones the panel repeats, which may not be representative of all 123 drugs; the variance-model estimate in section 3 covers all drugs.")

    w("\n## 5. What this means for Experiment 2\n")
    p4 = json.loads((path(cfg, "processed", "phase4_stats.json")).read_text())
    m2 = next(r for r in p4["table"] if r["model"] == "M2 combined")["centred_mean"]
    frac = m2 / st["central_ceiling"] if st["central_ceiling"] > 0 else float("nan")
    w(f"- No method can be expected to exceed a per-sample correlation of roughly {st['central_ceiling']:.2f} on observed residuals. For reference, Phase 4's "
      f"combined model reached a drug-centred Spearman of {m2:.3f}, i.e. {frac*100:.0f}% of this ceiling (Phase 4 centred residuals by training-fold drug "
      "mean and SD, close to but not identical to the Experiment 2 target).")
    w("- Results in Step 4 will be reported as a fraction of this ceiling.")

    w("\n## Figure\n")
    w(f"![ceiling](figures/{F['main'].name})")
    top = sigma2_drug.sort_values(ascending=False)
    w("\n## Appendix: single-curve noise variance by drug (largest and smallest)\n")
    w(md(pd.concat([top.head(8), top.tail(8)]).rename("noise_var").to_frame(), index=True))
    path(cfg, "reports", "exp2_step1_ceiling.md").write_text("\n".join(L) + "\n", encoding="utf-8")
