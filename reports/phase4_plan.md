# Phase 4 pre-registered analysis plan

Written 2026-10-14, before any validation analysis was run and before any reversion score was compared with drug response.
This file is committed to git on its own, so its timestamp precedes the validation code and results.
Any deviation from this plan will be listed in the Phase 4 report as a deviation.

## Question

Do reversion scores (how strongly a drug's LINCS signature points along a patient's push-to-maturity vector)
carry information about that patient's ex vivo drug sensitivity in Beat AML, and in particular
patient-specific information beyond knowing which drugs are generally potent?

## Data

- **Samples:** Beat AML 2.0 AML RNA-seq samples (healthy controls excluded) that have ex vivo drug data.
  All such samples are used; several patients contribute more than one sample, so every split and
  every bootstrap is done at the patient level (`dbgap_subject_id`).
- **Drugs:** Beat AML panel drugs that are matched to a LINCS compound, have a consensus signature, and have
  at least 50 RNA-seq samples with data (expected about 120 drugs; the exact number will be reported).
  If fewer than 15 drugs qualify, power is too low and the report says so.
- **Response:** sensitivity = minus the probit-fit AUC (higher = more sensitive).
- **Reversion score:** Phase 3 cosine score (primary) computed from the Phase 2 bulk stuck signature.
- **Features for baselines:**
  - Mutations: genes mutated in at least 5% of samples with DNA data, plus clinical FLT3-ITD and NPM1 calls
    and an indicator for the capture assay (exome vs targeted panel). Samples without DNA data get zeros
    plus a missing-DNA indicator.
  - Differentiation fractions: the six deconvolved malignant shares from Phase 2 plus the bulk maturity index.

## Split

5-fold GroupKFold by patient, fixed seed 0. Every model that is fitted uses training folds only; every metric
is computed on held-out folds. Drug-level centring statistics (drug means and SDs used in P2) are computed
on training folds and applied to test folds. The reversion score itself uses no response data.

## Models (all evaluated identically)

| id | model | description |
|---|---|---|
| B0 | random | random scores; metrics averaged over 100 draws |
| B1 | global drug ranking | mean training-fold sensitivity of each drug, same for every patient |
| B2 | mutation-only | per-drug ridge regression (alpha by inner CV) from mutation features to sensitivity |
| B3 | differentiation-fraction-only | per-drug ridge from the seven fraction features |
| M1 | ours, raw | the reversion score itself, no fitting |
| M1c | ours, calibrated | per-drug linear regression of sensitivity on the reversion score |
| M2 | combined | per-drug ridge on reversion score + fractions + mutations |

All fitted models predict sensitivity for each (patient, drug) pair. Per-drug models include the drug mean
through an intercept.

## Metrics

1. **Per-patient ranking (raw):** Spearman correlation, across the drugs a sample was tested on, between the
   model score and observed sensitivity. Samples with fewer than 20 tested drugs are excluded. Summarized as
   the mean over samples with a patient-level bootstrap 95% CI (2,000 resamples).
2. **Per-patient ranking (patient-specific):** the same, after centring both the score and the sensitivity of
   each drug (subtracting the training-fold drug mean, dividing by the training-fold SD). This removes the
   "which drugs are generally potent" signal; B1 is exactly 0 here by construction.
3. **Per-drug:** for each drug, Spearman across held-out samples between score and sensitivity; summarized as
   the mean over drugs with bootstrap CI over drugs, plus the number of drugs with BH-adjusted p < 0.05.

## Primary tests

- **P1 (as specified in the project brief):** M1, metric 1. Test: permutation of drug labels within each sample
  (10,000 permutations), statistic = mean per-sample Spearman, one-sided (greater).
- **P2 (patient-specific value):** M1, metric 2, same permutation scheme.
- Significance level 0.025 for each (Bonferroni over the two primary tests).

## Claim rules

- "Reversion scores carry information about drug response" requires P1 significant.
- "Patient-specific value" requires P2 significant **and** M1 or M2 beating B1 on metric 1 with a paired
  Wilcoxon signed-rank test over samples (p < 0.05), or M2 beating B1 on metric 1 and B2/B3 on metric 2.
- If neither holds, the result is reported as null. A null here is plausible and informative, because the
  ex vivo assay measures 3-day viability, not differentiation.

## Secondary analyses (reported, not used for claims)

- Same metrics with the connectivity score, landmark-only genes, the full signature (with ribosomal genes),
  and only drugs whose signature comes from AML/monocytic lines.
- **IDH subgroup check:** AGI-5198 sensitivity in IDH1 R132-mutant vs IDH1 wild-type samples, and AGI-6780 in IDH2
  R140/R172-mutant vs wild-type, one-sided Mann-Whitney (mutant more sensitive). Also whether reversion scores
  for these drugs differ between the groups.
- **Biological sanity check of the deconvolution:** monocytic AML is known to be less sensitive to venetoclax
  (Pei et al. 2020, Cancer Discovery). Expected: negative Spearman between mature share and venetoclax sensitivity.

## Caveat to state in the report

Beat AML ex vivo assays measure cell viability after about 3 days. Differentiation agents act by maturation,
often over a longer time. Ex vivo sensitivity is therefore a proxy for our target, and mismatch is expected,
especially for slow-acting differentiation drugs. Also, the panel contains no clinical differentiation agents
(no ATRA, enasidenib or ivosidenib).
