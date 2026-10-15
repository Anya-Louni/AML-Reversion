# Experiment 2 pre-registered analysis plan: patient-specific drug matching

Written 2026-10-15 and committed to git together with the split file `reports/exp2_split.csv`, **before any Experiment 2 model was fitted and before the locked test set was evaluated**. The test set is evaluated exactly once (enforced in code by a lock file). Any deviation will be reported as a deviation.

## Why a second experiment

The first pre-registered validation (Phase 4, `reports/phase4_plan.md`) found that reversion scores carry weak information about ex vivo response but no patient-specific value. A ceiling analysis (`reports/exp2_step1_ceiling.md`) then showed that Beat AML ex vivo data does contain patient-specific signal: drug x patient interaction is 44% of AUC variance, residual signal about 27% of total, and the best achievable per-sample correlation is estimated at 0.67 (bias-corrected 0.77; range 0.56 to 0.79). Experiment 2 was designed after seeing both results and is a separate test with its own held-out data.

## Data and split

- Samples: the same 520 Beat AML AML samples with RNA-seq and ex vivo data as Phase 4 (479 patients).
- **Split** (`reports/exp2_split.csv`, rule in `pipeline/exp2/split.py`): the open clinical summary labels each sample with its data-lock cohort. Patients whose samples are all labelled `Waves3+4` form the **locked test set** (169 patients, 178 samples). All other patients (`Waves1+2`, `Both`, or mixed labels) form the **development set** (310 patients, 342 samples). No patient is on both sides. The split is balanced for IDH hotspot mutations (19% vs 19%), FLT3-ITD (25% vs 25%) and specimen type (53% vs 54% bone marrow).
- **Drugs:** the 104 Phase 4 drugs with at least 50 development samples and at least 20 test samples.
- **Prior exposure of the test set.** Phase 4 already used every patient, including the test patients, with the earlier models; the Step 1 ceiling analysis used all samples descriptively (no model fitting, no tuning). The test set is therefore untouched only with respect to the Experiment 2 models. Nothing in Experiment 2 is tuned on test data. The upstream inputs (landscape, maturation signature, deconvolution reference, LINCS signatures, reversion scores) were built from van Galen and LINCS data and never used Beat AML response.

## Target

- Response = sensitivity = minus the probit AUC (higher = more sensitive), as in Phase 4.
- **Residual** for sample p and drug d = response(p, d) minus the mean response of drug d over **development samples only**. The same development drug means are used for test samples.

## Models

All fitting and tuning use development data only, with 5-fold patient-level cross-validation inside the development set (folds assigned by seeded shuffle of patients, seed 0).

| id | role | model |
|---|---|---|
| A | baseline | **Global ranking**: predicted residual = 0 for every pair. |
| B | baseline | **Mutation model**: per-drug ridge regression of the residual on mutation and clinical features: genes mutated in at least 5% of development samples with DNA data, clinical FLT3-ITD and NPM1, a targeted-panel indicator, a missing-DNA indicator, specimen type (bone marrow, blood, leukapheresis) and age at specimen collection. |
| H | **PRIMARY** | **Hierarchy-weighted model**: predicted response(p, d) = sum over cell types k of fraction(p, k) x sensitivity(k, d). Fractions are the 10 deconvolved fractions from Phase 2 (six malignant van Galen like-types plus T/NK, B, erythroid, pDC; they sum to 1). Per drug, the residual is regressed on the fractions with ridge regularization and no separate intercept, so sensitivity(k, d) = development drug mean + coefficient(k, d). |
| HR | secondary | H plus the Phase 3 reversion score of (p, d) as an extra feature. |
| HM | secondary | HR with reversion scores recomputed using only signatures from AML/monocytic lines (THP1, HL60, U937, SKM1, NOMO1); drugs without such signatures are excluded from this model and the number remaining is reported. |

- **Regularization**: for each model, one global ridge penalty, chosen from {0.01, 0.1, 1, 10, 100, 1000} by maximizing the primary metric in development cross-validation, then refitted on the whole development set. Features are standardized on development data.
- **The primary model is H.** Its choice is fixed here and will not change based on development or test results. All other models are secondary.

## Metrics

1. **Primary metric:** for each test sample with at least 20 tested drugs, the Spearman correlation across drugs between predicted and observed residuals. Samples of the same patient are averaged first, so each patient counts once. The metric is the mean over test patients, with a patient-level bootstrap 95% CI (10,000 resamples).
2. Overall Pearson correlation between predicted and observed residuals, pooled over all test pairs.
3. **Top-5 hit rate:** for each test sample, the share of its 5 drugs with the highest predicted residual that fall in the top 20% of its observed residuals (chance = 0.20), averaged over patients.

Model A predicts a constant, so its per-patient Spearman is defined as 0 and its hit rate as chance (0.20).

## Primary comparison and decision rule

- **Primary test:** model H vs model A on the primary metric. Paired sign-flip permutation test over test patients (10,000 permutations) of the per-patient difference H minus A, one-sided (H greater).
- **"Patient-specific value supported"** if and only if the primary test gives p < 0.05 **and** the 95% bootstrap CI of H's primary metric excludes 0.
- **Secondary comparisons** (reported, not used for the decision): H vs B, HR vs H, HM vs H on the same metric and test (two-sided), with the same CI method; within-patient drug-label permutation test for H as a robustness check.
- **Interpretation:** H's primary metric will also be reported as a fraction of the Step 1 ceiling, against both the pre-specified estimate (0.67) and the bias-corrected one (0.77).

## Biological sanity check (reported, not used for the decision)

The learned sensitivities (sensitivity(k, d), drugs x cell types) will be shown as a heatmap. Expected from prior work: venetoclax more effective (higher sensitivity) in stem/progenitor-like cells than in monocyte-like cells (Pei et al. 2020).

## Amendment 1 (2026-10-15, committed before any Experiment 2 model was fitted)

Two secondary analyses are added. Neither changes the primary model, metric, test or decision rule.

1. **Excluding low-signal samples.** All metrics and the H vs A comparison are repeated after removing the 74 samples defined in Step 1 as having total residual variance below the single-curve measurement noise variance (775.3 AUC units squared, residuals from the additive fit on the full Phase 4 matrix). The exact list is fixed in `reports/exp2_low_signal_samples.csv` (59 development, 15 test samples). It is excluded from development fitting and test scoring in this secondary analysis only. Disclosure: the list was derived from response data of all samples, including test samples, in the descriptive Step 1 analysis; it uses no model and no tuning, but it is not independent of test responses, so this analysis is secondary.
2. **Descriptive Phase 4 residual correlation.** The Phase 4 raw reversion score (Phase 3 cosine, model M1 in Phase 4, no fitting) is scored on the test set with the Experiment 2 target and primary metric (per-patient Spearman with development-mean residuals, bootstrap CI). It is reported descriptively, without a hypothesis test, so the earlier score can be read on the same scale as the new models.

## Caveat

Ex vivo response is a 3-day viability readout. Experiment 2 asks whether patient-specific viability response can be predicted; it does not test whether drugs induce differentiation.
