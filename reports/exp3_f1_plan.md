# Experiment 3, F1 pre-registered plan: does the linear simulation predict real drug-induced change?

Written 2026-10-15 and committed to git **before the analysis was run**. Data availability was established in [exp3_step0_data.md](exp3_step0_data.md). Any deviation will be reported as a deviation.

## Headline limitation, stated in advance

**This is a cell-line-only validation.** Step F0 found no open patient before/after transcriptome dataset treated with a drug that has a signature in our LINCS set. The question of whether the simulation transfers to patient cells therefore **cannot be answered here and will remain untested**. The F1 report must lead with this.

**The main test is deliberately favourable.** Our tretinoin consensus signature is built from HL60, SKM1, THP1 and U937, and **HL-60 is one of those lines**. Both primary datasets are ATRA in HL-60. So D1 and D2 are close to a best case for our pipeline: same drug, same cell line, same direction of perturbation. A good result there is the ceiling of what this approach can do, not evidence that it generalizes. D4 below is the deliberately unfavourable comparison.

## Comparisons

| id | Dataset | Cells | Contrast | Drug in our set | Tier | Role |
|---|---|---|---|---|---|---|
| D1 | GSE227542 | HL-60 | ATRA 4 h, 2 d, 4 d vs 0 h | tretinoin | AML/monocytic, includes HL60 | primary, favourable |
| D2 | GSE276457 | HL-60 | ATRA 48 h vs untreated | tretinoin | AML/monocytic, includes HL60 | primary, favourable |
| D3 | GSE276457 | HL-60 | palbociclib 48 h vs untreated | Palbociclib | AML/monocytic (THP1) | primary, different drug class |
| D4 | GSE228321 | MOLM13 | OTX015 16 h vs control | OTX-015 | non-myeloid only | secondary, unfavourable transfer |

Combination arms are held back for F4 and are not used in F1.

## Quantities

- **Real change**, per gene: counts to CPM, log2(CPM + 1), then mean(treated replicates) minus mean(control replicates). Genes with zero counts in every sample of a contrast are dropped.
- **Simulated change**, per gene: exactly the vector the Phase 5 simulation adds, `normalize(s_g * sd_g) * target_norm`, where `s` is the drug's LINCS consensus z-score, `sd_g` the gene's standard deviation in the normal reference, and the normalization and `target_norm` are the values already fixed in `config.yaml`. Restricted to landmark and best-inferred L1000 genes, as in Phase 5.
- Comparisons are made over genes shared by the dataset and the simulated vector. The shared-gene count is reported for every comparison.
- **Secondary variant:** the same tests using the raw LINCS z-score `s_g` without the `sd_g` weighting, to separate the signature itself from our implementation choice.

## Test 1, direction (confirmatory)

- Metrics: cosine similarity and Spearman correlation between the simulated and real change vectors.
- **Primary statistic:** the rank of the true drug's cosine among **all 1,691 compounds with an AML/monocytic-tier consensus signature**, each scored identically against the same real change vector. Empirical p = rank / 1691, one-sided.
- Baselines, all scored identically:
  1. the full 1,691-drug null distribution (above);
  2. **eight pre-specified negative controls** with no established myeloid differentiation mechanism: ciprofloxacin, famotidine, loratadine, paracetamol, lisinopril, furosemide, glipizide, lidocaine. These were chosen before running on the basis of drug class, and all have AML/monocytic-tier signatures so the comparison is like for like. Compounds with reported differentiation or AML activity (statins, metformin, methylxanthines, allopurinol) were deliberately excluded from this list;
  3. gene-label permutation of the simulated vector, 1,000 draws;
  4. where a dataset has two or more single-agent arms, the mean real change across its other drugs.
- **Decision rule, fixed now.** For each comparison, the simulation "predicts the direction of real change" if **cosine > 0 and the drug ranks in the top 5% of the 1,691-drug null**. Each comparison is judged separately; D1 and D2 are the primary claims.

## Test 2, maturity (confirmatory)

Two scores, both pre-registered, because our landscape mapper was trained on single-cell data and these datasets are bulk:

- **(a) Landscape score:** the Phase 1 mapper applied to the bulk profile. Platform transfer is untested, so this is reported with that caveat.
- **(b) Maturation signature score:** mean standardized expression of the up-with-maturation genes minus the down genes, from `maturation_signature.csv`. Platform-robust and the primary of the two.

For each comparison, compute each score for the real control, the real treated sample, and the simulated profile (control plus the simulated shift). Report:
- **direction:** the sign of the real change and of the simulated change;
- **magnitude:** the ratio simulated change over real change. Magnitude is expected to be poorly calibrated; it is reported, not tested.

**Decision rule, fixed now.** Maturity direction is "correct" for a differentiation agent if the real change and the simulated change are both positive on score (b). For palbociclib and OTX015 no direction is predicted in advance, so these are descriptive.

## Test 3, trajectory (exploratory, D1 only)

Our simulation is a straight line in expression space. Real cells may curve.
- Take real change vectors at 4 h, 2 d and 4 d relative to 0 h.
- Project each onto the unit simulated direction: parallel component `a_t`, residual `r_t - a_t * u`.
- Report the off-axis fraction `||residual|| / ||r_t||` at each timepoint, and whether `a_t` increases monotonically across the three timepoints.
- Descriptive only, no decision rule. Note the real timepoints are 4 h, 2 d and 4 d, not the days 1, 3 and 6 assumed in the brief.

## Dropped from the brief

- **Responder test:** dropped. GSE216730 is PDX with no responder labels; no open dataset supports it.
- **Cell-line versus patient comparison:** cannot be run, see the headline limitation. The report will state this rather than substituting a weaker proxy.

## What a null result means

If the simulation fails these tests, that is a publishable negative result about the method and will be the headline of the F1 report and of the project write-ups. In particular, a failure on D1 or D2, the favourable cases, would mean the linear-shift simulation should not be presented as predictive of real differentiation at all.
