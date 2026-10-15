# Experiment 3, F2 pre-registered plan: can a learned model beat the linear shift?

Written 2026-10-15 and committed to git **before any F2 model was fitted**. Follows [exp3_f1_results.md](exp3_f1_results.md), which found that the linear shift predicts the direction of real change above chance (cosine 0.28 to 0.62 at 16 h to 4 d) but misses 95% or more of the real movement and points the wrong way at 4 h.

## Headline constraint: the open data cannot support a deep dynamics model

Step F0 found every usable external dataset. Counting the **single-agent conditions with a measured before/after transcriptome**:

| # | Drug | Cell line | Time | Drug has a LINCS signature |
|---|---|---|---|---|
| 1-3 | ATRA | HL-60 | 4 h, 2 d, 4 d | yes |
| 4 | ATRA | HL-60 | 48 h | yes |
| 5 | palbociclib | HL-60 | 48 h | yes |
| 6 | OTX015 | MOLM13 | 16 h | yes (non-myeloid tier) |
| 7-9 | VTP50469 | MOLM13, MV4-11, OCI-AML3 | 16 h | **no** |

That is **9 conditions and 4 distinct drugs, of which 3 have a usable signature**, across 2 timepoint series.

**Pre-registered decision: we will not fit a neural ODE or a flow-matching model.** With three conditionable drugs, a continuous-time neural model in latent space could be made to fit, but it could not be honestly trained or evaluated, and any held-out number it produced would be noise. The brief's own rule is to start simple and stop when added complexity stops helping; here the data runs out before the complex models become meaningful, so we stop in advance rather than report an unsupportable result. If this project later gains a large perturbation dataset, that is the point at which these models become testable.

What we can do honestly is test **two specific, pre-specified hypotheses about why the linear shift underperforms**, each with enough data behind it.

## Hypothesis A: the consensus signature is the wrong context

Our simulation uses a drug's consensus signature averaged over AML/monocytic lines. F1's best case (ATRA in HL-60) still only reached cosine 0.41. If cell context is what is missing, a signature from the matching cell line should do better.

LINCS has an **HL-60-specific tretinoin signature**, and **46 drugs have both an HL-60 signature and at least one other myeloid-line signature**, which is enough to train a regularized map.

## Hypothesis B: the model has no time axis

Our simulation predicts the same change at 4 h, 2 d and 4 d. The real response at 4 h points a different way entirely.

## Models

Evaluated on the external GEO datasets, all with identical metrics.

| id | Model | Fitted on | Note |
|---|---|---|---|
| B0 | no change (zero vector) | nothing | floor |
| B1 | **current linear shift** (consensus signature) | nothing | **the model to beat** |
| B2 | mean real change of the other drugs in the same dataset | external data | only where a dataset has 2 or more single-agent arms |
| M1 | cell-line-matched LINCS signature instead of the consensus | nothing | ATRA/HL-60 only, the one case where a matched signature exists |
| M2 | **learned linear map** from consensus signature to HL-60-specific signature, ridge with reduced rank, trained on the 46 LINCS drugs with both | LINCS only | the test drug is always removed from training; leave-one-drug-out used to pick the penalty and rank |
| M3 | time-interpolation on the ATRA course: fit the real change at 4 h and 4 d, predict 2 d | external data at other timepoints | answers Hypothesis B |

**M3 is not a fair competitor to B1 and will not be reported as one.** It uses real measurements of the same drug at other timepoints, which B1 does not have. It answers a different question: is the real path smooth enough in time that interpolation works, i.e. is a time axis worth adding at all?

## Splits

- **M2 training:** LINCS only, never external data. The evaluated drug is removed from the training set, so the external comparison is a held-out drug in every case. The penalty and rank are chosen by leave-one-drug-out **within the 46 training drugs**, never using external data.
- **M3:** hold out the 2 d timepoint, fit on 4 h and 4 d.
- No external measurement is ever used to fit a model that is then scored on that same measurement.

## Metrics

Per external condition, identical to F1:
1. **cosine** between predicted and real change (primary);
2. per-gene **Spearman**;
3. **maturity shift**: sign agreement of the maturation-signature direction, and its value;
4. rank of the prediction among the 1,691-drug null, as in F1, where the prediction is a drug signature.

## Decision rule, fixed now

A learned model (M1 or M2) **beats the linear shift** if, across the external conditions where it applies:
1. it wins on cosine in a **majority** of conditions, and
2. the **mean cosine improvement over B1 exceeds 0.05**.

Both must hold. With at most six external conditions this is **underpowered and indicative, not confirmatory**; the report must say so and give the per-condition numbers rather than only the summary.

**If the linear shift wins, that is the result.** F1 already showed the linear model has real but limited skill; it is entirely possible that nothing available here improves on it, and that finding will be reported as the headline with no further model search.

## Single-cell distribution metrics

The brief asks for energy distance or MMD between predicted and real treated cell distributions. The only open single-cell before/after AML dataset (GSE228322) is control versus **VTP50469, which has no LINCS signature**, so no model in this project can be conditioned on it and there is nothing to compare a prediction against. We will report the measured distribution shift descriptively and carry those cells into F3, rather than inventing a comparison.

## Reporting

`reports/exp3_f2_results.md` will state, in order: whether anything beat the linear shift, by how much, and the honest size of the evidence. If the gain is small, the report will say the linear model is a legitimate winner, as the brief requires.
