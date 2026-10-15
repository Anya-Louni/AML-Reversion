# Experiment 3, F4 pre-registered plan: can we predict drug combinations?

Written 2026-10-15 and committed to git **before any F4 analysis was run**. Follows [exp3_f2_results.md](exp3_f2_results.md), which found that no learned model beat the linear shift.

## Two constraints that limit what F4 can claim

**1. Additive and sequential composition are the same thing in a linear model.** The brief asks us to contrast additive composition with sequential application through the learned dynamics, which differ only when the dynamics are nonlinear. F2 tested the available learned alternatives and none beat the linear shift, so the model we carry forward is linear, and for it applying drug A then drug B gives exactly the sum of their shifts. **This contrast therefore cannot be tested**, and we will not manufacture a difference by introducing an untested nonlinear model purely to make the comparison look possible.

**2. The model cannot represent synergy at all.** An additive composition predicts that two drugs together do what they do separately. Any genuine synergy or antagonism is outside what it can express. This matters for the positive control below: a good score there would show that two drugs each point toward maturity, not that we predicted their interaction.

## Data, verified

| Source | Combination | Both drugs in our set | Use |
|---|---|---|---|
| GSE276457 | ATRA + palbociclib, HL-60, 48 h | **yes** | the one real validation |
| GSE228321 | VTP50469 + OTX015, MOLM13, 16 h | no (VTP50469 has no signature) | descriptive only, cannot be predicted |
| Beat AML | **none** | n/a | verified: the "Combination" panel contains 3 single agents on 21 samples, and all 63,395 probit fits are `type = single-agent` |

## Test 0, sanity check, runs first

Compute the cosine between the **real** combination change and each **real** single-agent change in GSE276457. If the real combination is nearly identical to one single agent, the combination is dominated by that drug and Test 1 cannot distinguish any composition rule. This is reported before anything else and governs how Test 1 is read.

## Test 1, main (confirmatory)

Predictions of the ATRA + palbociclib change, each scored by cosine against the real combination change over the simulation's genes:

| id | Rule |
|---|---|
| ATRA only | the single-drug shift |
| palbociclib only | the single-drug shift |
| **additive** | normalized sum of the two shifts |
| elementwise max | per gene, the larger-magnitude of the two shifts |

**Decision rule, fixed now:** additive composition is useful if its cosine **beats both single-drug predictions**. We also report the rank of the additive prediction among all 1,691 single AML/monocytic drugs scored against the same real change: if some unrelated single drug predicts the combination better, that is worth knowing.

## Test 2, positive control ranking (weak by construction)

Score every pair of the 1,691 AML/monocytic-tier drugs by the maturation-signature shift of its additive combination, and report where **tretinoin + tranylcypromine** falls. The evidence for that pair is verified: Schenk et al. 2012 (*Nature Medicine* 18:605-611, doi:10.1038/nm.2661) report that LSD1 inhibition, specifically with tranylcypromine, unlocks an ATRA differentiation response in non-APL AML, with the combination superior to either drug alone in vivo.

**This test cannot confirm synergy**, per constraint 2. A high rank means only that both drugs push the same way. It is reported as a consistency check, not as validation of combination prediction.

## Test 3, descriptive

Report the measured VTP50469 + OTX015 behaviour in MOLM13 against its single agents, stating plainly that we cannot predict it because one drug has no signature.

## What we will not claim

If additive composition fails Test 1, or if Test 0 shows the combination is dominated by one drug, we will report that combination prediction is **not validated** and present any ranked combination output as untested exploratory hypotheses, clearly labelled, exactly as the brief requires. One combination in one cell line is in any case a single data point, and no outcome here can establish that the approach works in general.
