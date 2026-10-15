# Reversion Landscape: final report

An in-silico hypothesis generator for differentiation therapy research in AML, built from public data. Course project. All results are computational hypotheses; nothing here is medical advice or a treatment recommendation.

## Summary

- A normal myeloid differentiation landscape learned from single cells places AML patients' malignant cells consistently with independent labels.
- Reversion scores (LINCS signatures pointing toward maturity) recover known differentiation agents among 1,691 compounds, but a first pre-registered validation (Phase 4) found no patient-specific value for ex vivo drug response.
- A ceiling analysis showed that the ex vivo data does contain reproducible patient-specific signal (best achievable per-patient correlation about 0.67).
- **A second pre-registered experiment with a locked held-out test set (Experiment 2) found modest patient-specific value for a cell-hierarchy model: per-patient Spearman 0.20 (95% CI 0.17 to 0.24), about 30% of the ceiling.** Adding reversion scores gave only a small, suggestive gain.
- **Experiment 3 tested the simulation against real treated cells for the first time.** It predicts the *direction* of
  drug-induced change in AML cell lines (beating 1,691 alternatives in 5 of 6 tests) but captures only 8 to 39% of the
  real movement, is wrong at 4 hours, has never been tested in patient cells, and cannot represent resistance or synergy.

## 1. Introduction

Acute myeloid leukaemia (AML) is a cancer of myeloid progenitors that fail to complete maturation. Differentiation therapy, used clinically with all-trans retinoic acid (ATRA) in acute promyelocytic leukaemia and with IDH inhibitors in IDH-mutant AML, works by letting malignant cells resume maturation. We asked whether public data can be combined to (1) describe where each patient's malignant cells sit on the normal maturation path, (2) rank existing compounds by how strongly their expression signatures point toward maturity, and (3) test whether patient-level information predicts which drugs a patient responds to ex vivo.

## 2. Related work

- **Connectivity Map repurposing.** CMap and LINCS L1000 (Lamb et al. 2006; Subramanian et al. 2017) score compounds by how their expression signatures oppose or mimic a disease signature. Our reversion score follows this idea, using a maturation direction learned from single cells instead of a case-control signature.
- **AML hierarchies.** van Galen et al. 2019 classified malignant AML cells into six normal-like types. Their annotations are the backbone of our reference and checks.
- **Hierarchy-based drug response.** Two studies already established that AML differentiation state shapes drug response:
  - Bottomly et al. 2022 (Beat AML 2.0) deconvolved bulk transcriptomes with van Galen cell types and showed that ex vivo drug sensitivity is broadly governed by differentiation state.
  - Zeng et al. 2022 deconvolved more than 1,000 bulk AML transcriptomes into leukaemia stem, progenitor and mature components. Hierarchy composition was associated with chemotherapy response and with targeted-drug sensitivity profiles, and a seven-gene Primitive-vs-Mature biomarker was associated with response to 105 investigational drugs.
  - Our hierarchy model is built on the same idea; we do not claim it as new.
- **What we add.** (1) A pre-registered test of per-patient drug ranking on a locked held-out test set (patients from later Beat AML data waves), with a ceiling analysis from technical replicates that says how much patient-specific signal exists. (2) Integration with the reversion pipeline, testing whether drug-signature reversion scores add to the hierarchy model. (3) A transparent pipeline and an interactive visualization.
- **Perturbation modelling.** Methods such as scGen or CellOracle model expression responses to perturbations. Our simulation uses a far simpler linear shift, chosen for transparency.

## 3. Methods

Full methods and every parameter choice, with justifications, are in the phase reports.

- **Data (Phase 0).** van Galen GSE116256 (41,090 cells, 43 samples); Beat AML 2.0 open processed data (707 RNA-seq samples, 166 drugs, mutations, raw dose-response); CMap LINCS 2020 Level 5 signatures read by HTTP range requests (12,414 signatures). All genes harmonized to HGNC symbols. 140 of 166 Beat AML drugs matched to LINCS.
- **Landscape (Phase 1).** 7,956 normal myeloid cells, Harmony by sample, diffusion pseudotime checked against Palantir and 10 alternative roots, PHATE embedding, a potential surface, and a fixed PCA + kNN map for new cells. Maturation signature: the top 200 genes up and down with pseudotime, excluding ribosomal and mitochondrial genes.
- **Patients (Phase 2).** Maturity threshold fixed on normal cells. Malignant cells mapped with an off-manifold flag. Beat AML bulk samples deconvolved by NNLS into 10 cell types, validated on pseudo-bulk samples with the patient held out.
- **Drugs (Phase 3).** Consensus signatures per compound; reversion score = cosine between a drug signature and the patient's push-to-maturity vector, with KS connectivity as a second method. Positive controls fixed before scoring.
- **Phase 4 (initial validation).** Pre-registered plan committed before the analysis; 520 samples, 479 patients, 123 drugs; patient-level 5-fold CV against random, global-ranking, mutation and cell-fraction baselines.
- **Experiment 2 (main validation).**
  - **Ceiling:** additive decomposition of AUC and noise from 8,496 technical replicate pairs. Curves were refitted from the raw data; our refits match the published AUCs at Pearson 0.9995.
  - **Pre-registration:** plan and split committed before any model was fitted, plus one amendment committed before fitting.
  - **Split:** development = patients with any sample from waves 1+2 (310 patients); locked test = patients seen only in waves 3+4 (169 patients).
  - **Target:** sensitivity minus the drug's development mean.
  - **Primary model (fixed in the plan):** predicted response = sum over 10 cell types of fraction x learned sensitivity, fitted by per-drug ridge regression on development data.
  - **Primary metric:** mean per-patient Spearman over test samples with at least 20 drugs (this coverage rule was pre-registered), with a patient-level bootstrap CI.
  - **Primary test:** one-sided paired sign-flip permutation against the global ranking. The test set was evaluated once, enforced by a lock file.
- **Simulation (Phase 5).** Expression shifted linearly along each drug's signature, rescaled to a common length so that on average at most 20% of cells leave the normal manifold, and each frame re-mapped.

## 4. Results

**Landscape and placement.**
- Pseudotime orders the normal stages at Spearman 0.920 (Palantir 0.903; root sensitivity minimum 0.990).
- Held-out mapping error is 0.028 on a 0-1 scale.
- Malignant cells' mapped scores follow the authors' labels at Spearman 0.84.
- The maturity threshold separates normal immature from mature cells with 98.7% balanced accuracy.
- Deconvolved mature share is recovered at Pearson 0.895 on held-out pseudo-bulk samples.

**Drug ranking.** Among 1,691 compounds in AML/monocytic lines:
- Tier-1 positive controls rank high: AM-580 1st, tretinoin 14th, AUROC 0.81 (p = 0.003). The LSD1 inhibitor tranylcypromine ranks 2nd.
- Calcitriol scores low, a clear miss.
- Reversion rankings are very similar across patients (between-patient Spearman 0.93 in van Galen).

**Initial analysis: Phase 4 (first pre-registered test), reported in full.**
- Reversion scores carried weak information about which drugs a sample responds to (P1: mean within-sample Spearman 0.145, p = 0.0001).
- They carried no patient-specific signal once each drug's average potency was removed (P2: -0.008, p = 0.95).
- The global drug ranking reached 0.696. A combined model (0.720) beat it, but almost entirely through mutations and cell-type fractions: removing the reversion score changed it by +0.0005 raw and +0.002 drug-centred.
- Exploratory: per-drug calibration of the reversion score gave a drug-centred Spearman of 0.094.
- Mature share correlated with venetoclax resistance (Spearman -0.35), matching published biology.
- IDH-mutant samples showed no extra ex vivo sensitivity to tool IDH inhibitors.
- Pre-registered verdict: information about response yes, patient-specific value no.

**Ceiling analysis.**
- Drug, patient and drug x patient components explain 42%, 13.5% and 44% of AUC variance.
- Two independent replicate curves of the same sample and drug agree at Pearson 0.78.
- The drug x patient residual is reproducible enough that a perfect model could reach a per-patient correlation of about 0.67 (pre-specified estimate; 0.77 after a bias correction; range 0.56 to 0.79 across noise estimates).

**Experiment 2 (main validation), locked test set** (169 patients; 151 had at least 20 of the 104 drugs and were scored):

| Model | Per-patient Spearman (95% CI) | Top-5 hit rate (chance 0.20) |
|---|---|---|
| Global ranking | 0 | 0.20 |
| Mutation + clinical | 0.129 (0.095 to 0.165) | 0.33 |
| **Hierarchy (primary)** | **0.202 (0.165 to 0.238)** | **0.40** |
| Hierarchy + reversion | 0.213 (0.177 to 0.249) | 0.42 |
| Hierarchy + myeloid-line reversion | 0.226 (0.188 to 0.263) | 0.41 |

- **Pre-registered result:** the primary model beat the global ranking (p = 0.0001; CI excludes 0), so the pre-registered criterion for patient-specific value is met.
- **Size:** the effect is modest, about 30% of the pre-specified ceiling (26% of the bias-corrected one).
- **Mutation model:** the hierarchy model beat it on the test set (+0.073, p = 0.004), and held up from development (0.192) to test (0.202) while the mutation model dropped (0.195 to 0.129).
- **Reversion scores:** adding them raised the score by 0.011 (p = 0.03), or 0.013 for the myeloid-line version on its drug subset (p = 0.02). Across three secondary tests, neither survives a Bonferroni correction, so this gain is only suggestive. On its own, the raw reversion score still showed no patient-specific signal (0.002).
- **Robustness:** excluding 74 pre-defined low-signal samples did not change the result (0.209).
- **Biology check:** in the learned cell-type sensitivities, venetoclax has the largest stem/progenitor-over-monocyte preference of all 104 drugs, as expected from Pei et al. 2020; MEK inhibitors lean toward monocyte-like cells.

**Simulation (hypothesis generation).**
- With common-length shifts, 16% of cells end off the normal manifold on average.
- Outcomes vary widely by patient; tretinoin moves a mean of 17% of stuck cells past the threshold.
- These numbers describe the visualization, not predicted clinical effects.

## 5. Experiment 3: is the simulation real?

Everything above treats the 3D simulation as a way of generating hypotheses. Experiment 3 asks whether it corresponds
to anything measured. It was run in four parts, each pre-registered before its analysis, against open transcriptomes of
AML cells actually treated with these drugs ([data audit](exp3_step0_data.md)).

**The constraint that shapes all of it.** No public dataset treats *patient* AML cells with a drug whose signature we
hold. Every result below is a cell line or a xenograft, and the transfer to patients remains untested.

### F1: does the predicted direction match reality? ([plan](exp3_f1_plan.md), [results](exp3_f1_results.md))

Six comparisons against measured responses. In **5 of 6** the predicted direction beat almost all of 1,691 alternative
compounds scored identically (ranks 1, 4, 4, 11 and 37). The eight pre-specified negative controls sat near zero.

But the straight line captures only **8 to 39%** of the variance in what the cells really do, and **95 to 99.9%** of the
ATRA response is perpendicular to our axis. At 4 hours the prediction points the wrong way entirely (cosine -0.046,
rank 1338) and the real maturation shift is negative while we predict positive: the early response is not a small
version of the late one. Our model has no time axis, so it predicts the same change at 4 hours and 4 days.

The maturity score behaves correctly where it matters: for ATRA at 2 d, 4 d and 48 h, real and simulated maturation
shifts are both positive, and the landscape score rises in the real data.

**Claim constrained:** "predicts the direction of drug-induced expression change in AML cell lines" is supported.
"Simulates cell reversion" is not.

### F2: can a learned model do better? ([plan](exp3_f2_plan.md), [results](exp3_f2_results.md))

Using the drug's signature from the **matching cell line** beat the consensus on all 4 testable conditions (+0.079
cosine). All four are the same drug, so this is a lead rather than an established gain, but it points at a concrete fix.
A map learned from 46 LINCS drugs did **not** generalize (mean -0.037, badly wrong on palbociclib). The linear model stands.

No neural ODE or flow-matching model was fitted, and that was pre-registered: three conditionable drugs and nine
external conditions cannot support one, and any held-out number would be noise.

**The clearest finding is about time.** Interpolating between the real 4 h and 4 d responses predicts the held-out 2 d
response at cosine **0.93 against 0.28** for our static prediction, on identical genes. The real path is smooth in
time; our model simply has no time axis. (This is not a fair contest, since interpolation sees real measurements of
the same drug; it measures what a time axis would be worth.)

### F3: can we predict which cells resist? ([plan](exp3_f3_plan.md), [results](exp3_f3_results.md))

**The question cannot be asked properly of open data.** Single-cell sequencing destroys the cell it measures, so no
cell is observed before and after treatment; per-cell resistance is not measurable. We substituted a state-level
question and labelled it as such.

The pre-registered bar was **not met**. Cell state beat the baselines in the xenograft arm (AUROC 0.76) but that arm
confounds treatment with sequencing run, so it cannot be read as biology. Across three patients, every model scored
*below* chance when pooled, meaning a rule learned on two patients actively misleads on the third: nothing transfers.

One result survives the confound: under a menin inhibitor, mean maturity rises **+0.051**, the direction a
differentiation agent should produce on an axis fixed in Phase 1. Under chemotherapy, maturity does not move
(-0.015): surviving malignant cells are not more mature, just fewer.

### F4: can we predict combinations? ([plan](exp3_f4_plan.md), [results](exp3_f4_results.md))

Two structural limits first. With a linear model, applying drug A then B is **identical** to summing their shifts, so
the additive-versus-sequential contrast has no content here. And an additive model **cannot express synergy**, which is
what makes combinations interesting.

One real combination exists in open data (ATRA + palbociclib, HL-60). Predicting from both drugs beat either alone
(0.493 against 0.414 and 0.313). Two caveats: an elementwise-max rule scored marginally higher (indistinguishable on
one combination), and the additive prediction ranks only **3rd** among 1,691 single drugs, so an unrelated compound
predicts this combination better than we do.

The published ATRA + LSD1-inhibitor pair ranks 102 of 1,428,895 pairs, but this does not validate combination
prediction: our score is a normalized sum, so pairs rank high exactly when both drugs already point toward maturity,
and tranylcypromine was already the second-highest single compound in Phase 3.

Beat AML contributes nothing: its "Combination" panel holds 3 single agents, and all 63,395 fits are single-agent.

### What Experiment 3 changes

The simulation is now labelled in the web app with what it has been validated to do. The honest summary is that it
predicts a *direction* in cell lines, accounts for a minority of real movement, has no time axis, has never been tested
in patients, and cannot represent resistance or synergy.

## 6. Limitations

- Experiment 2 was designed after the Phase 4 null. Its test patients were already used, with different models, in Phase 4 and the descriptive ceiling analysis. The test set is untouched only with respect to the Experiment 2 models.
- Ex vivo response is 3-day viability, not maturation; the Beat AML panel has no clinical differentiation agents. Neither validation tests whether drugs induce differentiation.
- The hierarchy model explains only about 30% of the achievable patient-specific signal. Its predictions are strongly shrunk (heavy ridge penalty), so they rank drugs but are not calibrated AUC values.
- LINCS signatures come from cell lines after short exposures; most AML/monocytic signatures are from THP-1 alone.
- The web app applies the hierarchy model, trained on bulk-deconvolved Beat AML fractions, to single-cell counts from van Galen patients. This domain shift is untested.
- The linear expression shift does not model dose, time, toxicity, resistance or cell fate dynamics.
- Only 12 single-cell patients are available, some with few malignant cells.
- **The simulation has no time axis**, which Experiment 3 identified as its largest single gap: knowing the time course
  raises prediction of a held-out timepoint from cosine 0.28 to 0.93.
- **Resistance and synergy are outside the model.** F3 could not validate resistant-state prediction at all, and an
  additive combination model cannot express interaction by construction.
- Every Experiment 3 result is a cell line or xenograft. No open dataset allows the patient test.

## 7. Future work

- A perturbation dataset large enough to fit a time-conditioned dynamics model: Experiment 3 showed the gain is there
  but that three conditionable drugs cannot support the model.
- Lineage-barcoded single-cell experiments, which would make per-cell resistance measurable for the first time.
- Wet-lab testing of top-ranked compounds with differentiation readouts (CD11b, CD14, morphology) rather than viability.
- Hierarchy models that share strength across drugs or use drug target information.
- Perturbation response models trained on myeloid cells and longer exposures.
- Ex vivo datasets that include differentiation agents and maturation endpoints.

## References

- Bottomly D. et al. Integrative analysis of drug response and clinical outcome in acute myeloid leukemia. *Cancer Cell* 2022;40:850-864. doi:10.1016/j.ccell.2022.07.002
- Lamb J. et al. The Connectivity Map. *Science* 2006.
- Pei S. et al. Monocytic subclones confer resistance to venetoclax-based therapy in patients with acute myeloid leukemia. *Cancer Discovery* 2020;10:536-551. doi:10.1158/2159-8290.CD-19-0710
- Subramanian A. et al. A next generation Connectivity Map: L1000 platform and the first 1,000,000 profiles. *Cell* 2017.
- Tyner J.W. et al. Functional genomic landscape of acute myeloid leukaemia. *Nature* 2018.
- van Galen P. et al. Single-cell RNA-seq reveals AML hierarchies relevant to disease progression and immunity. *Cell* 2019.
- Zeng A.G.X. et al. A cellular hierarchy framework for understanding heterogeneity and predicting drug response in acute myeloid leukemia. *Nature Medicine* 2022;28:1212-1223. doi:10.1038/s41591-022-01819-x
