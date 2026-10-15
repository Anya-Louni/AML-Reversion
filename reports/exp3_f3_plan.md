# Experiment 3, F3 pre-registered plan: can we predict which cells stay stuck?

Written 2026-10-15 and committed to git **before any F3 model was fitted**. Constraints were established first in [exp3_f3_feasibility.md](exp3_f3_feasibility.md), which should be read alongside this.

**The whole of F3 is exploratory.** The two constraints below mean no confirmatory claim is available, and the report will lead with them rather than with a number.

## Constraint 1: the stated target is not observable

The brief defines a resistant cell as "a cell whose maturity score fails to cross the threshold after treatment". Single-cell RNA-seq destroys the cell it measures, so no cell is seen both before and after treatment. Per-cell persistence is therefore **not measurable** in any open dataset here.

**Substitute target, used throughout and labelled as such:** instead of asking whether a given cell resists, we ask which untreated **cell states** are over-represented after treatment. Formally, we train a classifier to separate treated from untreated malignant cells using only interpretable state features, and read its output for an untreated cell as that cell's **state-enrichment score**. A state enriched after treatment is one whose cells disproportionately persist.

## Constraint 2: treatment and batch are confounded

GSE228322 is one control 10x run and one treated run, so any treated-versus-control difference includes run-to-run batch effects. **A high AUROC there cannot be attributed to treatment.** We mitigate, without solving, by (a) restricting features to interpretable quantities rather than raw genome-wide expression, and (b) resting the interpretation on the *direction of the shift along the maturity axis*, which is a specific prediction, rather than on separability as such. The van Galen arm has the same issue across timepoints, plus real biological change over weeks.

## Data

| Arm | Source | Units | Split |
|---|---|---|---|
| A | GSE228322, PDX AML, control vs VTP50469 at 16 h | 4,522 control and 2,236 treated cells, **1 PDX, 1 drug** | 5-fold over cells; **no held-out sample is possible** |
| B | van Galen, malignant cells, day 0 vs follow-up | **3 patients** (AML328, AML329, AML707B) | leave-one-patient-out, 3 folds |

Arm B is cytotoxic chemotherapy, so it measures **survival**, not failure to differentiate. These are different phenomena and the report will not merge them.

## Features (identical across models, all from the untreated state)

- Phase 1 landscape: maturity score, embedding x and y, and distance to the normal reference (the off-manifold measure from Phase 2).
- Cell type: the authors' malignant like-type for van Galen; for the PDX, the type predicted by our landscape mapping.
- Markers: CD34, PROM1, CD14, LYZ, CSF1R, MPO, ELANE, log-normalized.

## Models and the bar to beat

| id | Features |
|---|---|
| S1 | maturity score alone |
| S2 | cell type alone |
| **F** | **full: landscape + cell type + markers** |

All are gradient-boosted trees with identical settings; only the feature set differs. Metric: **AUROC for separating treated from untreated cells on held-out data**.

**Decision rule, fixed now.** Cell state carries resistance information beyond the simple baselines if **F beats both S1 and S2 by more than 0.02 AUROC** in **both** arms. Anything less is reported as no improvement.

## Interpretive analysis, less vulnerable to batch

For each arm, report the **direction of the maturity shift**: the difference in mean maturity score between treated and untreated cells, and the share of cells below the Phase 2 maturity threshold in each. A differentiation agent should raise maturity; if instead the cells that persist are the less mature ones, that is the resistance signal the project predicts. This is a directional prediction about one axis and is reported as the primary biological readout.

Feature importances from model F will be listed and compared with published resistance biology, with sources cited and anything unverifiable flagged.

## Link to our own simulation

Phase 5 labels cells "resistant" when the simulated shift never carries them past the maturity threshold. For the 3 van Galen patients we compare that **simulated** resistant fraction with the **measured** persistence of the same cell states. This is a weak test and is labelled as such: Phase 5 simulates differentiation agents while these patients received chemotherapy, so a null result would be expected and would not condemn the simulation.

## Dropped from the brief

- **Revumenib responder link:** no open dataset carries responder labels (F0).
- **Beat AML ex vivo link at the patient level:** van Galen and Beat AML are different cohorts, so predicted resistant fractions cannot be matched to ex vivo responses per patient.

## What a null result means

If F does not beat the baselines, the honest conclusion is that cell state beyond maturity and type adds nothing detectable here, on data this thin. That will be reported as the result.
