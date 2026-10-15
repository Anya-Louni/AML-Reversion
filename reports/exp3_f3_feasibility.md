# Experiment 3, F3: data feasibility notes (written before the F3 plan)

Verified 2026-10-15 against the actual files. F3 was **not run**; this records what the data can and cannot support, so the plan can be written without re-deriving it.

## Four constraints that shape any honest F3 design

**1. The brief's definition of a resistant cell is not directly observable.** It asks for "a cell whose maturity score fails to cross the threshold after treatment". Single-cell RNA-seq destroys the cell it measures, so the same cell is never seen before and after. Without lineage barcoding, per-cell persistence cannot be measured. Any F3 must either substitute a state-level quantity (which untreated cell states are over-represented after treatment) and label it clearly as a substitute, or say the question cannot be answered.

**2. Treatment and batch are perfectly confounded in the only usable single-cell dataset.** GSE228322 is one control 10x run and one treated 10x run. A classifier separating treated from control cells cannot distinguish treatment biology from run-to-run batch effects. A high AUROC there would not license a biological claim. Restricting features to interpretable ones (landscape position, cell type, maturity score) limits but does not remove this.

**3. There is almost no held-out-sample power.**

| Source | Usable units | Note |
|---|---|---|
| GSE228322 | **1 PDX, 1 drug** (4,522 control + 2,236 treated cells) | no held-out sample possible; random cell splits only |
| van Galen | **3 patients** (AML328, AML329, AML707B) | the only ones with >= 50 malignant cells at both day 0 and a follow-up |

Full van Galen counts (malignant cells, day 0 then follow-ups): AML328 693 then 1,145/153/161; AML329 253 then 259/16; AML707B 1,370 then 61/13/1/5. Every other patient either starts with none or loses nearly all of them (AML556 2,062 to 0; AML475 308 to 1; AML420B 100 to 1). Leave-one-patient-out over three patients is three folds, which is very weak.

**4. The drug in the single-cell data has no LINCS signature.** VTP50469 is absent from our set (confirmed in F0), so no model in this project can be conditioned on it. Resistance can be described for those cells but not predicted from a drug signature.

## Consequences for the brief's F3 steps

- **Step 2 (per-cell classifier):** feasible only against a state-level substitute target, not the stated definition.
- **Step 3 (held-out patients or cell lines):** feasible only in van Galen, with three folds. Report as underpowered.
- **Step 4 (revumenib responder link):** **not possible.** F0 established GSE216730 is PDX with no responder labels, and no open dataset carries them.
- **Step 4 (Beat AML ex vivo link):** not possible per patient. The van Galen and Beat AML cohorts are different people, so predicted resistant fractions cannot be matched to ex vivo responses at the patient level.
- **Step 5 (interpretation against published biology):** feasible, but the van Galen arm measures survival under cytotoxic chemotherapy, not failure to differentiate. Those are different phenomena and must not be conflated.

## Suggested honest scope, if F3 proceeds

Ask whether an untreated cell's state predicts that its state is over-represented after treatment, with maturity score alone and cell type alone as the baselines to beat, evaluated on held-out cells in GSE228322 and leave-one-patient-out in van Galen. Treat the whole of F3 as exploratory, and lead the report with constraints 1 and 2 rather than with an AUROC.
