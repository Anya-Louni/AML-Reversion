# Reversion Landscape

An in-silico hypothesis generator for differentiation therapy research in acute myeloid leukaemia (AML), built from public data.

AML cells are stuck in an immature state. Differentiation therapy (ATRA, IDH inhibitors) works by letting malignant cells finish maturing. This project:
- learns the normal myeloid maturation path from healthy bone marrow and places patient cells on it;
- ranks existing compounds by how strongly their expression signatures point toward maturity;
- tests, with pre-registered analyses, whether patient-level information predicts which drugs a patient responds to ex vivo;
- visualizes the idea as a 3D Waddington landscape.

**Research simulation only.** All outputs are computational hypotheses. Nothing here is medical advice or a treatment recommendation.

## Pipeline

```
van Galen 2019 scRNA-seq ──► Phase 1: normal landscape (pseudotime, PHATE map, surface, mapper)
                                   │
                                   ▼
                          Phase 2: place patient cells, deconvolve Beat AML bulk, "stuck" signatures
                                   │
LINCS L1000 (CMap 2020) ──► Phase 3: reversion scores for 1,727 compounds, positive controls
                                   │
Beat AML 2.0 ex vivo  ────► Phase 4: first pre-registered validation (initial analysis)
                                   │
                                   ▼
                          Experiment 2: ceiling analysis, then pre-registered locked-test validation
                                   │
                                   ▼
                          Phase 5: simulation export ──► web app (web/)
```

Phase 0 downloads and audits all sources. Each step writes a report in `reports/`.

## Results summary

**Main validation: Experiment 2 (pre-registered, locked held-out test set).**
- A cell-hierarchy model (cell-type composition x learned per-cell-type drug sensitivity) ranks a patient's drugs better than the global drug ranking on patients from later Beat AML data waves that were never used for fitting.
- Mean per-patient Spearman is **0.20 (95% CI 0.17 to 0.24, p = 0.0001)**, and the top-5 predictions hit the patient's observed top 20% of drugs 40% of the time (chance 20%).
- **The effect is modest:** a ceiling analysis from technical replicates puts the best achievable correlation at about 0.67, so the model reaches about **30% of the ceiling**.
- **Reversion scores add only a small, suggestive gain** (+0.011 to +0.013), which does not survive correction for multiple tests.

**Initial analysis: Phase 4 (first pre-registered test), reported in full.**
- Reversion scores carried weak information about drug response (within-patient Spearman 0.145, p = 0.0001).
- They carried no patient-specific value (p = 0.95).
- This null result motivated the ceiling analysis and Experiment 2.

| Question | Result |
|---|---|
| Does the landscape capture normal maturation? | Pseudotime orders HSC < Prog < GMP < ProMono < Mono at Spearman 0.92; stem markers fall and monocyte markers rise. |
| Are malignant cells placed sensibly? | Mapped scores follow the authors' independent labels at Spearman 0.84. |
| Do known differentiation agents score high? | Among 1,691 compounds in AML/monocytic lines: AM-580 1st, tranylcypromine 2nd, tretinoin 14th (AUROC 0.81, p = 0.003). Calcitriol is a clear miss. |
| How much patient-specific signal does ex vivo data contain? | Drug x patient interaction is 44% of AUC variance; best achievable per-patient correlation about 0.67 (range 0.56 to 0.79). |
| Can drugs be ranked per patient (Experiment 2)? | Yes, modestly: 0.20 on the locked test set, about 30% of the ceiling. A mutation model reached 0.13. |
| Does biology check out? | Learned venetoclax sensitivity is highest in stem/progenitor-like cells, and mature share predicts venetoclax resistance (Spearman -0.35), matching Pei et al. 2020. |

**Related work.** Bottomly et al. 2022 and Zeng et al. 2022 already showed that AML differentiation hierarchy is associated with drug response. Our contribution is a pre-registered, locked-test validation of per-patient drug ranking, plus its integration with the reversion pipeline.

Reports: [phase 0](reports/phase0_data_report.md), [1](reports/phase1_landscape_report.md), [2](reports/phase2_patients_report.md), [3](reports/phase3_drug_ranking_report.md), [4 plan](reports/phase4_plan.md), [4 results](reports/phase4_validation_report.md), [Experiment 2: ceiling](reports/exp2_step1_ceiling.md), [plan](reports/exp2_plan.md), [development](reports/exp2_step3_development.md), [test results](reports/exp2_results.md), [5](reports/phase5_simulation_report.md), [final report](reports/final_report.md).

## The simulation

The 3D landscape animates cells moving downhill when their expression is shifted along a drug's LINCS signature. It is a tool for **hypothesis generation**, not a prediction: a linear expression shift does not model dose, time, toxicity, resistance or cell fate. Use it to explore which drugs point in a maturation-like direction for which patient, then check the evidence in the reports.

## Limitations

- Experiment 2 was designed after the Phase 4 null. Its test patients were already used with different models in Phase 4, so the test set is untouched only with respect to the Experiment 2 models.
- Ex vivo assays measure 3-day viability, not maturation, and the panel has no clinical differentiation agents.
- The hierarchy model captures only about 30% of the achievable patient-specific signal, and its predictions are for ranking, not calibrated response values.
- LINCS signatures come from cell lines after short exposures (mostly 6 to 24 h).
- The web app applies the hierarchy model, trained on Beat AML bulk deconvolution, to single-cell counts; that transfer is untested.

## Reproduce

Requires Python 3.12 and Node 20+.

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt   # .venv/bin/python on Linux/macOS
make all     # phases 0-5, then builds the web app
make test
```

Experiment 2 runs separately:

```bash
make exp2
```

`make exp2` runs the ceiling analysis and the development-set models. The locked test evaluation (`EXP2_CONFIRM_TEST=yes python -m pipeline.run --phase exp2-4`) can run only once per data directory. `make all` downloads about 0.9 GB of raw data and reads about 0.6 GB of LINCS signatures by HTTP range requests. Configuration and the random seed are in `config.yaml`.

### Web app

```bash
cd web
npm install
npm run dev
```

The app is fully static (React, Vite, react-three-fiber) and reads precomputed files in `web/public/data/`. `web/netlify.toml` deploys it to Netlify (base directory `web`, publish directory `web/dist`).

## Data sources and references

- van Galen P. et al. *Cell* 2019. GEO GSE116256.
- Bottomly D. et al. *Cancer Cell* 2022. Beat AML 2.0 open processed data, CC-BY-4.0, https://biodev.github.io/BeatAML2/
- Tyner J.W. et al. *Nature* 2018.
- Subramanian A. et al. *Cell* 2017. CMap LINCS 2020 release.
- Zeng A.G.X. et al. *Nature Medicine* 2022.
- Pei S. et al. *Cancer Discovery* 2020.
- HGNC complete gene set.

## Team

_Team members and roles to be added._
