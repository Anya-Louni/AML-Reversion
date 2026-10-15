# Reversion Landscape

An in-silico hypothesis generator for differentiation therapy research in acute myeloid leukaemia (AML), built from public data.

AML cells are stuck in an immature state. Differentiation therapy (ATRA, IDH inhibitors) works by letting malignant cells finish maturing. This project learns the normal myeloid maturation path from healthy bone marrow, places patient cells on it, ranks existing compounds by how strongly their expression signatures point toward maturity, visualises the idea as a 3D Waddington landscape, and tests the rankings against ex vivo drug response data.

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
Beat AML 2.0 ex vivo  ────► Phase 4: pre-registered validation vs baselines
                                   │
                                   ▼
                          Phase 5: simulation export ──► Phase 6: web app (web/)
```

Phase 0 downloads and audits all sources. Each phase writes a report in `reports/`.

## Results summary

| Question | Result |
|---|---|
| Does the landscape capture normal maturation? | Pseudotime orders HSC < Prog < GMP < ProMono < Mono at Spearman 0.92 (Palantir 0.90). Stem markers fall and monocyte markers rise. |
| Are malignant cells placed sensibly? | Mapped scores follow the authors' independent malignant cell-type labels at Spearman 0.84. |
| Do known differentiation agents score high? | Among 1,691 compounds profiled in AML/monocytic lines: AM-580 ranks 1st, tranylcypromine 2nd, tretinoin 14th. Tier-1 controls AUROC 0.81 (p = 0.003). Calcitriol is a clear miss. |
| Do reversion scores relate to ex vivo response? | Weakly: mean within-sample Spearman 0.145 (pre-registered P1, p = 0.0001). |
| Do they add patient-specific value? | No, by the pre-registered test (P2: -0.008, p = 0.95). A global drug ranking is far stronger (0.70). In exploratory analyses, per-drug calibration shows a small patient-level signal, and the score adds almost nothing on top of mutations and cell-type fractions. |
| Sanity check | Deconvolved mature share predicts venetoclax resistance (Spearman -0.35, p = 3e-12), matching published biology. |

Full details: [phase 0](reports/phase0_data_report.md), [1](reports/phase1_landscape_report.md), [2](reports/phase2_patients_report.md), [3](reports/phase3_drug_ranking_report.md), [4 plan](reports/phase4_plan.md), [4 results](reports/phase4_validation_report.md), [5](reports/phase5_simulation_report.md), and the [final report](reports/final_report.md).

## Limitations

- LINCS signatures come from cell lines after short exposures (mostly 6 to 24 h); transfer to patient cells is uncertain.
- The animation is a linear expression-shift model, a visualisation of a hypothesis, not a model of cell fate.
- Beat AML ex vivo assays measure 3-day viability, not maturation, and the panel has no clinical differentiation agents.
- Rankings are dominated by drug-level effects; patient-specific differences are small.

## Reproduce

Requires Python 3.12 and Node 20+.

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt   # .venv/bin/python on Linux/macOS
make all     # phases 0-5, then builds the web app
make test
```

`make all` downloads about 0.9 GB of raw data and reads about 0.6 GB of LINCS signatures by HTTP range requests. Configuration and the random seed are in `config.yaml`. Raw and processed data are not stored in the repository.

### Web app

```bash
cd web
npm install
npm run dev
```

The app is fully static (React, Vite, react-three-fiber) and reads precomputed files in `web/public/data/`. `web/netlify.toml` deploys it to Netlify.

## Data sources

- van Galen P. et al. Single-cell RNA-seq reveals AML hierarchies relevant to disease progression and immunity. *Cell* 2019. GEO GSE116256.
- Bottomly D. et al. Integrative analysis of drug response and clinical outcome in acute myeloid leukemia. *Cancer Cell* 2022. Open processed data, CC-BY-4.0, https://biodev.github.io/BeatAML2/
- Tyner J.W. et al. Functional genomic landscape of acute myeloid leukaemia. *Nature* 2018.
- Subramanian A. et al. A next generation Connectivity Map: L1000 platform and the first 1,000,000 profiles. *Cell* 2017. CMap LINCS 2020 release.
- HGNC complete gene set, for gene symbol harmonisation.

## Team

_Team members and roles to be added._
