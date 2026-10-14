# AML differentiation reversion simulator

An in-silico hypothesis generator for differentiation therapy in acute myeloid leukemia (AML).
It learns the normal myeloid differentiation landscape from healthy bone marrow single-cell data,
places AML patients on it, ranks existing drugs by how strongly their expression signatures push
patient cells toward maturity, and validates the rankings against ex vivo drug responses.

**Research simulation only. Nothing here is a treatment recommendation.**

## Status

| Phase | Content | Status |
|---|---|---|
| 0 | Data acquisition and audit | done, see [reports/phase0_data_report.md](reports/phase0_data_report.md) |
| 1 | Normal differentiation landscape | in progress |
| 2-7 | Patient placement, drug ranking, validation, simulation export, web app, packaging | planned |

## Data sources

- van Galen et al. 2019, *Cell*. Single-cell RNA-seq of AML and healthy bone marrow. GEO GSE116256.
- Bottomly et al. 2022, *Cancer Cell* (Beat AML 2.0). Open processed data, CC-BY-4.0: https://biodev.github.io/BeatAML2/
- CMap LINCS 2020 L1000 compound signatures (Subramanian et al. 2017, *Cell*), public S3 release used by clue.io.
- HGNC complete gene set, for gene symbol harmonization.

Raw and processed data are not stored in this repository; the pipeline downloads them.

## Reproduce

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt   # .venv/bin/python on Linux/macOS
make phase0      # or: python -m pipeline.run --phase 0
make test
```

Configuration lives in `config.yaml`; the random seed is fixed there.
