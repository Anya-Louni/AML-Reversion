"""Export the Experiment 3 validation results for the web app."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from pipeline.common import ROOT, get_logger

log = get_logger("exp3.export")
WEB = ROOT / "web" / "public" / "data"


def build(proc: Path) -> dict:
    f1 = json.loads((proc / "exp3_f1_stats.json").read_text())
    f2 = json.loads((proc / "exp3_f2_stats.json").read_text())
    f3 = json.loads((proc / "exp3_f3_stats.json").read_text())
    f4 = json.loads((proc / "exp3_f4_stats.json").read_text())

    direction = [{"id": r["id"], "drug": r["drug"], "role": r["role"], "cosine": r["cosine"],
                  "rank": r["rank_in_null"], "n_null": r["n_null"], "passes": r["passes"],
                  "variance_explained": r["cosine"] ** 2}
                 for r in f1["direction"]]
    f2res = pd.DataFrame(f2["results"])
    models = f2res.pivot(index="id", columns="model", values="cosine").reset_index().to_dict(orient="records")
    return {
        "f1": {
            "direction": direction,
            "n_pass": int(sum(d["passes"] for d in direction)),
            "n_total": len(direction),
            "trajectory": f1["trajectory"],
            "maturity": [{k: r[k] for k in ("id", "drug", "real_shift", "sim_shift", "same_direction")}
                         for r in f1["maturity"]],
            "n_sim_genes": f1["n_sim_genes"],
            "n_null_drugs": f1["n_null_drugs"],
        },
        "f2": {"models": models, "verdict": f2["verdict"], "m3": f2["m3"],
               "n_map_training_drugs": f2["n_map_training_drugs"]},
        "f3": {"verdict": f3["verdict"], "direction": f3["direction"], "overall": f3["overall"],
               "cells": f3["cells"]},
        "f4": {"test0": f4["test0"], "test1": f4["test1"],
               "test2": {k: v for k, v in f4["test2"].items() if k != "top10"},
               "top_pairs": f4["test2"]["top10"], "test3": f4["test3"]},
        "validation_label": validation_label(f1, f2),
    }


def validation_label(f1: dict, f2: dict) -> dict:
    """One honest sentence about what the simulation has been validated to do."""
    d = [r for r in f1["direction"]]
    n_pass = sum(r["passes"] for r in d)
    ve = [r["cosine"] ** 2 for r in d if r["passes"]]
    return {
        "headline": "Direction validated in cell lines; trajectory not validated",
        "detail": (f"Tested against measured drug responses in AML cell lines: the predicted direction of change beat "
                   f"{f1['n_null_drugs']:,} alternative drugs in {n_pass} of {len(d)} comparisons, but captures only "
                   f"{min(ve) * 100:.0f} to {max(ve) * 100:.0f}% of the real movement and is wrong at 4 hours. "
                   "Never tested in patient cells; no open data allows it."),
        "level": "partial",
    }


def run(cfg: dict) -> None:
    proc = ROOT / cfg["paths"]["processed"]
    WEB.mkdir(parents=True, exist_ok=True)
    from pipeline.phase5_simulation_export.run import clean_nan
    payload = clean_nan(build(proc))
    (WEB / "exp3.json").write_text(json.dumps(payload, indent=1, allow_nan=False))
    log.info("wrote exp3.json (%.1f kB)", (WEB / "exp3.json").stat().st_size / 1e3)
