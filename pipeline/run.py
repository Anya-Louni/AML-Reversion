"""Entry point: python -m pipeline.run --phase N"""
from __future__ import annotations

import argparse
import importlib
import random

import numpy as np

from pipeline.common import load_config

PHASES = {0: "pipeline.phase0_data.run", 1: "pipeline.phase1_landscape.run", 2: "pipeline.phase2_patients.run", 3: "pipeline.phase3_drugs.run",
          4: "pipeline.phase4_validation.run", 5: "pipeline.phase5_simulation_export.run",
          "exp2-1": "pipeline.exp2.step1_run", "exp2-3": "pipeline.exp2.step3_run"}
PHASES = {str(k): v for k, v in PHASES.items()}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", type=str, required=True, choices=list(PHASES))
    ap.add_argument("--config", default=None)
    args = ap.parse_args()
    cfg = load_config(args.config)
    random.seed(cfg["seed"])
    np.random.seed(cfg["seed"])
    importlib.import_module(PHASES[args.phase]).run(cfg)


if __name__ == "__main__":
    main()
