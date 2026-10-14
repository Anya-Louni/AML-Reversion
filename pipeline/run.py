"""Entry point: python -m pipeline.run --phase N"""
from __future__ import annotations

import argparse
import importlib
import random

import numpy as np

from pipeline.common import load_config

PHASES = {0: "pipeline.phase0_data.run"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", type=int, required=True, choices=sorted(PHASES))
    ap.add_argument("--config", default=None)
    args = ap.parse_args()
    cfg = load_config(args.config)
    random.seed(cfg["seed"])
    np.random.seed(cfg["seed"])
    importlib.import_module(PHASES[args.phase]).run(cfg)


if __name__ == "__main__":
    main()
