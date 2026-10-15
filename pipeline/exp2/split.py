"""Experiment 2 split: locked test set = patients seen only in Beat AML waves 3+4.

Every patient is assigned to exactly one side. A patient goes to the test set only if all of their
samples are labelled 'Waves3+4' in the clinical summary; patients labelled 'Waves1+2', 'Both', or with
mixed labels go to development. The rule needs no randomness, so the split is fully reproducible.
"""
from __future__ import annotations

import datetime as dt
from pathlib import Path

import pandas as pd

TEST_COHORT = "Waves3+4"
LOCK_FILE = "exp2_test_evaluated.lock"


def make_split(samples: pd.Index, subject: pd.Series, cohort: pd.Series) -> pd.DataFrame:
    df = pd.DataFrame({"sample": samples, "subject": subject.loc[samples].astype(str).values,
                       "cohort": cohort.reindex(samples).astype(str).values})
    only_test = df.groupby("subject").cohort.agg(lambda c: set(c) == {TEST_COHORT})
    df["split"] = df.subject.map(only_test).map({True: "test", False: "development"})
    return df.sort_values(["split", "subject", "sample"]).reset_index(drop=True)


def check_split(df: pd.DataFrame) -> None:
    sides = df.groupby("subject").split.nunique()
    if (sides > 1).any():
        raise ValueError("a patient appears on both sides of the split")
    if set(df.split) != {"development", "test"}:
        raise ValueError("split must contain development and test")


def assert_test_not_evaluated(proc: Path) -> None:
    if (proc / LOCK_FILE).exists():
        raise RuntimeError(f"the locked test set was already evaluated ({(proc / LOCK_FILE).read_text().strip()}); "
                           "Experiment 2 allows a single test evaluation")


def mark_test_evaluated(proc: Path) -> None:
    (proc / LOCK_FILE).write_text(f"evaluated {dt.datetime.now().isoformat(timespec='seconds')}\n")
