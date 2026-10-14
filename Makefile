PY ?= .venv/Scripts/python.exe

.PHONY: all phase0 phase1 test clean-processed

all: phase0 phase1

phase0:
	$(PY) -m pipeline.run --phase 0

phase1:
	$(PY) -m pipeline.run --phase 1

test:
	$(PY) -m pytest -q tests

clean-processed:
	rm -rf data/processed/*
