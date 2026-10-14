PY ?= .venv/Scripts/python.exe

.PHONY: all phase0 test clean-processed

all: phase0

phase0:
	$(PY) -m pipeline.run --phase 0

test:
	$(PY) -m pytest -q tests

clean-processed:
	rm -rf data/processed/*
