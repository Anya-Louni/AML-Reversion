PY ?= .venv/Scripts/python.exe

.PHONY: all phase0 phase1 phase2 phase3 phase4 phase5 exp2 exp3 web test clean-processed

all: phase0 phase1 phase2 phase3 phase4 phase5 web

phase0:
	$(PY) -m pipeline.run --phase 0

phase1:
	$(PY) -m pipeline.run --phase 1

phase2:
	$(PY) -m pipeline.run --phase 2

phase3:
	$(PY) -m pipeline.run --phase 3

phase4:
	$(PY) -m pipeline.run --phase 4

phase5:
	$(PY) -m pipeline.run --phase 5

exp2:
	$(PY) -m pipeline.run --phase exp2-1
	$(PY) -m pipeline.run --phase exp2-3

exp3:
	$(PY) -m pipeline.run --phase exp3-f1
	$(PY) -m pipeline.run --phase exp3-f2
	$(PY) -m pipeline.run --phase exp3-f3
	$(PY) -m pipeline.run --phase exp3-f4
	$(PY) -m pipeline.run --phase exp3-export

web:
	cd web && npm ci && npm run build

test:
	$(PY) -m pytest -q tests

clean-processed:
	rm -rf data/processed/*
