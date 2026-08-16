PY      ?= python3
PYTHONPATH := src
export PYTHONPATH

CC      ?= cc
CFLAGS  ?= -O3 -std=c11 -Wall -Wextra

NPROC   ?=
NPROCFLAG := $(if $(NPROC),--nproc $(NPROC),)

ENGINES := build/hsmc build/hsmc2d build/edmd

.PHONY: all help engines test test-fast \
        gate-g0 gate-g0b gate-g1 gate-g1-2d gate-g2 gate-g3 gate-g4 gate-g5 gate-t8 \
        gates-light gates-heavy geometry degeneracy validate-pilot \
        sweep-2d-smoke campaign-2d campaign-3d isoconfig-smoke paper clean clean-data

help:
	@echo "Light (laptop, seconds-to-minutes; all executed in this repository):"
	@echo "  make engines        compile hsmc, hsmc2d, edmd"
	@echo "  make test-fast      test suite without the slow scans (~40 s)"
	@echo "  make test           full test suite incl. certificates (~1 min)"
	@echo "  make geometry       gates G0 + G0b (theory only)"
	@echo "  make degeneracy     T4 theory run -> results/degeneracy.md"
	@echo "  make validate-pilot analysis ports vs the recorded pilot numbers"
	@echo "  make gate-g1-2d     2D EOS + phi=0.80 audits (T5 gate)"
	@echo "  make gate-g2        Voronoi pipeline on perfect lattices"
	@echo "  make gate-g5        refscore pre-registration proof (git history)"
	@echo "  make gate-t8        topology classifier gate"
	@echo "  make gate-g1        3D EOS vs Carnahan-Starling (~3 min on 9 cores)"
	@echo "  make gates-light    all of the above, blocking"
	@echo ""
	@echo "Heavy (implemented, NOT executed here -- hours to cluster-scale):"
	@echo "  make gate-g3        equilibration-length gate (1e5-sweep points)"
	@echo "  make gate-g4        finite-size gate (N up to 10976)"
	@echo "  make campaign-2d    T6 replication (~355 core-hours) -- the decision node"
	@echo "  make campaign-3d    T9 (~4800 core-hours) -- BLOCKED until T6 passes"
	@echo "  (use scripts/run_sweep.py --preset ... --dry-run for the cost first)"
	@echo ""
	@echo "  make paper          build the manuscript (needs a TeX installation;"
	@echo "                      falls back to the article-class preview if revtex4-2"
	@echo "                      is missing; paper/paper_preview.pdf is checked in)"

engines: $(ENGINES)

build/hsmc: src/hsga/engine/hsmc.c
	@mkdir -p build
	$(CC) $(CFLAGS) -o $@ $< -lm

build/hsmc2d: src/hsga/engine/hsmc2d.c
	@mkdir -p build
	$(CC) $(CFLAGS) -o $@ $< -lm

build/edmd: src/hsga/engine/edmd.c
	@mkdir -p build
	$(CC) $(CFLAGS) -o $@ $< -lm

test:
	$(PY) -m pytest -q

test-fast:
	$(PY) -m pytest -q -m "not slow"

geometry gate-g0 gate-g0b:
	$(PY) scripts/run_geometry.py

degeneracy:
	$(PY) scripts/run_degeneracy.py

validate-pilot:
	$(PY) scripts/validate_pilot.py

gate-g1: engines
	$(PY) -m hsga.gates.gate_s1_eos $(NPROCFLAG)

gate-g1-2d: engines
	$(PY) -m hsga.gates.gate_g1_2d $(NPROCFLAG)

gate-g2:
	$(PY) -m hsga.gates.gate_s2_lattices

gate-g3: engines
	$(PY) -m hsga.gates.gate_s3_equilibration $(NPROCFLAG)

gate-g4: engines
	$(PY) -m hsga.gates.gate_s4_finitesize $(NPROCFLAG)

gate-g5:
	$(PY) -m hsga.gates.gate_g5_preregistration

gate-t8:
	$(PY) -m hsga.gates.gate_t8_topology

# Blocking chain of everything executable on a laptop.
gates-light: geometry gate-g2 gate-g5 gate-t8 gate-g1-2d gate-g1
	@echo "all light gates passed"

# The heavy chain; run only with the compute to back it.
gates-heavy: gate-g3 gate-g4
	@echo "heavy gates passed"

sweep-2d-smoke: engines
	$(PY) scripts/run_sweep.py --preset 2d-smoke $(NPROCFLAG)
	$(PY) scripts/run_analysis.py --data-dir data/2d-smoke --dim 2 --structural $(NPROCFLAG)

campaign-2d: engines
	$(PY) scripts/run_sweep.py --preset 2d-replication $(NPROCFLAG)
	$(PY) scripts/run_analysis.py --data-dir data/2d-replication --dim 2 --structural $(NPROCFLAG)

campaign-3d: engines
	@echo "T9 is gated on the T6 decision (see results/campaign_2d-replication/report.md)."
	@echo "Run scripts/run_sweep.py --preset 3d-campaign explicitly once T6 has passed."

isoconfig-smoke: engines
	$(PY) scripts/run_isoconfig.py --smoke $(NPROCFLAG)

paper:
	@cd paper && ( \
	  if command -v pdflatex >/dev/null 2>&1; then \
	    if kpsewhich revtex4-2.cls >/dev/null 2>&1; then \
	      pdflatex -interaction=nonstopmode paper.tex && pdflatex -interaction=nonstopmode paper.tex; \
	    else \
	      echo "revtex4-2 missing -> building the article-class preview (landmine 8)"; \
	      pdflatex -interaction=nonstopmode paper_preview.tex && pdflatex -interaction=nonstopmode paper_preview.tex; \
	    fi \
	  else \
	    echo "no TeX installation found; paper/paper_preview.pdf is the checked-in build"; \
	    exit 1; \
	  fi )

clean:
	rm -rf build .pytest_cache src/**/__pycache__ tests/__pycache__

clean-data:
	rm -rf data
