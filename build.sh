#!/usr/bin/env bash
# Rerun every certificate script (optimised Python: checks never rely on assert) and refill the manuscript.
# Usage: PY=/path/to/python bash build.sh   (needs numpy, scipy, sympy; see requirements.txt)
set -euo pipefail
cd -- "$(dirname -- "$0")"
PY="${PY:-python3}"
mkdir -p results
"$PY" -I -B -O code/algebra_checks.py results
"$PY" -I -B -O code/universality_checks.py results
"$PY" -I -B -O code/looks_and_mixture_controls.py results
"$PY" -I -B -O code/finite_horizon_mc.py results 1        # ~7 CPU-minutes
"$PY" -I -B -O code/repair_efficiency.py results 4000 100000   # ~4 CPU-minutes
"$PY" -I -B -O code/ou_spectral.py results 0.08          # ~25 CPU-minutes
"$PY" -I -B code/fill_tables.py results manuscript_template.md manuscript.md references.md
