#!/usr/bin/env bash
# Fast blocking checks (< 2 minutes): exact algebra, universality expansion, looks/mixture controls.
set -euo pipefail
cd -- "$(dirname -- "$0")"
PY="${PY:-python3}"
mkdir -p quick_out
"$PY" -I -B -O code/algebra_checks.py quick_out
"$PY" -I -B -O code/universality_checks.py quick_out
"$PY" -I -B -O code/looks_and_mixture_controls.py quick_out
"$PY" -I -B code/compare_results.py quick_out results algebra_checks.json universality_checks.json looks_and_mixture_controls.json
