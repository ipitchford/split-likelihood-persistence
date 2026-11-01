# Clean replay, round 0 (2026-11-01)

Fresh copy of bundle_r0 (read-only original) at replay/fresh_r0; every script run with
env/bin/python -I -B -O (Python 3.12, numpy 2.1.3, scipy 1.14.1, sympy 1.13.3) from replay/ as cwd,
outputs to replay/out. Results compared with bundle_r0/results by compare_results.py (JSON identical apart
from the _environment block):

- algebra_checks.json IDENTICAL, ALL CHECKS PASSED
- universality_checks.json IDENTICAL, ALL CHECKS PASSED
- finite_horizon_mc.json IDENTICAL, ALL CHECKS PASSED
- repair_efficiency.json IDENTICAL, ALL CHECKS PASSED
- ou_spectral.json IDENTICAL, ALL CHECKS PASSED

All Monte Carlo is seeded, so identity is expected; the replay confirms the outputs are produced by the shipped
code under optimised Python (checks do not rely on assert).
