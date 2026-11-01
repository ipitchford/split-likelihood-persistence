# Verification record

| Claim | Evidence |
|---|---|
| Thm 1(i) identity, eigenvalues of A | `algebra_checks.py` (exact symbolic, plus 5 numeric data sets) |
| Thm 3 c(r), its limits and monotonicity | `algebra_checks.py` (symbolic) |
| Thm 4 fold matrices, A_K spectrum, min-aggregator value 1/K | `algebra_checks.py` (symbolic for K = 2..7, plus numeric upper-bound checks) |
| Thm 2 sandwich constant | `algebra_checks.py` (random J, H) |
| Thm 2 quadratic expansion | `universality_checks.py` (Poisson, Bernoulli, n up to 1e8; error / rate bounded) |
| Thm 5 pathwise domination; stitched remainder O(1) | `repair_efficiency.py` |
| Thm 5 / Prop 5A validity, widths (Table 5) | `repair_efficiency.py` (4,000 replications, n up to 1e5) |
| Reproduction of Pace et al. Tables 2-3; Poisson Table 5 | `finite_horizon_mc.py` (200,000 / 100,000 sequences) |
| Prop 7 continuum: exponents, survival, calibrated c* (Tables 1, 3) | `ou_spectral.py` (Richardson in h, convergence and direct checks) |
| Looks table, mixture-row diagnostics (Table 4, Sec. 7.5) | `looks_and_mixture_controls.py` |
| Lemma 1, Prop 0, Thms 1-6, Prop 5A, Prop 7 proofs | written proofs; AI-assisted proof audit and Astra audit (see REVIEW_SUMMARY.md) |

Negative controls (must be detected): a wrong constant r/(1+r); an oracle statistic; a non-dominating direction; the
split-naive excess at alpha = 0.20; a wrong Fisher information; the Markov bound for single looks.
Proposition 0, Lemma 1 and Theorem 6 are analytic and have no script.
