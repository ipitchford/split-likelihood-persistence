# Agent navigation

## Identity

- Release slug: `split-likelihood-persistence`. Version `0.1.0-candidate`, tag `v0.1.0-candidate`. DOI 10.5281/zenodo.23856728; concept DOI 10.5281/zenodo.23856727.
- Repository: https://github.com/ipitchford/split-likelihood-persistence. Release page: https://evidencepress.org/releases/split-likelihood-persistence/
- Scholarly creator: Anonymous. Status: unrefereed candidate ([STATUS.md](STATUS.md)).
- Manuscript: [PDF](split_likelihood_persistence.pdf) (15 pages), Markdown [manuscript.md](manuscript.md) generated from [manuscript_template.md](manuscript_template.md) and [references.md](references.md) by [code/fill_tables.py](code/fill_tables.py); generated LaTeX [split_likelihood_persistence.tex](split_likelihood_persistence.tex); builder [build_pdf.sh](build_pdf.sh) with [paper_header.tex](paper_header.tex); bibliography [refs.bib](refs.bib) with registry records in [crossref/](crossref/).
- Machine-readable claims: [CLAIMS.json](CLAIMS.json). Overview: [README.md](README.md).
- Assurance limits: [ASSURANCE.md](ASSURANCE.md). Review record (rating-free): [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). Provenance: [PROVENANCE.md](PROVENANCE.md).
- Prior art: [RESEARCH_GATES.json](RESEARCH_GATES.json), [novelty_report.md](novelty_report.md).
- Citation: [CITATION.cff](CITATION.cff), [.zenodo.json](.zenodo.json). Licences: [LICENSES.md](LICENSES.md), [LICENSE](LICENSE), [LICENSE-MIT.txt](LICENSE-MIT.txt), [PUBLIC_DOMAIN.md](PUBLIC_DOMAIN.md).

## Scientific map

Object: the universal-inference split likelihood ratio E_n(theta) (validation stream over a model fitted on a separate training stream), recomputed after each new pair of observations, with the fixed-threshold set C_n = {theta : log E_n(theta) < log(1/epsilon)} (the "split-naive" sequence of Pace, Salvan and Sartori, 2023, Section 3). Question: does it have asymptotic persistence at least 1 - epsilon, i.e. lim_n P(theta in C_m for all m >= n) >= 1 - epsilon? Answer: no, it is 0.

- At a glance (before Section 1): plain-English problem, answer, and a severity table by level (80%, 95%, 99%).
- Section 1: the conjecture verbatim, contributions, what is classical, prior-art search.
- Section 2: two-stream scheme; Proposition 0 (Hewitt-Savage refutation in any model); Lemma 1 (vector LIL for 2-homogeneous functionals).
- Section 3: Theorem 1 (Gaussian: Q_n = (S_n T_n - T_n^2/2)/n exactly; limsup/log log n = (sqrt5 - 1)/2).
- Section 4: Theorem 2 (regular models; sandwich constant; information identity R4).
- Section 5: Theorem 3 (split ratio, c(r)); Theorem 4 (fixed-K cross-fitting, aggregation classification).
- Section 6: Theorem 5 (threshold repair lambda_1 g(n)^2/n), Theorem 6 (sharpness, kappa > 3), Proposition 5A (model-free checkpoint repair).
- Section 7: Proposition 7 (exact OU representation); Table 1 (failure exponents, crossing horizons); Tables 2-3 (reproduction of Pace et al.; calibrated thresholds); Table 4 (finite looks); Section 7.5 (Poisson discrepancy diagnosis).
- Section 8: Table 5 (efficiency of repaired and existing anytime-valid intervals). Section 9: discussion and limitations. Section 10: scripts.

## Evidence and replay

- Claim-to-evidence map and negative controls: [VERIFICATION.md](VERIFICATION.md). Fast checks: [quick_checks.sh](quick_checks.sh) (under 2 minutes). Full replay: [build.sh](build.sh) (about 40 CPU-minutes). Exact comparer: [code/compare_results.py](code/compare_results.py); cross-platform comparer for CI: [tools/ci_compare.py](tools/ci_compare.py). Environment: [requirements.txt](requirements.txt), [ENVIRONMENT.txt](ENVIRONMENT.txt). CI: [.github/workflows/replay.yml](.github/workflows/replay.yml). Records: [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), [REPLAY_RECEIPT_round0.md](REPLAY_RECEIPT_round0.md). Integrity: [SHA256SUMS](SHA256SUMS).
- Scripts and outputs: [code/algebra_checks.py](code/algebra_checks.py) -> [results/algebra_checks.json](results/algebra_checks.json); [code/universality_checks.py](code/universality_checks.py) -> [results/universality_checks.json](results/universality_checks.json); [code/looks_and_mixture_controls.py](code/looks_and_mixture_controls.py) -> [results/looks_and_mixture_controls.json](results/looks_and_mixture_controls.json); [code/finite_horizon_mc.py](code/finite_horizon_mc.py) -> [results/finite_horizon_mc.json](results/finite_horizon_mc.json); [code/repair_efficiency.py](code/repair_efficiency.py) -> [results/repair_efficiency.json](results/repair_efficiency.json); [code/ou_spectral.py](code/ou_spectral.py) -> [results/ou_spectral.json](results/ou_spectral.json); shared helpers [code/common.py](code/common.py). Run logs: `results/*_run.log`.

## Trust boundary

All review was internal and AI-assisted (five Claude Code reviewer roles, a separate-model adversarial review and proof audit, a separate proof audit, a clean replay, a re-review). There was no peer review, no independent reproduction and no formal verification. Theorems rest on written proofs; Proposition 0, Lemma 1 and Theorem 6 have no script. The OU continuum numbers are numerical with empirical error estimates, and the 99% crossing horizon is extrapolated. The LIL, zero-one law, Ville's inequality, normal mixtures and stitching are classical. Priority holds only within the documented bounded searches.
