# The recomputed split likelihood-ratio interval is not a confidence sequence

Evidence Press research release, 1 November 2026. Version 0.1.0-candidate. DOI: [10.5281/zenodo.23856728](https://doi.org/10.5281/zenodo.23856728) (concept DOI, all versions: [10.5281/zenodo.23856727](https://doi.org/10.5281/zenodo.23856727)). Release page: <https://evidencepress.org/releases/split-likelihood-persistence/>. Author: Anonymous. Not externally peer reviewed. UK English. Text, data and diagrams CC0-1.0; code MIT.

**In one sentence:** a statistics method that is valid if you look once silently fails if you keep checking it as data arrive.

## In plain words

Universal inference (Wasserman, Ramdas and Balakrishnan, 2020) gives confidence intervals that are valid at one fixed sample size in almost any statistical model: split the sample in two, fit on one half, test candidate values on the other. Pace, Salvan and Sartori (Canadian Journal of Statistics 51, 2023, Section 3) conjectured that recomputing these intervals after every new pair of observations keeps a long-run guarantee ("asymptotic persistence at least 1 − ε"). This package shows that the long-run guarantee is zero, gives the exact constants, two repairs, and measures how bad the problem is in practice.

The damage is real but slow. The error first exceeds its nominal level at a last look about 100 times later than the first look for an 80% interval, about 10^4 times later at 95% and about 10^9 times later at 99%.

## Main results

- **Proposition 0 (refutation).** In any model, with any permutation-symmetric training estimator, the asymptotic persistence is 0 whenever the fixed-sample non-coverage does not vanish (Hewitt–Savage zero–one law).
- **Theorem 1 (exact Gaussian rate).** log E_n = (S_n T_n − T_n²/2)/n exactly for two independent standardised random walks, so limsup log E_n / log log n = (√5 − 1)/2 almost surely.
- **Theorem 2 (universality).** The same constant in every regular model satisfying the information identity; ((√5 − 1)/2)·μ_max(H⁻¹J) for a working likelihood.
- **Theorems 3–4 (splits and aggregation).** Split ratio r gives c(r) = (√(r² + 4r) − r)/2; fixed-K cross-fitting gives c(1/(K − 1)) for arithmetic averaging and 1/K for equal-weight geometric, harmonic or minimum aggregation.
- **Theorems 5–6, Proposition 5A (repairs).** The threshold λ₁g(n)²/n, with λ₁ = (√5 − 1)/4 and any valid scalar Gaussian boundary g, gives a genuine confidence sequence; within thresholds λ₁{2 log log n + κ log log log n + O(1)} validity holds iff κ > 3. A checkpoint repair works in every model.
- **Proposition 7 and Section 7 (severity).** An exact Ornstein–Uhlenbeck representation gives failure exponents 0.0421 (80%), 0.0054 (95%) and 0.0005 (99%), crossing horizons, calibrated thresholds and the effect of finitely many looks; Pace et al.'s Gaussian tables are reproduced.

**Classical versus new.** The law of the iterated logarithm, the zero–one law, Ville's inequality, normal mixtures and stitching are classical and credited. New, within a bounded search: the refutation of the conjecture, the constants, the aggregation classification, the transfer of the sharp threshold and the severity calculation.

**Not claimed.** Fixed-n universal inference, the running-MLE ("split-exact") construction and Pace et al.'s mixture intervals remain valid; their Section 6 conjecture is not addressed.

## Files

| File | Contents |
|---|---|
| `split_likelihood_persistence.pdf` | The paper (15 pages), complete proofs |
| `manuscript.md` | Source of the paper, generated from `manuscript_template.md` and `references.md` by `code/fill_tables.py`; `build_pdf.sh` and `paper_header.tex` rebuild the PDF; `split_likelihood_persistence.tex` is the generated LaTeX |
| `code/` | Certificate scripts (no `assert`; run under `python -O`; each has negative controls), `compare_results.py`, `fill_tables.py` |
| `results/` | JSON outputs and run logs |
| `quick_checks.sh`, `build.sh` | Fast blocking checks (under 2 minutes); full replay (about 40 CPU-minutes) |
| `tools/ci_compare.py` | Cross-platform comparison of regenerated results (used by CI) |
| `REVIEW_SUMMARY.md`, `VERIFICATION.md`, `novelty_report.md` | Rating-free review record, claim-to-evidence map, prior-art record |
| `STATUS.md`, `ASSURANCE.md`, `PROVENANCE.md`, `CLAIMS.json`, `RESEARCH_GATES.json`, `AI_INDEX.md` | Status, assurance boundary, provenance, machine-readable claims and research gates, agent index |
| `REPLAY_RECEIPT.md`, `REPLAY_RECEIPT_round0.md`, `ENVIRONMENT.txt` | Replay records and environment |
| `refs.bib`, `crossref/` | Bibliography of exactly the cited works, with registry records |

## Reproduce

```
python3 -m venv env && env/bin/pip install -r requirements.txt
PY=env/bin/python bash quick_checks.sh      # fast
PY=env/bin/python bash build.sh             # full (rewrites results/ and manuscript.md)
```

## Cite

Anonymous (2026). *The recomputed split likelihood-ratio interval is not a confidence sequence* (version 0.1.0-candidate). Evidence Press. https://doi.org/10.5281/zenodo.23856728
