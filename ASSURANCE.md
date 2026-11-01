# Assurance boundary

## What was done

1. **Written proofs** in the manuscript ([split_likelihood_persistence.pdf](split_likelihood_persistence.pdf); source [manuscript.md](manuscript.md), generated from [manuscript_template.md](manuscript_template.md) by [code/fill_tables.py](code/fill_tables.py) so that no tabulated number is typed by hand).
2. **Computational checks.** Six certificate scripts ([VERIFICATION.md](VERIFICATION.md)): exact symbolic identities and eigenvalues (Theorems 1, 3 and 4, the sandwich constant); the quadratic expansion of Theorem 2 on Poisson and Bernoulli streams up to n = 10^8; exact Monte Carlo of Pace et al.'s Gaussian and Poisson experiments; an Ornstein-Uhlenbeck spectral and survival-curve solver with Richardson extrapolation, finer-mesh, time-step and direct checks; exact simulation of finitely many looks and of the mixture-row controls; validity and widths of the repaired and comparison intervals. All Monte Carlo is seeded. Every check raises an error on failure (no `assert`), so it remains active under `python -O`.
3. **Negative controls** inside the scripts must be detected on every run: a wrong constant r/(1+r); an oracle statistic; a non-dominating direction; the split-naive excess at alpha = 0.20; a wrong Fisher information; the Markov bound for single looks.
4. **Internal AI-assisted review.** Five reviewer roles (editor, methodology, domain, perspective, devil's advocate) run as Claude Code subagents; an adversarial review by a separate model (OpenAI Codex, "Astra", through a local companion tool, read-only); a clean replay of the round-0 package by the lead (all outputs identical); one revision round; a re-review of the revision; a separate proof audit of the new and rewritten proofs; and a second Astra audit of those proofs. A rating-free summary is in [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md).
5. **Release replay.** Every script was rerun with `-I -B -O` from a fresh copy of the release package; every regenerated JSON output was compared with the shipped one (identical apart from the `_environment` block), and the manuscript was regenerated from the replayed results and compared byte for byte; see [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md). CI reruns the quick checks on Linux as a blocking job and the full suite as a non-blocking job.
6. **Citations.** Every reference with a DOI or arXiv identifier was copied from its Crossref or arXiv record ([crossref/](crossref/), [refs.bib](refs.bib)); DOIs were checked against the registries at release. Quotations and table values of Pace et al. were checked against the publisher's open-access full text (CC BY 4.0).

## What was not done

- No human specialist or journal peer review, no independent reproduction by an unaffiliated party, no proof-assistant formalisation.
- The theorems rest on the written proofs; scripts corroborate formulas, constants and implementations. Proposition 0, Lemma 1 and Theorem 6 are analytic and have no script.
- The continuum (OU) results of Section 7 are numerical with empirical, not certified, error estimates; the 99% crossing horizon is a power-law extrapolation.
- **Priority is bounded** by the documented searches in [novelty_report.md](novelty_report.md). The law of the iterated logarithm, the Hewitt-Savage zero-one law, Ville's inequality, Robbins's normal mixtures, stitching, the Kolmogorov-Erdős-Feller test, Vovk-Wang e-merging and the Ramdas-Wang time mixture are classical and credited where used. Ramdas and Wang (2025, Section 7.7) remarked without proof that split likelihood-ratio e-variables are not e-processes.
- The Poisson discrepancy with Pace et al.'s Table 5 is diagnosed (it also affects their non-split mixture row) but not resolved.
- Agreement between AI reviewers is not independent evidence.
