# Review summary (public, AI-assisted internal review; not external peer review)

This candidate was checked on 1 November 2026 by an AI-assisted internal process. None of it is human peer review
or independent reproduction.

## Process
- Five reviewer roles (editor, methodology, domain, perspective, devil's advocate), each run as an independent
  Claude Code subagent on the round-0 manuscript, scripts and outputs.
- An adversarial mathematical review by a different model (Codex "Astra", read-only).
- A clean replay of every round-0 script from a fresh copy under optimised Python; all outputs were identical.
- One revision round, followed by a re-review of the changes, a separate proof audit of the new and rewritten proofs,
  and a second Astra audit of those proofs.

## Main objections and how they were resolved
1. *The universal constant needs the information identity* (Astra). (R1)-(R3) do not imply J = H; Astra gave an
   explicit mixture counterexample. Fix: an explicit assumption (R4); the sandwich constant is stated for the general case.
2. *The refutation does not need the LIL* (devil's advocate). Fix: added Proposition 0, a Hewitt-Savage zero-one-law
   refutation valid in any model with a symmetric training estimator. Theorems 1-4 are now presented as exact rate results.
3. *The headline overstated practical damage* (editor, perspective, devil's advocate). Fix: severity is now reported by level,
   with the horizon at which the error first exceeds nominal (about 100x at 80%, 10^3.9x at 95%, 10^8.8x at 99%),
   a finite-number-of-looks table, and an at-a-glance table.
4. *Calibrated thresholds at eps = 0.01 were not converged* (methodology, Astra). Fix: Richardson extrapolation of the
   survival curves, root-finding in c, finer-mesh and time-step checks, and two direct verifications. The corrected values
   agree with an independent solver used by the methodology reviewer.
5. *The repair was recommended beyond its proof* (Astra, devil's advocate). Fix: added Proposition 5A, a model-free checkpoint
   repair (dominated by the Ramdas-Wang time-mixture e-process). The Gaussian repair is now recommended only within its
   assumptions.
6. *The Poisson discrepancy was unexplained* (all reviewers). Fix: the same paper's Gaussian mixture row reproduces, but its
   Poisson mixture row does not. The discrepancy therefore also affects a non-split method; its source remains unresolved.
7. *Missing context* (domain, perspective). Fix: added group-sequential (Armitage, McPherson & Rowe; Pocock), boundary-crossing
   (Robbins & Siegmund; Lai) and fixed-n split-LRT literature, and corrected the Ramdas-Wang section references.
8. *Smaller items.* Theorem 3 algebra; aggregation scope (F = 1 is a symmetric e-merging function with persistence 1); the
   joint sub-Gaussian condition in Theorem 5; the mode of convergence of the asymptotic width ratio; several wording
   corrections.

## Verification outcome
- The re-review found all correctness items resolved and every number consistent with the result files. It flagged a few
  minor wording items, which were fixed.
- The proof audit found no open correctness issue.
- The second Astra audit found no CRITICAL or MAJOR defect, and its minor wording items were fixed.
- The refutation of the conjecture as literally stated (the source's equation (6), its definition of asymptotic persistence
  and its sampling scheme) was judged airtight by Astra and the devil's advocate.

## Limits
- The prior-art search was bounded.
- The numerical continuum results have empirical, not certified, error estimates.
- Growing-K or re-randomised cross-fitting and irregular-model rate constants are open.
