"""Numerical check of the quadratic expansion in Theorem 2 (Poisson and Bernoulli streams).

Exact split log e-value versus [U' I^{-1} V - V' I^{-1} V / 2]/n along long simulated paths;
the theorem says the difference is O((log log n)^{3/2}/sqrt n) almost surely.
Run: python -I -B universality_checks.py OUTDIR
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import check, expect_failure, finish  # noqa: E402

out = sys.argv[1] if len(sys.argv) > 1 else "."
rng = np.random.default_rng(20261101)
NS = np.array([10**3, 10**4, 10**5, 10**6, 10**7, 10**8])
PATHS = 400
res = {"n": NS.tolist(), "paths": PATHS}


def poisson(lam, wrong_info=False):
    errs = np.zeros((PATHS, len(NS)))
    for p in range(PATHS):
        inc = np.diff(np.concatenate([[0], NS]))
        S0 = np.cumsum(rng.poisson(lam * inc)).astype(float)
        T1 = np.cumsum(rng.poisson(lam * inc)).astype(float)
        lhat = T1 / NS
        exact = S0 * np.log(lhat / lam) - NS * (lhat - lam)
        U = (S0 - NS * lam) / lam; V = (T1 - NS * lam) / lam
        Iinv = 1.0 if wrong_info else lam
        quad = (U * Iinv * V - 0.5 * V * Iinv * V) / NS
        errs[p] = np.abs(exact - quad)
    return errs


def bernoulli(pr):
    errs = np.zeros((PATHS, len(NS)))
    for p in range(PATHS):
        inc = np.diff(np.concatenate([[0], NS]))
        S0 = np.cumsum(rng.binomial(inc, pr)).astype(float)
        T1 = np.cumsum(rng.binomial(inc, pr)).astype(float)
        ph = T1 / NS
        exact = S0 * np.log(ph / pr) + (NS - S0) * np.log((1 - ph) / (1 - pr))
        v = pr * (1 - pr)
        U = (S0 - NS * pr) / v; V = (T1 - NS * pr) / v
        quad = (U * v * V - 0.5 * V * v * V) / NS
        errs[p] = np.abs(exact - quad)
    return errs


for name, errs in [("poisson_lambda_1", poisson(1.0)), ("poisson_lambda_0.3", poisson(0.3)),
                   ("bernoulli_p_0.5", bernoulli(0.5)), ("bernoulli_p_0.1", bernoulli(0.1))]:
    q99 = np.quantile(errs, 0.99, axis=0)
    mx = errs.max(axis=0)
    scale = np.log(np.log(NS)) ** 1.5 / np.sqrt(NS)
    res[name] = {"q99_abs_error": q99.tolist(), "max_abs_error": mx.tolist(), "max_over_rate": (mx / scale).tolist()}
    print(name, np.round(mx, 6))
    check(mx[-1] < 1e-2, f"{name}: max |exact - quadratic| at n=1e8 below 1e-2", f"{mx[-1]:.2e}")
    check(np.all(np.diff(q99) < 0), f"{name}: 99% error quantile decreasing in n")
    check(np.max(mx / scale) < 50, f"{name}: error / ((loglog n)^1.5/sqrt n) bounded", f"{np.max(mx/scale):.2f}")
bad = poisson(3.0, wrong_info=True)
expect_failure(bad.max(axis=0)[-1] < 1e-2, "Poisson(3) with wrong Fisher information does not match the expansion",
               f"{bad.max(axis=0)[-1]:.3f}")
finish(res, os.path.join(out, "universality_checks.json"))
