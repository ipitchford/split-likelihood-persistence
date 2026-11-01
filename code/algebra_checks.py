"""Exact algebra behind the LIL coefficients (Theorems 1, 3, 4 and Proposition 5).

Run: python -I -B algebra_checks.py OUTDIR
Everything symbolic is exact (sympy); numerical maximisations are cross-checks.
"""
import os
import sys
import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import check, expect_failure, finish  # noqa: E402

out = sys.argv[1] if len(sys.argv) > 1 else "."
res = {}

# ---------------------------------------------------------------- 1. Gaussian identity
n, a, b, th, sig = sp.symbols("n a b theta sigma", positive=True)
S, T = sp.symbols("S T", real=True)  # whitened-free sums: S = sum (X_i - theta)/sigma
# validation: n*a obs with sum of (X-theta)/sigma = S ; training mean deviation d = sigma*T/(n*b)
d = sig * T / (n * b)
# log-likelihood ratio of the validation sample at theta+d versus theta, sum over obs:
# sum[(x-theta) d - d^2/2]/sigma^2 = (sigma*S*d - n*a*d^2/2)/sigma^2
logE = sp.simplify((sig * S * d - n * a * d**2 / 2) / sig**2)
target_equal = (S * T - T**2 / 2) / n
check(sp.simplify(logE.subs({a: 1, b: 1}) - target_equal) == 0,
      "Gaussian split log e-value equals (S T - T^2/2)/n for equal groups")
# whitened unequal split: s = S/sqrt(n a), t = T/sqrt(n b)
s_, t_, r = sp.symbols("s t r", positive=True)
logE_w = sp.simplify(logE.subs({S: s_ * sp.sqrt(n * a), T: t_ * sp.sqrt(n * b)}))
Ar = sp.Matrix([[0, sp.sqrt(r) / 2], [sp.sqrt(r) / 2, -r / 2]])
quad = (sp.Matrix([[s_, t_]]) * Ar * sp.Matrix([s_, t_]))[0]
check(sp.simplify(logE_w - quad.subs(r, a / b)) == 0,
      "unequal split: whitened quadratic matrix is A_r with r=a/b")

# numerical check of the identity on random data
rng = np.random.default_rng(20261101)
for trial in range(5):
    na, nb = int(rng.integers(5, 60)), int(rng.integers(5, 60))
    theta, sigma = rng.normal(), float(rng.uniform(0.5, 2))
    x = rng.normal(theta, sigma, na)
    y = rng.normal(theta, sigma, nb)
    thhat = y.mean()
    direct = np.sum(-(x - thhat) ** 2 + (x - theta) ** 2) / (2 * sigma**2)
    Sx, Ty = np.sum(x - theta) / sigma, np.sum(y - theta) / sigma
    formula = Sx * Ty / nb - na * Ty**2 / (2 * nb**2)
    check(abs(direct - formula) < 1e-9 * max(1, abs(direct)), f"identity numeric trial {trial}",
          f"{direct:.12g} vs {formula:.12g}")

# ---------------------------------------------------------------- 2. eigenvalues and c(r)
lam = sp.symbols("lam")
cp = sp.factor((Ar - lam * sp.eye(2)).det())
lam_max = sp.simplify((-r + sp.sqrt(r**2 + 4 * r)) / 4)
check(sp.simplify(cp.subs(lam, lam_max)) == 0, "lambda_max(A_r) = (sqrt(r^2+4r)-r)/4")
c_r = sp.simplify(2 * lam_max)
check(sp.simplify(c_r - (sp.sqrt(r**2 + 4 * r) - r) / 2) == 0, "c(r) = (sqrt(r^2+4r)-r)/2")
golden = sp.nsimplify(c_r.subs(r, 1))
check(sp.simplify(golden - (sp.sqrt(5) - 1) / 2) == 0, "c(1) = (sqrt5-1)/2", str(golden))
res["c(1)"] = float(golden)
check(sp.limit(c_r, r, 0) == 0 and sp.limit(c_r, r, sp.oo) == 1, "c(r)->0 (r->0), ->1 (r->oo)")
dc = sp.simplify(sp.diff(c_r, r))
check(all(float(dc.subs(r, v)) > 0 for v in [1e-3, 0.1, 1, 10, 1e3]), "c(r) increasing (spot values)")
lam_min = sp.simplify((-r - sp.sqrt(r**2 + 4 * r)) / 4)
check(sp.simplify(cp.subs(lam, lam_min)) == 0, "lambda_min(A_r) = -(r+sqrt(r^2+4r))/4 < 0")
res["c_table"] = {str(v): float(c_r.subs(r, sp.Rational(v))) for v in ["1/9", "1/4", "1/3", "1/2", "1", "2", "3", "4", "9"]}
# wrong formula must be rejected
wrong = r / (1 + r)
expect_failure(sp.simplify(wrong.subs(r, 1) - golden) == 0, "wrong coefficient r/(1+r) differs from c(1)")

# ---------------------------------------------------------------- 3. two-fold cross-fitting
A = Ar.subs(r, 1)
Aswap = sp.Matrix([[-sp.Rational(1, 2), sp.Rational(1, 2)], [sp.Rational(1, 2), 0]])
G = (A + Aswap) / 2
ev = sorted(G.eigenvals().keys(), key=lambda z: float(z))
check(ev == [-sp.Rational(3, 4), sp.Rational(1, 4)], "geometric 2-fold eigenvalues -3/4, 1/4", str(ev))


# ---------------------------------------------------------------- 4. K-fold aggregation
def fold_matrices(K):
    """Whitened quadratic matrix of fold j: validate stream j against mean of the others."""
    mats = []
    for j in range(K):
        M = sp.zeros(K, K)
        others = [i for i in range(K) if i != j]
        for i in others:  # s_j * sum_others / (K-1), symmetrised
            M[j, i] += sp.Rational(1, 2 * (K - 1))
            M[i, j] += sp.Rational(1, 2 * (K - 1))
        for i in others:
            for k in others:
                M[i, k] -= sp.Rational(1, 2 * (K - 1) ** 2)
        mats.append(M)
    return mats


def max_min_on_sphere(mats, trials=4000, seed=1):
    """Numerically maximise min_j x^T M_j x over the unit sphere (cross-check only)."""
    from scipy.optimize import minimize
    Ms = [np.array(M, dtype=float) for M in mats]
    K = Ms[0].shape[0]
    g = np.random.default_rng(seed)
    best = -np.inf
    f = lambda x: -min(float(x @ M @ x) for M in Ms) / float(x @ x)
    for _ in range(trials // 100):
        x0 = g.normal(size=K)
        rr = minimize(f, x0, method="Nelder-Mead", options={"xatol": 1e-12, "fatol": 1e-14, "maxiter": 20000})
        best = max(best, -rr.fun)
    return best


kres = {}
for K in range(2, 8):
    mats = fold_matrices(K)
    rK = sp.Rational(1, K - 1)
    # each fold has whitened ratio r=1/(K-1): its positive eigenvalue
    lam_fold = max(float(v) for v in mats[0].eigenvals().keys())
    check(abs(2 * lam_fold - float(c_r.subs(r, rK))) < 1e-12, f"K={K}: fold coefficient = c(1/(K-1))",
          f"{2*lam_fold:.12f}")
    Gk = sum(mats, sp.zeros(K, K)) / K
    one = sp.ones(K, 1)
    AK = one * one.T / (2 * (K - 1) ** 2) - (2 * K - 1) * sp.eye(K) / (2 * K * (K - 1) ** 2)
    check(sp.simplify(Gk - AK) == sp.zeros(K, K), f"K={K}: geometric-mean matrix equals closed form A_K")
    evs = Gk.eigenvals()
    check(sp.Rational(1, 2 * K) in evs and -sp.Rational(2 * K - 1, 2 * K * (K - 1) ** 2) in evs,
          f"K={K}: A_K eigenvalues 1/(2K) and -(2K-1)/(2K(K-1)^2)", str(evs))
    # all-ones direction: every fold equals 1/(2K) on the unit vector
    u = one / sp.sqrt(K)
    vals = [sp.simplify((u.T * M * u)[0]) for M in mats]
    check(all(v == sp.Rational(1, 2 * K) for v in vals), f"K={K}: each fold = 1/(2K) on all-ones unit vector")
    mm = max_min_on_sphere(mats)
    # exact value 1/K: lower bound from the all-ones vector (above), upper bound because
    # min_j <= mean_j (geometric aggregator, eigenvalue 1/(2K)). Numerics must not exceed it.
    check(2 * mm <= 1.0 / K + 1e-9, f"K={K}: numeric max of min-aggregator never exceeds 1/K", f"{2*mm:.10f}")
    Ms = [np.array(M, dtype=float) for M in mats]
    xs = np.random.default_rng(K).normal(size=(200000, K)); xs /= np.linalg.norm(xs, axis=1, keepdims=True)
    mins = np.min(np.stack([np.einsum("ij,jk,ik->i", xs, M, xs) for M in Ms]), axis=0)
    check(2 * mins.max() <= 1.0 / K + 1e-12, f"K={K}: 2e5 random directions never exceed 1/K", f"{2*mins.max():.10f}")
    kres[K] = {"arithmetic_or_any_p_gt_0": float(c_r.subs(r, rK)), "geometric_or_p_le_0": 1.0 / K,
               "min_numeric": 2 * mm}
res["K_fold"] = kres

# ---------------------------------------------------------------- 5. Godambe / misspecified extension
# coefficient = c(1) * mu_max(H^{-1} J): check spectrum of A (x) J^{1/2} H^{-1} J^{1/2}
g = np.random.default_rng(7)
for trial in range(4):
    d_ = 3
    L1 = g.normal(size=(d_, d_)); J = L1 @ L1.T + d_ * np.eye(d_)
    L2 = g.normal(size=(d_, d_)); H = L2 @ L2.T + d_ * np.eye(d_)
    Hi = np.linalg.inv(H)
    Aflt = np.array(A, dtype=float)
    Mq = np.kron(Aflt, Hi)  # quadratic form in (U,V), U,V each Cov J per step
    Sig = np.kron(np.eye(2), J)
    w, Vv = np.linalg.eigh(Sig); Sh = Vv @ np.diag(np.sqrt(w)) @ Vv.T
    coef = 2 * max(0, np.linalg.eigvalsh(Sh @ Mq @ Sh).max())
    mu = np.linalg.eigvals(Hi @ J).real.max()
    check(abs(coef - float(golden) * mu) < 1e-10, f"Godambe extension trial {trial}", f"{coef:.12f} vs {float(golden)*mu:.12f}")
# correctly specified (J=H) gives c(1) regardless of dimension
check(abs(2 * np.linalg.eigvalsh(np.kron(np.array(A, dtype=float), np.eye(4))).max() - float(golden)) < 1e-14,
      "A (x) I_d has the same coefficient for any d")

# ---------------------------------------------------------------- 6. negative control: oracle statistic
# replacing the training estimate by the true theta gives log E = 0 exactly -> coefficient 0
check(sp.simplify(logE.subs(T, 0)) == 0, "oracle (T=0) statistic is identically zero")
expect_failure(2 * max(float(v) for v in A.eigenvals().keys()) <= 0, "split statistic would be harmless (it is not)")

finish(res, os.path.join(out, "algebra_checks.json"))
