"""Repaired split interval versus existing anytime-valid intervals (Gaussian mean, sigma = 1).

Methods (two groups of n observations each, 2n in total):
  naive      : split-naive (Pace et al. eq. 9), fixed threshold log(1/alpha). Not time-uniform.
  rep-mix    : split statistic with boundary b(n) = l1 * g(n)^2 / n, g = Robbins normal-mixture
               boundary for a standard Gaussian walk (Theorem 6 reduction).
  rep-stitch : same reduction with the stitched (Howard et al. 2021 style) boundary.
  run-mle    : running-MLE / split-exact sequence of Wasserman, Ramdas & Balakrishnan (2020),
               pooled stream, plug-in = mean of previous observations (Pace et al. eq. 10 form).
  mix-pool   : Robbins normal-mixture confidence sequence on the pooled 2n observations.
  stitch-pool: stitched LIL confidence sequence on the pooled 2n observations.
  checkpoint : model-free repair: split-naive set evaluated only at n_k = 2^k with threshold
               log(1/(alpha w_k)) (Markov + union bound; equivalently a checkpoint-weighted case of the
               time-mixture e-process of Ramdas and Wang 2025, Sec. 7.9), held constant between checkpoints.
Reports half-widths and empirical time-uniform non-coverage over 1 <= n <= NMAX.

Run: python -I -B repair_efficiency.py OUTDIR [reps] [nmax]
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import check, expect_failure, finish  # noqa: E402

out = sys.argv[1] if len(sys.argv) > 1 else "."
REPS = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
NMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 100000
L1 = (np.sqrt(5) - 1) / 4
RHO = 1.0  # prior precision in sum units: N(0,1) prior on theta as in Pace et al.


def g2_mixture(m, alpha, rho=RHO):
    """Squared Robbins normal-mixture boundary for |W_m| (two-sided, valid for all m>=1)."""
    return (m + rho) * (2 * np.log(1 / alpha) + np.log1p(m / rho))


def g2_stitched(m, alpha, eta=2.0, delta=0.5):
    """Squared stitched boundary for |W_m| (r = 1): 2 m (1+1/c_k) L_k on eta^k <= m < eta^{k+1}."""
    m = np.asarray(m, dtype=float)
    k = np.floor(np.log(m) / np.log(eta) + 1e-12)
    ck = 1 + np.log(k + 1)
    wk = np.log(k + np.e) ** (-delta) - np.log(k + 1 + np.e) ** (-delta)
    Lk = np.log(1 / (alpha * wk)) + 0.5 * np.log1p(eta * ck)
    return 2 * m * (1 + 1 / ck) * Lk


def b_split(n, alpha, kind):
    g2 = g2_mixture(n, alpha) if kind == "mix" else g2_stitched(n, alpha)
    return L1 * g2 / n


res = {"rho": RHO, "reps": REPS, "nmax": NMAX}
# ---- exact checks of the boundary algebra
w = np.log(np.arange(0, 10**6) + np.e) ** (-0.5) - np.log(np.arange(1, 10**6 + 1) + np.e) ** (-0.5)
check(np.all(w > 0) and abs(w.sum() - (1 - np.log(10**6 + np.e) ** -0.5)) < 1e-9, "stitching weights positive, telescoping")
# asymptotic form b(n) = 2 l1 loglog n + l1 (3 + 2 delta) logloglog n + O(1): bounded remainder
ns = np.array([1e4, 1e8, 1e16, 1e32, 1e64, 1e128])
rem = b_split(ns, 0.05, "stitch") - (2 * L1 * np.log(np.log(ns)) + L1 * (3 + 2 * 0.5) * np.log(np.log(np.log(ns))))
res["stitched_remainder"] = rem.tolist()
check(np.ptp(rem[2:]) < 2.0, "stitched boundary remainder bounded (O(1)) over n = 1e16..1e128", str(np.round(rem, 3)))

# ---- pathwise reduction Q_n <= l1 W_n^2 / n with W = projection on the top eigenvector
rng = np.random.default_rng(5)
A = np.array([[0, 0.5], [0.5, -0.5]])
wv, Vv = np.linalg.eigh(A); u = Vv[:, 1]
check(abs(wv[1] - L1) < 1e-15, "top eigenvalue of A is (sqrt5-1)/4")
Z = rng.normal(size=(20000, 2)).cumsum(axis=0)
nn = np.arange(1, 20001)
Q = (Z[:, 0] * Z[:, 1] - 0.5 * Z[:, 1] ** 2) / nn
W = Z @ u
check(np.all(Q <= L1 * W**2 / nn + 1e-9), "pathwise Q_n <= l1 W_n^2/n")
Wbad = Z @ np.array([1.0, 0.0])
expect_failure(np.all(Q <= L1 * Wbad**2 / nn + 1e-9), "wrong direction (first coordinate) does not dominate Q")

# ---- simulation
ALPHAS = [0.05, 0.20]
checkpoints = np.array([10, 100, 1000, 10000, 100000])
checkpoints = checkpoints[checkpoints <= NMAX]
methods = ["naive", "rep-mix", "rep-stitch", "run-mle", "mix-pool", "stitch-pool", "checkpoint"]
miss = {a: {mm: 0 for mm in methods} for a in ALPHAS}
widths = {a: {mm: {int(c): [] for c in checkpoints} for mm in methods} for a in ALPHAS}
n = np.arange(1, NMAX + 1, dtype=float)
m2 = 2 * n
# latest checkpoint 2^k <= n and its stitching weight (delta = 1/2)
kk = np.floor(np.log2(n) + 1e-12).astype(int)
nk = (2 ** kk).astype(int)
wk_cp = np.log(kk + np.e) ** (-0.5) - np.log(kk + 1 + np.e) ** (-0.5)
check(abs(sum(np.log(k + np.e) ** (-0.5) - np.log(k + 1 + np.e) ** (-0.5) for k in range(200000)) - (1 - np.log(200000 + np.e) ** (-0.5))) < 1e-9, "checkpoint weights telescope to <= 1")
pre = {a: {"mix": b_split(n, a, "mix"), "stitch": b_split(n, a, "stitch"),
           "gmix": np.sqrt(g2_mixture(m2, a)), "gst": np.sqrt(g2_stitched(m2, a))} for a in ALPHAS}
for rep in range(REPS):
    x = rng.standard_normal(NMAX); y = rng.standard_normal(NMAX)
    Sx, Ty = np.cumsum(x), np.cumsum(y)
    xb, yb = Sx / n, Ty / n
    # pooled stream interleaved x1, y1, x2, y2, ...
    z = np.empty(2 * NMAX); z[0::2] = x; z[1::2] = y
    cz = np.cumsum(z); cz2 = np.cumsum(z * z)
    t = np.arange(1, 2 * NMAX + 1, dtype=float)
    prev_mean = np.concatenate([[np.nan], cz[:-1] / t[:-1]])  # plug-in for i>=2
    dterm = np.concatenate([[0.0], np.cumsum((z[1:] - prev_mean[1:]) ** 2)])  # sum_{i=2}^t
    # run-mle at even totals t = 2n: log R(theta) = [ (t-1) theta^2 - 2 theta s + q ]/2 - D/2, over i = 2..t
    te = t[1::2]; s = cz[1::2] - z[0]; q = cz2[1::2] - z[0] ** 2; D = dterm[1::2]; tt = te - 1
    centre_rm = s / tt
    for a in ALPHAS:
        c = np.log(1 / a)
        half = {}
        diff2 = (xb - yb) ** 2
        half["naive"] = np.sqrt(diff2 + 2 * c / n)
        half["rep-mix"] = np.sqrt(diff2 + 2 * pre[a]["mix"] / n)
        half["rep-stitch"] = np.sqrt(diff2 + 2 * pre[a]["stitch"] / n)
        # run-mle: (tt theta^2 - 2 theta s + q - D)/2 < c  ->  theta in centre +- sqrt(...)
        rad = (s / tt) ** 2 - (q - D - 2 * c) / tt
        half["run-mle"] = np.sqrt(np.maximum(rad, 0.0))
        half["mix-pool"] = pre[a]["gmix"] / m2
        half["stitch-pool"] = pre[a]["gst"] / m2
        ck = np.log(1 / (a * wk_cp))
        half["checkpoint"] = np.sqrt(diff2[nk - 1] + 2 * ck / nk)
        centres = {"naive": xb, "rep-mix": xb, "rep-stitch": xb, "run-mle": centre_rm,
                   "mix-pool": cz[1::2] / m2, "stitch-pool": cz[1::2] / m2, "checkpoint": xb[nk - 1]}
        for mm in methods:
            cover = np.abs(centres[mm]) < half[mm]
            if mm == "run-mle":
                cover = cover & (rad > 0)
            if not cover.all():
                miss[a][mm] += 1
            for cpt in checkpoints:
                widths[a][mm][int(cpt)].append(2 * half[mm][cpt - 1])
out_tab = {}
for a in ALPHAS:
    out_tab[str(a)] = {}
    for mm in methods:
        out_tab[str(a)][mm] = {"noncoverage_1_to_nmax": miss[a][mm] / REPS,
                               "median_width": {str(k): float(np.median(v)) for k, v in widths[a][mm].items()}}
    print(a, {mm: out_tab[str(a)][mm]["noncoverage_1_to_nmax"] for mm in methods})
    for cpt in checkpoints:
        print("  n=", cpt, {mm: round(out_tab[str(a)][mm]["median_width"][str(int(cpt))], 5) for mm in methods})
res["table"] = out_tab


def ub(p, m):  # one-sided 99.9% upper Wilson-type allowance
    return p + 3.1 * np.sqrt(max(p * (1 - p), 1.0 / m) / m)


for a in ALPHAS:
    for mm in ["rep-mix", "rep-stitch", "run-mle", "mix-pool", "stitch-pool", "checkpoint"]:
        p = out_tab[str(a)][mm]["noncoverage_1_to_nmax"]
        check(p <= ub(a, REPS), f"alpha={a}: {mm} time-uniform non-coverage <= alpha (MC)", f"{p:.4f}")
    p = out_tab[str(a)]["naive"]["noncoverage_1_to_nmax"]
    res.setdefault("naive_noncoverage", {})[str(a)] = p
    if a == 0.20:  # at alpha=0.05 the finite-horizon excess is small (reported, not a control)
        check(p > a + 3.1 * np.sqrt(a * (1 - a) / REPS), f"negative control alpha={a}: split-naive exceeds alpha over 1..{NMAX}", f"{p:.4f}")
# asymptotic width ratio repaired-split / pooled-stitched -> sqrt(4 l1) * sqrt(2)/sqrt(2) ... computed exactly
# leading terms: rep: 2*b/n ~ 4 l1 loglog n / n ; pooled: g^2/(2n)^2 ~ 2(2n) loglog(2n)/(4n^2) = loglog/n
res["asymptotic_halfwidth_ratio_rep_over_pooled_LIL"] = float(np.sqrt(4 * L1))  # stitched repair vs pooled stitched, in probability
check(abs(np.sqrt(4 * L1) - np.sqrt(np.sqrt(5) - 1)) < 1e-15, "asymptotic ratio sqrt(sqrt5-1) ~ 1.112")
finish(res, os.path.join(out, "repair_efficiency.json"))
