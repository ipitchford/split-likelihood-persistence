"""Exact Monte Carlo of the split-naive sequence at Pace, Salvan & Sartori's horizons.

Gaussian (known variance 1, theta = 0) and Poisson (lambda = 1) two-stream schemes:
each group grows by one fresh observation per step and groups never mix.
Records uncoverage (theta outside C_n for some n0<=n<=N) and incompatibility
(empty intersection of C_n over n0<=n<=N), plus uncoverage at intermediate horizons.

Run: python -I -B finite_horizon_mc.py OUTDIR [scale]
scale multiplies the number of replications (default 1).
"""
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import check, finish  # noqa: E402

out = sys.argv[1] if len(sys.argv) > 1 else "."
scale = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
EPS = np.array([0.20, 0.10, 0.05, 0.01])
CUT = np.log(1.0 / EPS)


def gauss_run(n0, N, reps, seed, batch=50000, ratios=(2, 5, 10, 20, 50, 100, 200, 400)):
    rng = np.random.default_rng(seed)
    unc = np.zeros(len(EPS)); inc = np.zeros(len(EPS))
    unc_at = {R: np.zeros(len(EPS)) for R in ratios if n0 * R <= N}
    done = 0
    while done < reps:
        B = min(batch, reps - done)
        # sums at n0 drawn exactly, then one step at a time
        Sx = rng.normal(0.0, np.sqrt(n0), B)  # validation group sum (sigma=1, theta=0)
        Ty = rng.normal(0.0, np.sqrt(n0), B)  # training group sum
        maxQ = np.full(B, -np.inf)
        lo = np.full((len(EPS), B), -np.inf); hi = np.full((len(EPS), B), np.inf)
        n = n0
        while True:
            xb, yb = Sx / n, Ty / n
            Q = (Sx * Ty - 0.5 * Ty * Ty) / n
            np.maximum(maxQ, Q, out=maxQ)
            half = np.sqrt((xb - yb) ** 2 + 2.0 * CUT[:, None] / n)
            np.maximum(lo, xb - half, out=lo); np.minimum(hi, xb + half, out=hi)
            if (n // n0) * n0 == n and (n // n0) in unc_at:
                unc_at[n // n0] += (maxQ[None, :] >= CUT[:, None]).sum(axis=1)
            if n == N:
                break
            z = rng.standard_normal((2, B))
            Sx += z[0]; Ty += z[1]; n += 1
        unc += (maxQ[None, :] >= CUT[:, None]).sum(axis=1)
        inc += (lo > hi).sum(axis=1)
        done += B
    return {"n0": n0, "N": N, "reps": reps,
            "uncoverage": (unc / reps).tolist(), "incompatibility": (inc / reps).tolist(),
            "uncoverage_by_ratio": {str(R): (v / reps).tolist() for R, v in unc_at.items()}}


def poisson_run(n0, N, reps, seed, lam=1.0, batch=50000):
    rng = np.random.default_rng(seed)
    unc = np.zeros(len(EPS)); done = 0
    while done < reps:
        B = min(batch, reps - done)
        S0 = rng.poisson(lam * n0, B).astype(float)  # validation sum
        T1 = rng.poisson(lam * n0, B).astype(float)  # training sum
        maxQ = np.full(B, -np.inf); n = n0
        while True:
            lhat = np.where(T1 > 0, T1 / n, 0.5 / n)
            Q = S0 * np.log(lhat / lam) - n * (lhat - lam)
            np.maximum(maxQ, Q, out=maxQ)
            if n == N:
                break
            S0 += rng.poisson(lam, B); T1 += rng.poisson(lam, B); n += 1
        unc += (maxQ[None, :] >= CUT[:, None]).sum(axis=1)
        done += B
    return {"n0": n0, "N": N, "reps": reps, "lambda": lam, "uncoverage": (unc / reps).tolist()}


def se(p, m):
    return np.sqrt(np.maximum(p * (1 - p), 1e-12) / m)


res = {"eps": EPS.tolist()}
t0 = time.time()
R1 = int(200000 * scale)
g1 = gauss_run(100, 40000, R1, 11)
res["gauss_table2_n100_40000"] = g1
print("gauss table2", g1["uncoverage"], g1["incompatibility"], f"{time.time()-t0:.0f}s")
g2 = gauss_run(400, 80000, int(100000 * scale), 12)
res["gauss_table3_n400_80000"] = g2
print("gauss table3", g2["uncoverage"], g2["incompatibility"], f"{time.time()-t0:.0f}s")
p1 = poisson_run(100, 10000, int(100000 * scale), 13)
p2 = poisson_run(50, 5000, int(100000 * scale), 14)
res["poisson_n100_10000"] = p1; res["poisson_n50_5000"] = p2
print("poisson", p1["uncoverage"], p2["uncoverage"], f"{time.time()-t0:.0f}s")

# Pace et al. (2023) published Monte Carlo percentages (Tables 2, 3, 5), 10,000 / 5,000 sequences
pace = {"t2_unc": [23.66, 9.15, 3.35, 0.37], "t2_inc": [7.82, 2.09, 0.54, 0.03],
        "t3_unc": [21.67, 8.30, 3.36, 0.28], "t3_inc": [5.77, 1.50, 0.41, 0.02],
        "t5_unc": [43.00, 20.34, 8.76, 1.34]}
res["pace_published_percent"] = pace


def agree(ours, m_ours, theirs_pct, m_theirs, label):
    for e, po, pt in zip(EPS, ours, theirs_pct):
        pt = pt / 100.0
        z = abs(po - pt) / np.sqrt(se(po, m_ours) ** 2 + se(pt, m_theirs) ** 2)
        check(z < 4.0, f"{label} eps={e}: ours {100*po:.2f}% vs published {100*pt:.2f}%", f"z={z:.2f}")


agree(g1["uncoverage"], R1, pace["t2_unc"], 10000, "Table 2 uncoverage")
agree(g1["incompatibility"], R1, pace["t2_inc"], 10000, "Table 2 incompatibility")
agree(g2["uncoverage"], g2["reps"], pace["t3_unc"], 10000, "Table 3 uncoverage")
agree(g2["incompatibility"], g2["reps"], pace["t3_inc"], 10000, "Table 3 incompatibility")
# Poisson: equation (13) as printed does NOT reproduce Table 5 under either horizon reading.
# Record the z-scores (a documented discrepancy, not a pass condition), and check instead the
# universality prediction: Poisson uncoverage is close to the Gaussian value at the same
# horizon ratio (Gaussian n0=100, ratio 100).
zs = []
for p in (p1, p2):
    zs.append(max(abs(po - pt / 100) / np.sqrt(se(po, p["reps"]) ** 2 + se(pt / 100, 5000) ** 2)
                  for po, pt in zip(p["uncoverage"], pace["t5_unc"])))
res["poisson_max_z_vs_table5"] = {"n100_10000": zs[0], "n50_5000": zs[1]}
print("Poisson vs published Table 5, max z:", zs)
g100 = g1["uncoverage_by_ratio"]["100"]
res["gauss_ratio100_from_n100"] = g100
for e, pp, gg in zip(EPS, p1["uncoverage"], g100):
    check(abs(pp - gg) < 0.02 + 4 * se(gg, R1), f"Poisson (n 100..10000) close to Gaussian ratio-100 value, eps={e}",
          f"{100*pp:.2f}% vs {100*gg:.2f}%")
# monotone in horizon (sanity) and above eps at 20% for the shorter burn-in (Pace's anti-conservativeness)
byR = g1["uncoverage_by_ratio"]
ks = sorted(byR, key=int)
check(all(byR[ks[i]][0] <= byR[ks[i + 1]][0] for i in range(len(ks) - 1)), "uncoverage nondecreasing in horizon")
check(g1["uncoverage"][0] > 0.20, "eps=0.20 split-naive is anti-conservative at Pace's horizon")
# Negative control: a valid fixed-n statement. At a single n the uncoverage must be <= eps (Markov).
rng = np.random.default_rng(99)
m = 400000; n = 1000
Sx = rng.normal(0, np.sqrt(n), m); Ty = rng.normal(0, np.sqrt(n), m)
Q = (Sx * Ty - 0.5 * Ty ** 2) / n
fixed = [(Q >= c).mean() for c in CUT]
res["fixed_n_noncoverage"] = fixed
check(all(f <= e for f, e in zip(fixed, EPS)), "negative control: fixed-n non-coverage <= eps", str(fixed))
finish(res, os.path.join(out, "finite_horizon_mc.json"))
