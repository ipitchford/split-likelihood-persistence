"""(a) Split-naive uncoverage with a finite number of looks (group-sequential view), exact Gaussian.
(b) Diagnostic controls for the Poisson discrepancy: reproduce Pace et al.'s *mixture* rows
    (Gaussian Table 2 normal mixture; Poisson Table 5 Gamma(1,1) mixture) by exact simulation.

Run: python -I -B looks_and_mixture_controls.py OUTDIR
"""
import os
import sys
import numpy as np
from scipy.special import gammaln

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import check, finish  # noqa: E402

out = sys.argv[1] if len(sys.argv) > 1 else "."
EPS = np.array([0.20, 0.10, 0.05, 0.01]); CUT = np.log(1 / EPS)
res = {"eps": EPS.tolist()}
rng = np.random.default_rng(777)


def looks_uncoverage(n_first, R, L, reps=1_000_000):
    looks = np.unique(np.round(n_first * R ** (np.arange(L) / (L - 1))).astype(int))
    S = np.zeros(reps); T = np.zeros(reps); prev = 0
    mq = np.full(reps, -np.inf)
    for n in looks:
        S += rng.normal(0, np.sqrt(n - prev), reps); T += rng.normal(0, np.sqrt(n - prev), reps); prev = n
        np.maximum(mq, (S * T - 0.5 * T * T) / n, out=mq)
    return [(mq >= c).mean() for c in CUT], looks.tolist()


tab = {}
for (nf, R) in [(100, 100), (10, 1000)]:
    for L in (1, 2, 3, 5, 10, 20, 50):
        if L == 1:
            S = rng.normal(0, np.sqrt(nf), 1_000_000); T = rng.normal(0, np.sqrt(nf), 1_000_000)
            u = [((S * T - 0.5 * T * T) / nf >= c).mean() for c in CUT]; lk = [nf]
        else:
            u, lk = looks_uncoverage(nf, R, L)
        tab[f"first{nf}_R{R}_L{L}"] = {"uncoverage": u, "looks": lk}
        print(nf, R, L, np.round(np.array(u) * 100, 2), flush=True)
res["looks"] = tab
for key in ("first100_R100", "first10_R1000"):
    seq = [tab[f"{key}_L{L}"]["uncoverage"][0] for L in (1, 2, 3, 5, 10, 20, 50)]
    check(all(seq[i] <= seq[i + 1] + 0.002 for i in range(len(seq) - 1)), f"{key}: uncoverage increases with number of looks (MC tol)")
    check(seq[0] <= 0.20, f"{key}: single look respects eps=0.20 (Markov)", f"{seq[0]:.4f}")


# (b1) Gaussian normal-mixture CS (Pace eq. 8 with mu0=0, tau0^2=1), total sizes 200..80000 (Table 2)
def gauss_mixture(n0, N, reps):
    S = rng.normal(0, np.sqrt(2 * n0), reps); m = 2 * n0; mq = np.full(reps, -np.inf)
    while True:
        lr = -0.5 * np.log1p(m) + S * S / (2 * (1 + m))  # log q/p at theta=0, tau0^2=1
        np.maximum(mq, lr, out=mq)
        if m == 2 * N:
            break
        S += rng.normal(0, np.sqrt(2.0), reps); m += 2
    return [(mq >= c).mean() for c in CUT]


g = gauss_mixture(100, 40000, 40000)
res["gauss_mixture_table2"] = {"ours": g, "published_percent": [3.22, 1.79, 0.94, 0.21]}
print("gauss mixture", np.round(np.array(g) * 100, 2))
se = lambda p, m: np.sqrt(max(p * (1 - p), 1e-9) / m)
for e, ours, pub in zip(EPS, g, [3.22, 1.79, 0.94, 0.21]):
    z = abs(ours - pub / 100) / np.sqrt(se(ours, 40000) ** 2 + se(pub / 100, 10000) ** 2)
    check(z < 4, f"control: Gaussian mixture row of Table 2 reproduced, eps={e}", f"{100*ours:.2f} vs {pub} (z={z:.2f})")


# (b2) Poisson Gamma(1,1) mixture CS (Pace eq. 12), lambda = 1, both horizon readings of Table 5
def poisson_mixture(n0, N, reps):
    Sm = rng.poisson(2 * n0, reps).astype(float); m = 2 * n0; mq = np.full(reps, -np.inf)
    while True:
        lr = gammaln(1 + Sm) - (1 + Sm) * np.log(1 + m) + m  # log q/p at lambda=1, a=b=1
        np.maximum(mq, lr, out=mq)
        if m == 2 * N:
            break
        Sm += rng.poisson(2.0, reps); m += 2
    return [(mq >= c).mean() for c in CUT]


p1 = poisson_mixture(100, 10000, 20000); p2 = poisson_mixture(50, 5000, 20000)
res["poisson_mixture"] = {"n100_10000": p1, "n50_5000": p2, "published_percent": [9.20, 5.44, 3.34, 1.00]}
print("poisson mixture", np.round(np.array(p1) * 100, 2), np.round(np.array(p2) * 100, 2))
# Diagnostic (not a pass condition): does either reading reproduce the published Poisson mixture row?
res["poisson_mixture_reproduces_published"] = bool(min(abs(p1[0] - 0.092), abs(p2[0] - 0.092)) < 0.01)
check(max(p1[0], p2[0]) < 0.092 - 0.01, "diagnostic: Poisson mixture row also differs from published Table 5 (ours lower)",
      f"{100*p1[0]:.2f}, {100*p2[0]:.2f} vs 9.20")
finish(res, os.path.join(out, "looks_and_mixture_controls.json"))
