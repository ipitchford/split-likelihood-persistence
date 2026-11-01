"""Continuum (Ornstein-Uhlenbeck) analysis of the split-naive sequence.

X(t) = e^{-t/2} (B1(e^t), B2(e^t)) is a stationary 2-d OU process with generator
L f = (1/2) Laplacian f - (1/2) x . grad f, and Q_n = X(log n)^T A X(log n).
Rotate to the eigenbasis of A: Q = l1 u^2 + l2 v^2 with l1 = (sqrt5-1)/4, l2 = -(sqrt5+1)/4.
Persistence over n0 <= n <= R n0 tends (n0 -> infinity) to
    P(R; c) = P_stationary( Q(t) < c for all 0 <= t <= log R ),  c = log(1/eps).
We compute (i) the principal Dirichlet eigenvalue mu(c) of -L on D_c = {Q < c}, so that
P(R; c) ~ C(c) R^{-mu(c)}; (ii) P(R; c) itself by exponential time stepping; (iii) the
horizon-calibrated threshold c*(eps, R) solving 1 - P(R; c*) = eps.
Ground-state transform: f = e^{|x|^2/4} w turns -L into H = -1/2 Lap + |x|^2/8 - 1/2.
Curved Dirichlet boundary handled with Shortley-Weller differences (second order).

Run: python -I -B ou_spectral.py OUTDIR [h]
"""
import os
import sys
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import check, finish  # noqa: E402

out = sys.argv[1] if len(sys.argv) > 1 else "."
H0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.08
L1 = (np.sqrt(5) - 1) / 4
L2 = -(np.sqrt(5) + 1) / 4
VMAX = 8.0


def build(c, h, l1=L1, l2=L2, vmax=VMAX, full_plane=False):
    """Return (sparse operator for H, node coordinates, cell weight h^2)."""
    umax = np.sqrt((c + abs(l2) * vmax**2) / l1) + 2 * h if not full_plane else vmax
    us = np.arange(-umax, umax + h / 2, h)
    vs = np.arange(-vmax, vmax + h / 2, h)
    U, V = np.meshgrid(us, vs, indexing="ij")

    def inside(u, v):
        if full_plane:
            return (np.abs(u) < vmax) & (np.abs(v) < vmax)
        return (l1 * u * u + l2 * v * v < c) & (np.abs(v) < vmax)

    mask = inside(U, V)
    idx = -np.ones(U.shape, dtype=np.int64)
    idx[mask] = np.arange(mask.sum())
    rows, cols, vals = [], [], []
    nu, nv = U.shape

    def dist_u(u, v, sgn):
        # distance from (u,v) to boundary along +-u
        if full_plane:
            return vmax - abs(u) if sgn * u > 0 else h
        target = np.sqrt((c - l2 * v * v) / l1)  # c - l2 v^2 > 0 always
        return target - sgn * u

    def dist_v(u, v, sgn):
        cands = [vmax - sgn * v]
        if not full_plane and l1 * u * u - c > 0:
            vb = np.sqrt((l1 * u * u - c) / (-l2))
            for vv in (vb, -vb):
                dd = sgn * (vv - v)
                if dd > 0:
                    cands.append(dd)
        return min(cands)

    I, J = np.nonzero(mask)
    for i, j in zip(I, J):
        p = idx[i, j]; u, v = us[i], vs[j]
        diag = (u * u + v * v) / 8 - 0.5
        for axis in (0, 1):
            hm = hp = h
            nb = {}
            for sgn in (-1, 1):
                ii, jj = (i + sgn, j) if axis == 0 else (i, j + sgn)
                if 0 <= ii < nu and 0 <= jj < nv and mask[ii, jj]:
                    nb[sgn] = idx[ii, jj]
                else:
                    dd = dist_u(u, v, sgn) if axis == 0 else dist_v(u, v, sgn)
                    dd = min(max(dd, 1e-3 * h), h)
                    if sgn < 0:
                        hm = dd
                    else:
                        hp = dd
            # -1/2 second derivative, Shortley-Weller weights
            cm = -1.0 / (hm * (hm + hp)); cp_ = -1.0 / (hp * (hm + hp))
            diag += 1.0 / (hm * hp)
            if -1 in nb:
                rows.append(p); cols.append(nb[-1]); vals.append(cm)
            if 1 in nb:
                rows.append(p); cols.append(nb[1]); vals.append(cp_)
        rows.append(p); cols.append(p); vals.append(diag)
    n = int(mask.sum())
    Hm = sps.csr_matrix((vals, (rows, cols)), shape=(n, n))
    return Hm, U[mask], V[mask], h * h


def principal(c, h):
    Hm, u, v, w = build(c, h)
    val, vec = spla.eigs(Hm, k=1, sigma=0.0, which="LM")
    return float(np.real(val[0]))


TMAX = np.log(1e6)
TGRID = np.linspace(0.0, TMAX, 701)


DT = 0.01


def survival_curve(c, h, dt=DT):
    """P(no crossing on [0,T]) on TGRID, stationary start.

    Crank-Nicolson for w' = -H w after four backward-Euler quarter steps (Rannacher start),
    sampled on TGRID by linear interpolation in time.
    """
    Hm, u, v, w = build(c, h)
    g = np.exp(-(u * u + v * v) / 4)
    wt = g * w / (2 * np.pi)
    Id = sps.identity(Hm.shape[0], format="csc")
    Hc = Hm.tocsc()
    be = spla.splu(Id + (dt / 4) * Hc)
    cn = spla.splu(Id + (dt / 2) * Hc)
    rhs = (Id - (dt / 2) * Hc).tocsr()
    times, vals = [0.0], [float(g @ wt)]
    x = g.copy()
    for _ in range(4):
        x = be.solve(x)
    t = dt
    times.append(t); vals.append(float(x @ wt))
    while t < TMAX - 1e-12:
        x = cn.solve(rhs @ x)
        t += dt
        times.append(t); vals.append(float(x @ wt))
    return np.interp(TGRID, np.array(times), np.array(vals))


_CACHE = {}


def survival_rich(c):
    """Richardson-extrapolated (h = 2*H0/2... i.e. H0 and H0/2) survival curve on TGRID."""
    key = round(float(c), 12)
    if key not in _CACHE:
        Pc = survival_curve(c, H0)
        Pf = survival_curve(c, H0 / 2)
        _CACHE[key] = (Pf + (Pf - Pc) / 3.0, Pc, Pf)
    return _CACHE[key][0]


def R_cross(P, mu, target):
    """log10 of the horizon ratio at which uncoverage 1-P reaches target (C R^-mu beyond the grid)."""
    U = 1 - P
    if U[-1] >= target:
        k = int(np.argmax(U >= target))
        T = np.interp(target, U[k - 1:k + 1], TGRID[k - 1:k + 1])
        return T / np.log(10)
    return (TGRID[-1] + (np.log(P[-1]) - np.log(1 - target)) / mu) / np.log(10)


res = {"lambda1": L1, "lambda2": L2, "h": H0}
# --- control 1: full plane (truncated box) has bottom eigenvalue ~ 0 (constants are invariant)
Hf, _, _, _ = build(0.0, 0.1, full_plane=True)
v0 = float(np.real(spla.eigs(Hf, k=1, sigma=-0.25, which="LM")[0][0]))
check(abs(v0) < 1e-3, "control: untruncated OU operator has ground eigenvalue 0", f"{v0:.2e}")
# --- control 2: P(T=0) equals stationary probability of D_c (Gaussian quadrature check)
from scipy import integrate
EPS = [0.20, 0.10, 0.05, 0.01]
table = {}
for eps in EPS:
    c = np.log(1 / eps)
    mus = [principal(c, hh) for hh in (2 * H0, H0, H0 / 2)]
    rich = mus[2] + (mus[2] - mus[1]) / 3.0  # second-order Richardson
    conv = abs(mus[2] - mus[1]) / abs(mus[1] - mus[0]) if mus[1] != mus[0] else 0
    check(conv < 0.5, f"eps={eps}: eigenvalue converges under refinement", f"mus={mus}, ratio={conv:.3f}")
    Rs = np.array([2, 5, 10, 20, 50, 100, 200, 400, 1e3, 1e4, 1e6])
    Pr = survival_rich(c)
    Ps = np.interp(np.log(Rs), TGRID, Pr)
    P0 = float(Pr[0])
    # exact stationary mass of D_c
    f = lambda vv: (2 * __import__("scipy.stats", fromlist=["norm"]).norm.cdf(np.sqrt((c - L2 * vv * vv) / L1)) - 1) * np.exp(-vv * vv / 2) / np.sqrt(2 * np.pi)
    mass = integrate.quad(f, -12, 12, epsabs=1e-13)[0]
    check(abs(P0 - mass) < 2e-3, f"eps={eps}: P(T=0) matches stationary mass of D_c", f"{P0:.5f} vs {mass:.5f}")
    check(mass >= 1 - eps, f"eps={eps}: stationary mass >= 1-eps (Markov)", f"{mass:.5f}")
    check(all(Ps[i] >= Ps[i + 1] - 1e-9 for i in range(len(Ps) - 1)), f"eps={eps}: survival nonincreasing")
    # tail fit C R^{-mu} from the two longest horizons (relative agreement, no absolute slack)
    slope = -(np.log(Ps[-1]) - np.log(Ps[-2])) / (np.log(Rs[-1]) - np.log(Rs[-2]))
    check(abs(slope - rich) < 0.03 * rich, f"eps={eps}: survival decay rate equals eigenvalue (3% rel)", f"{slope:.6f} vs {rich:.6f}")
    Pc, Pf = _CACHE[round(float(c), 12)][1:]
    rel = float(np.max(np.abs((1 - Pf) - (1 - Pr)) / (1 - Pr[1:]).max()))
    table[str(eps)] = {"c": c, "mu_h": mus, "mu": rich, "P0_stationary_mass": mass,
                       "uncoverage_by_R": {str(int(R)): float(1 - p) for R, p in zip(Rs, Ps)},
                       "richardson_correction_rel_max": rel,
                       "log10_R_for_uncoverage_eps": float(R_cross(Pr, rich, eps)),
                       "log10_R_for_uncoverage_2eps": float(R_cross(Pr, rich, 2 * eps)),
                       "log10_R_for_uncoverage_half": float(R_cross(Pr, rich, 0.5))}
    print(eps, table[str(eps)]["mu"], table[str(eps)]["uncoverage_by_R"])
res["exponents"] = table


# --- independent convergence check of pi itself at the hardest case (eps = 0.01): finer mesh h = H0/4
c01 = np.log(100.0)
P4 = survival_curve(c01, H0 / 4)
Pr01 = survival_rich(c01)
Pf01 = _CACHE[round(float(c01), 12)][2]
P4r = P4 + (P4 - Pf01) / 3.0
U_r, U_4r = 1 - Pr01[-1], 1 - P4r[-1]
check(abs(U_r - U_4r) < 0.01 * U_4r, "eps=0.01: Richardson pi(1e6) agrees with finer-mesh Richardson (1% rel)",
      f"{U_r:.6f} vs {U_4r:.6f}")
res["convergence_check_eps0.01_R1e6"] = {"rich_h08_h04": U_r, "rich_h04_h02": U_4r, "h02_raw": float(1 - P4[-1])}
# time-step check at the coarse mesh
Pdt = survival_curve(c01, H0, dt=DT / 2)
Pc01 = _CACHE[round(float(c01), 12)][1]
check(abs((1 - Pdt[-1]) - (1 - Pc01[-1])) < 0.005 * (1 - Pc01[-1]), "time step dt=0.01 converged (0.5% rel)",
      f"{1-Pc01[-1]:.6f} vs {1-Pdt[-1]:.6f}")

# --- horizon-calibrated thresholds c*(eps, R): root of the Richardson uncoverage, PCHIP in c on a local grid
from scipy.interpolate import PchipInterpolator
cal = {}
for eps in EPS:
    L = np.log(1 / eps)
    cg = np.linspace(0.65 * L, L + 0.8, 10)
    curves = [survival_rich(cc) for cc in cg]
    cal[str(eps)] = {}
    for R in (10, 100, 400, 1e4, 1e6):
        T = np.log(R)
        lu = np.array([np.log(1 - np.interp(T, TGRID, cv)) for cv in curves])  # decreasing in c
        check(np.all(np.diff(lu) < 0), f"eps={eps}, R={int(R)}: uncoverage decreasing in c")
        check(lu[0] > np.log(eps) > lu[-1], f"eps={eps}, R={int(R)}: threshold bracketed")
        f = PchipInterpolator(cg, lu - np.log(eps))
        cs = brentq(lambda x: float(f(x)), cg[0], cg[-1], xtol=1e-8)
        cal[str(eps)][str(int(R))] = cs
    print("calibrated", eps, cal[str(eps)], flush=True)
    vals = [cal[str(eps)][k] for k in ("10", "100", "400", "10000", "1000000")]
    check(all(vals[i] < vals[i + 1] for i in range(4)), f"eps={eps}: calibrated threshold grows with horizon")
# direct verification of two calibrated thresholds (fresh Richardson solves at c*)
for eps, R in ((0.01, 1e6), (0.2, 100)):
    cs = cal[str(eps)][str(int(R))]
    u = 1 - np.interp(np.log(R), TGRID, survival_rich(cs))
    check(abs(u - eps) < 0.005 * eps, f"direct check: uncoverage at c*({eps},{int(R)}) equals eps (0.5% rel)", f"{u:.6f}")
# consistency with Table 1: c*(eps,R) < log(1/eps) iff pi(R, log(1/eps)) < eps
for eps in EPS:
    for R in ("10", "100", "400", "10000", "1000000"):
        u = table[str(eps)]["uncoverage_by_R"][R]
        check((cal[str(eps)][R] < np.log(1 / eps)) == (u < eps), f"eps={eps}, R={R}: Table 1 and calibration agree on conservativeness")
res["calibrated_thresholds"] = cal
finish(res, os.path.join(out, "ou_spectral.json"))
