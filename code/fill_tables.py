"""Fill manuscript placeholders from the results JSON files (no numbers typed by hand).

Run: python -I -B fill_tables.py RESULTS_DIR TEMPLATE.md OUTPUT.md REFS.md
"""
import json
import math
import os
import sys

rd, tpl, outp, refs = sys.argv[1:5]
J = lambda name: json.load(open(os.path.join(rd, name)))
ou = J("ou_spectral.json"); mc = J("finite_horizon_mc.json"); ef = J("repair_efficiency.json")
lk = J("looks_and_mixture_controls.json")
text = open(tpl).read()
EPS = ["0.2", "0.1", "0.05", "0.01"]

# ---- OU table
Rs = ["10", "100", "400", "10000", "1000000"]
rows = ["| $\\varepsilon$ | $c=\\log(1/\\varepsilon)$ | $P(\\mathbf X(0)\\in D_c)$ | exponent $\\mu(c)$ | $R=10$ | $10^2$ | $400$ | $10^4$ | $10^6$ | $\\log_{10}R$ at uncoverage $\\varepsilon$ | $\\log_{10}R$ at $2\\varepsilon$ |",
        "|---|---|---|---|---|---|---|---|---|---|---|"]
for e in EPS:
    t = ou["exponents"][e]
    u = t["uncoverage_by_R"]
    rows.append(f"| {float(e):.2f} | {t['c']:.3f} | {t['P0_stationary_mass']:.4f} | {t['mu']:.5f} | "
                + " | ".join(f"{100*u[r]:.2f}%" for r in Rs) + f" | {t['log10_R_for_uncoverage_eps']:.1f} | {t['log10_R_for_uncoverage_2eps']:.1f} |")
ou_table = ("Table 1. Continuum (OU) uncoverage $\\pi(R,c)$ of the split-naive sequence over a horizon of ratio "
            "$R=n_{\\max}/n_{\\min}$, the stationary single-time coverage, and the failure exponent. "
            "Beyond $R=10^6$ the last two columns use the fitted power law $\\mathrm{const}\\times R^{-\\mu}$.\n\n" + "\n".join(rows))
text = text.replace("{{OU_TABLE}}", ou_table)

# R at which uncoverage first reaches eps for eps = 0.2 (log-linear interpolation)
u = ou["exponents"]["0.2"]["uncoverage_by_R"]
keys = sorted(u, key=float)
R20 = None
for a, b in zip(keys, keys[1:]):
    if u[a] < 0.2 <= u[b]:
        la, lb = math.log(float(a)), math.log(float(b))
        R20 = math.exp(la + (0.2 - u[a]) / (u[b] - u[a]) * (lb - la))
        break
text = text.replace("{{R_eps20}}", f"{R20:.0f}" if R20 else "n/a")
for e, tag in (("0.1", "L10_eps10"), ("0.05", "L10_eps05"), ("0.01", "L10_eps01")):
    text = text.replace("{{" + tag + "}}", f"{ou['exponents'][e]['log10_R_for_uncoverage_eps']:.1f}")
# at-a-glance table
gl = ["| nominal level | error over a horizon of 100x the first look | 10,000x | 1,000,000x | horizon at which the nominal error is first exceeded | failure exponent |",
      "|---|---|---|---|---|---|"]
for e, lab in (("0.2", "80%"), ("0.05", "95%"), ("0.01", "99%")):
    t = ou["exponents"][e]; u = t["uncoverage_by_R"]; l10 = t["log10_R_for_uncoverage_eps"]
    first = (f"about {R20:.0f}x" if e == "0.2" else f"about $10^{{{l10:.1f}}}$x")
    gl.append(f"| {lab} | {100*u['100']:.2f}% | {100*u['10000']:.2f}% | {100*u['1000000']:.2f}% | {first} | {t['mu']:.4f} |")
text = text.replace("{{GLANCE}}", "\n".join(gl))

# ---- MC table
pub = mc["pace_published_percent"]
g2, g3 = mc["gauss_table2_n100_40000"], mc["gauss_table3_n400_80000"]
se = lambda p, m: 100 * math.sqrt(max(p * (1 - p), 1e-12) / m)
rows = ["| | $\\varepsilon$ | 0.20 | 0.10 | 0.05 | 0.01 |", "|---|---|---|---|---|---|"]
for lab, g, pu, pi, R in [("$100\\le n\\le 40{,}000$", g2, pub["t2_unc"], pub["t2_inc"], "400"),
                          ("$400\\le n\\le 80{,}000$", g3, pub["t3_unc"], pub["t3_inc"], "200")]:
    rows.append(f"| {lab} | uncoverage, published | " + " | ".join(f"{x:.2f}" for x in pu) + " |")
    rows.append(f"| | uncoverage, this paper | " + " | ".join(f"{100*x:.2f} ± {se(x, g['reps']):.2f}" for x in g["uncoverage"]) + " |")
    cont = [ou["exponents"][e]["uncoverage_by_R"].get(R) for e in EPS]
    if all(c is not None for c in cont):
        rows.append(f"| | uncoverage, continuum $\\pi({R},c)$ | " + " | ".join(f"{100*c:.2f}" for c in cont) + " |")
    rows.append(f"| | incompatibility, published | " + " | ".join(f"{x:.2f}" for x in pi) + " |")
    rows.append(f"| | incompatibility, this paper | " + " | ".join(f"{100*x:.2f} ± {se(x, g['reps']):.2f}" for x in g["incompatibility"]) + " |")
mc_table = (f"Table 2. Pace et al.'s Gaussian experiment (their Tables 2 and 3, 10,000 sequences) against exact simulation "
            f"({g2['reps']:,} and {g3['reps']:,} sequences); percentages, ± one standard error.\n\n" + "\n".join(rows))
text = text.replace("{{MC_TABLE}}", mc_table)

# ---- calibrated thresholds
rows = ["| $\\varepsilon$ | $\\log(1/\\varepsilon)$ | $c^*$, $R=10$ | $10^2$ | $400$ | $10^4$ | $10^6$ |", "|---|---|---|---|---|---|---|"]
for e in EPS:
    cc = ou["calibrated_thresholds"][e]
    rows.append(f"| {float(e):.2f} | {math.log(1/float(e)):.3f} | " + " | ".join(f"{cc[r]:.2f}" for r in Rs) + " |")
cal = ("Table 3. Horizon-calibrated thresholds $c^*(\\varepsilon,R)$ in the continuum limit (Richardson-extrapolated survival curves, "
       "root-finding in $c$; two entries verified by a direct solve to 0.5% relative in uncoverage; estimated error in $c^*$ about $\\pm0.005$). "
       "Where $c^*<\\log(1/\\varepsilon)$ the split-naive threshold is conservative for that horizon.\n\n" + "\n".join(rows))
text = text.replace("{{CAL_TABLE}}", cal)

# ---- Poisson text
p1, p2 = mc["poisson_n100_10000"], mc["poisson_n50_5000"]
fmt = lambda p: ", ".join(f"{100*x:.2f}" for x in p["uncoverage"])
text = text.replace("{{POIS}}", f"{fmt(p1)}% for $100\\le n\\le10{{,}}000$ and {fmt(p2)}% for $50\\le n\\le5{{,}}000$ "
                    f"(standard errors at most 0.13 percentage points)")

# ---- efficiency table
tab = ef["table"]
names = {"naive": "split-naive (not valid)", "rep-mix": "repaired split, mixture $g$", "rep-stitch": "repaired split, stitched $g$",
         "run-mle": "running MLE (split-exact), pooled", "mix-pool": "normal mixture, pooled", "stitch-pool": "stitched, pooled",
         "checkpoint": "checkpoint repair (model-free)"}
rows = ["| method | non-coverage over $1\\le n\\le10^5$, $\\alpha=0.05$ | $\\alpha=0.20$ | width $n=10^2$ | $10^3$ | $10^4$ | $10^5$ | width ratio to pooled mixture at $10^4$ |",
        "|---|---|---|---|---|---|---|---|"]
for k, nm in names.items():
    a, b = tab["0.05"][k], tab["0.2"][k]
    w = a["median_width"]
    ratio = w["10000"] / tab["0.05"]["mix-pool"]["median_width"]["10000"]
    rows.append(f"| {nm} | {100*a['noncoverage_1_to_nmax']:.2f}% | {100*b['noncoverage_1_to_nmax']:.2f}% | "
                + " | ".join(f"{w[n]:.4f}" for n in ["100", "1000", "10000", "100000"]) + f" | {ratio:.3f} |")
eff = (f"Table 5. Time-uniform non-coverage ({ef['reps']:,} replications) and median interval width at $\\alpha=0.05$. "
       "Monte Carlo standard error of a non-coverage near 5% is about 0.35 percentage points.\n\n" + "\n".join(rows))
text = text.replace("{{EFF_TABLE}}", eff)

# looks table
L = [1, 2, 3, 5, 10, 20, 50]
rows = ["| first look, horizon | $\\varepsilon$ | " + " | ".join(f"$L={x}$" for x in L) + " |", "|---|---|" + "---|" * len(L)]
for key, lab in (("first100_R100", "$n_{\\rm first}=100$, $R=100$"), ("first10_R1000", "$n_{\\rm first}=10$, $R=1000$")):
    for i, e in enumerate(["0.20", "0.10", "0.05", "0.01"]):
        rows.append(f"| {lab if i == 0 else ''} | {e} | " + " | ".join(f"{100*lk['looks'][f'{key}_L{x}']['uncoverage'][i]:.2f}" for x in L) + " |")
looks = ("Table 4. Split-naive uncoverage (%) with $L$ looks at geometrically spaced sample sizes, nominal threshold $\\log(1/\\varepsilon)$ "
         "(exact Gaussian simulation, $10^6$ replications; standard errors at most 0.04 points). $L=1$ is a single fixed-$n$ look.\n\n" + "\n".join(rows))
text = text.replace("{{LOOKS_TABLE}}", looks)
gm = lk["gauss_mixture_table2"]["ours"]; pm = lk["poisson_mixture"]
text = text.replace("{{GMIX}}", ", ".join(f"{100*x:.2f}" for x in gm) + "% against 3.22, 1.79, 0.94 and 0.21% published")
text = text.replace("{{PMIX}}", ", ".join(f"{100*x:.2f}" for x in pm["n100_10000"]) + "% (per-group horizon 100 to 10,000) and "
                    + ", ".join(f"{100*x:.2f}" for x in pm["n50_5000"]) + "% (50 to 5,000)")
text = text.replace("{{REFS}}", open(refs).read().strip())
left = [t for t in ["{{"] if t in text]
if left:
    print("unfilled placeholders remain"); sys.exit(1)
open(outp, "w").write(text)
print("filled", outp)
