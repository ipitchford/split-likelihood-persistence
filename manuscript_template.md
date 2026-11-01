---
title: "The recomputed split likelihood-ratio interval is not a confidence sequence"
subtitle: "A conjecture of Pace, Salvan and Sartori resolved, with sharp constants, repairs and finite-horizon severity"
author: "Anonymous"
date: "1 November 2026 · Version 0.1.0-candidate · [doi:10.5281/zenodo.23856728](https://doi.org/10.5281/zenodo.23856728) · Evidence Press research release; not externally peer reviewed"
abstract: |
  Universal inference (Wasserman, Ramdas and Balakrishnan, 2020) turns a split likelihood ratio into a confidence set that is valid at any fixed sample size. Pace, Salvan and Sartori (2023, Section 3) conjectured that the same set, recomputed after every new pair of observations with the same threshold, has asymptotic persistence at least $1-\varepsilon$: with probability at least $1-\varepsilon$ the true parameter would eventually stay inside every later set. We show that the conjecture is false. A zero--one law argument gives asymptotic persistence $0$ in every model, for any permutation-symmetric training estimator, whenever the fixed-sample non-coverage does not vanish (Proposition 0). For the Gaussian two-stream scheme the log split likelihood ratio at the truth is an exact quadratic form in two independent random walks, and its $\limsup$ divided by $\log\log n$ is $(\sqrt5-1)/2$ almost surely (Theorem 1). The same constant holds in every regular model that satisfies the information identity (Theorem 2). Unequal splits give $(\sqrt{r^2+4r}-r)/2$, and fixed-$K$ cross-fitting gives $c(1/(K-1))$ for arithmetic averaging and $1/K$ for equal-weight geometric, harmonic or minimum aggregation (Theorems 3 and 4). Replacing the threshold $\log(1/\varepsilon)$ by $\lambda_1g(n)^2/n$, with $\lambda_1=(\sqrt5-1)/4$ and $g$ any valid boundary for a scalar Gaussian random walk, gives a genuine confidence sequence (Theorem 5). In the iterated-logarithm scale this repair is sharp: thresholds $\lambda_1\{2\log\log n+\kappa\log\log\log n+O(1)\}$ are achievable if and only if $\kappa>3$ (Theorem 6). A checkpoint repair works in every model. An exact Ornstein--Uhlenbeck representation measures the practical damage, which is real but slow. Persistence decays like a power of the monitoring horizon with exponent $0.042$ at nominal level $0.80$, $0.0054$ at $0.95$ and $0.0005$ at $0.99$. The nominal $95\%$ error is first exceeded only at horizons of roughly $10^4$ times the first look, and the $99\%$ error only at about $10^{9}$. We reproduce Pace et al.'s Gaussian tables to Monte Carlo accuracy. The classical tools are credited where used. What is new is their application to this statistic: the constants, the aggregation classification, the transfer of the sharp threshold and the severity calculation.
---

# At a glance: what this paper is about

**The problem in plain English.** Many statistical methods come with a promise such as "this interval contains the true value 95% of the time". Usually that promise holds only if you look at the data *once*, at a sample size fixed in advance. Scientists rarely work that way: data arrive over time and people check the interval repeatedly. Methods built for repeated checking are called *confidence sequences*. Their promise is stronger: with probability at least 95%, the interval contains the truth *at every sample size at once*.

**The specific question.** *Universal inference* is a popular recipe that gives a valid single-look interval in almost any statistical model. You split one sample into two halves, fit the model on one half, and test candidate parameter values on the other. (The two halves are a device of the method; they are not treatment groups.) Pace, Salvan and Sartori (2023) asked whether that interval, recomputed every time a new pair of observations arrives, behaves like a confidence sequence in the long run. They conjectured that it does.

**Our answer.** No. If you keep checking forever, the interval is certain to exclude the true value again and again. The long-run guarantee is not 95% but zero (Proposition 0 and Theorem 1). Changing the split fraction or averaging over several splits does not change that (Theorems 3 and 4). A simple change of threshold repairs the method, and that repair cannot be improved in its growth rate (Theorems 5 and 6).

**How much does it matter in practice?** Less than "zero" suggests, which is why simulations did not reveal it. The chance of staying correct throughout decays only like a small power of how long you keep monitoring. The table gives the Gaussian continuum limit (Section 7); the 99% crossing horizon is an extrapolation.

{{GLANCE}}

("105x" means a last look 105 times later than the first, not $10^5$.) So the problem is serious for loose levels and very long monitoring, and mild for 95--99% intervals over horizons below about $10^3$--$10^4$ times the first look. A handful of interim looks costs much less than continuous monitoring (Table 4 in Section 7.4). Section 7.3 gives numerically calibrated thresholds for a horizon fixed in advance.

**How to read this paper.** Section 2 fixes notation, gives the simple zero--one-law refutation and states one probability lemma. Section 3 gives the exact Gaussian rate. Sections 4 and 5 show that the rate does not depend on the Gaussian model, the split fraction or the aggregation rule. Section 6 gives repairs and their sharpness. Section 7 measures the finite-horizon severity. Section 8 compares efficiency with existing anytime-valid intervals. Section 9 says what is and is not affected. Proofs are complete and appear with the results; the scripts are listed in Section 10.

**What is not claimed.** Single-look universal inference remains valid. So do the genuinely sequential "split-exact" (running maximum-likelihood) construction of Wasserman et al. (2020), also used by Pace et al., and Pace et al.'s mixture intervals. Their separate Section 6 conjecture is not addressed. The core mechanism, a law of the iterated logarithm, is classical; Section 1.2 says exactly what is new.

# 1. Introduction

## 1.1 The conjecture

Let $Y_1,Y_2,\dots$ be independent with density $p_\theta$, $\theta\in\Theta\subseteq\mathbb R^d$. Universal inference (Wasserman, Ramdas and Balakrishnan, 2020) splits a sample into a *validation* part $D_0$ and a *training* part $D_1$, computes an estimator $\hat\theta_1$ from $D_1$, and forms
$$
E(\theta)=\prod_{i\in D_0}\frac{p_{\hat\theta_1}(Y_i)}{p_\theta(Y_i)},\qquad
C=\{\theta:\ E(\theta)<1/\varepsilon\}.
$$
Conditionally on $D_1$, $E(\theta)$ has expectation at most one under $p_\theta$. So $E(\theta)$ is an *e-value*, and by Markov's inequality $C$ covers the truth with probability at least $1-\varepsilon$ at that one sample size.

Pace, Salvan and Sartori (2023) study confidence sequences with persistence level $1-\varepsilon$, meaning $P_\theta(\theta\in C_n\text{ for every }n\ge1)\ge1-\varepsilon$. In their Section 3 they take two groups of size $n$ that each grow by one fresh observation per step and never mix, and recompute the split set (their equation (6)) at every $n$ with the same threshold. They call the result *split-naive*. They note that it "is not proposed as a confidence sequence" by Wasserman et al. They then conjecture that it may be used as an approximate confidence sequence with *asymptotic persistence level* $1-\varepsilon$:
$$
\lim_{n\to\infty}P_\theta\left(\theta\in C_m\ \text{for every}\ m\ge n\right)\ \ge\ 1-\varepsilon. \tag{1.1}
$$
They single out a finite-horizon version as the property of practical interest: for $n_{\max}\gg n_{\min}$, coverage at every $n_{\min}\le n\le n_{\max}$ with probability at least $1-\varepsilon$. Their Tables 2 and 3 report Monte Carlo uncoverage over $100\le n\le 40{,}000$ and $400\le n\le 80{,}000$. At $\varepsilon=0.20$ it is $23.66\%$ and $21.67\%$, and at $\varepsilon=0.10,0.05,0.01$ it is below $\varepsilon$. Their Poisson example (Table 5) reports larger uncoverages; we return to it in Section 7.5. (Quotations and table values were checked against the publisher's open-access full text, CC BY 4.0, Early View version of 9 December 2022.)

## 1.2 Contributions and what is classical

1. **Refutation (Proposition 0).** In any model, with any permutation-symmetric training estimator, the split-naive sequence has asymptotic persistence $0$ whenever its single-time non-coverage does not tend to zero. The argument uses only the Hewitt--Savage zero--one law. It covers the irregular models for which universal inference was designed.
2. **Exact rate (Theorem 1).** In the Gaussian two-stream scheme, $\log E_n(\theta)=(S_nT_n-T_n^2/2)/n$ exactly, where $S_n,T_n$ are independent standardised random walks. Hence $\limsup_n\log E_n(\theta)/\log\log n=(\sqrt5-1)/2$ almost surely. The tail intersection $\bigcap_{m\ge n}C_m$ is empty almost surely, so Pace et al.'s compatibility criterion fails as well.
3. **Universality of the rate (Theorem 2).** The same constant holds in every regular parametric model of fixed dimension that satisfies the information identity, including Pace et al.'s Poisson example. For a working likelihood the constant is $\tfrac{\sqrt5-1}{2}\mu_{\max}(H^{-1}J)$.
4. **Splits and aggregation (Theorems 3 and 4).** A fixed validation-to-training ratio $r$ gives the constant $c(r)=(\sqrt{r^2+4r}-r)/2$. For fixed-$K$ cross-fitting we give exact constants for weighted power means and bounds for symmetric e-merging functions.
5. **Repairs (Theorems 5 and 6, Proposition 5A).** Any valid boundary $g$ for a scalar Gaussian walk yields a valid confidence sequence through $b(n)=\lambda_1g(n)^2/n$. Within thresholds $\lambda_1\{2\log\log n+\kappa\log\log\log n+O(1)\}$, validity is achievable if and only if $\kappa>3$. A checkpoint repair is valid in every model.
6. **Finite-horizon severity (Proposition 7, Section 7).** An exact Ornstein--Uhlenbeck (OU) representation gives the continuum limit of Pace et al.'s experiment. From it we compute failure exponents, the horizons at which the error exceeds nominal, horizon-calibrated thresholds and the effect of a finite number of looks. We also reproduce their Gaussian tables.
7. **Efficiency (Section 8).** At practical sample sizes the repaired split interval is about 13--14% wider than the pooled normal-mixture interval and about 10--11% wider than the running-MLE interval.

*Classical ingredients, credited where used.* We rely on the following standard results:

- the scalar law of the iterated logarithm (LIL) for i.i.d. finite-variance sums (Hartman and Wintner, 1941) and its functional form (Strassen, 1964);
- the Hewitt--Savage zero--one law (Hewitt and Savage, 1955);
- Ville's (1939) maximal inequality;
- Robbins's normal mixtures (Robbins, 1970; Darling and Robbins, 1967);
- the geometric-epoch "stitching" of Howard, Ramdas, McAuliffe and Sekhon (2021, Theorem 1);
- the Kolmogorov--Erdős--Feller integral test (Erdős, 1942; Feller, 1943; Motoo, 1959);
- Vovk and Wang's (2021) characterisation of symmetric e-merging functions;
- the time-mixture construction of Ramdas and Wang (2025, Section 7.9).

Two further points of credit. First, a second-order term of $\tfrac32\log\log\log$ (on the scale of the boundary itself) is classical for sequential boundaries (Darling and Robbins, 1967; Robbins, 1970), and Balsubramani (2026) recently re-derived it from a minimax viewpoint. Second, the inflation of error under repeated significance testing is the subject of group-sequential analysis (Armitage, McPherson and Rowe, 1969; Pocock, 1977); Section 7 is its counterpart for this statistic.

Fixed-sample properties of the split likelihood-ratio test (power, the choice of split ratio, Gaussian analysis) are studied by Tse and Davison (2022), Dunn, Ramdas, Balakrishnan and Wasserman (2023) and Strieder and Drton (2022). Boundary-crossing theory for confidence sequences goes back to Robbins and Siegmund (1970) and Lai (1976), and asymptotic confidence sequences are developed by Waudby-Smith et al. (2024). Ramdas and Wang (2025, Section 7.7) remark, without proof, that split likelihood-ratio e-variables are not e-processes. Our results make that remark exact for this statistic and show its consequence for (1.1).

*Prior-art search.* On 31 October and 1 November 2026 we searched Semantic Scholar, arXiv, OpenAlex and general web indexes for the refutation, the constants and the aggregation results; the query record is kept with the release. We also checked all works recorded as citing Pace et al. (2023) and Pace and Salvan (2020), including Di Caterina, Pace, Salvan and Sartori (2026). We found no earlier statement of the refutation, of the constants $(\sqrt5-1)/2$, $c(r)$ and $1/K$, or of the failure exponents. A bounded search does not prove priority.

# 2. Setting, a zero--one law, and one lemma

**Two-stream scheme.** Let $(X_i)$ (validation) and $(Y_i)$ (training) be independent i.i.d. sequences from $p_\theta$, observed in pairs. Pace et al.'s random once-and-for-all assignment of each new pair to the two groups gives exactly this. After $n$ pairs let $\hat\theta_n=t_n(Y_1,\dots,Y_n)$, and define
$$
E_n(\theta)=\prod_{i=1}^n\frac{p_{\hat\theta_n}(X_i)}{p_\theta(X_i)},\qquad
C_n=\{\theta:\ \log E_n(\theta)<\log(1/\varepsilon)\}.
$$
This is Pace et al.'s split-naive set. For each fixed $n$, $E_\theta[E_n(\theta)\mid Y_1,\dots,Y_n]\le1$ (with equality when the densities have common support), so $P_\theta(\theta\notin C_n)\le\varepsilon$. Write $L_2(n)=\log\log n$ and $L_3(n)=\log\log\log n$.

**Proposition 0 (zero--one law).** Suppose each $t_n$ is a symmetric function of its $n$ arguments, and let $A_n=\{\theta\notin C_n\}$. If $\limsup_nP_\theta(A_n)>0$, then $P_\theta(A_n\text{ infinitely often})=1$. Consequently $P_\theta(\theta\in C_m\ \forall m\ge n)=0$ for every $n$, so the left side of (1.1) is $0$, for any model and any $\varepsilon\in(0,1)$.

*Proof.* By the reverse Fatou inequality, $P(\limsup A_n)=\lim_nP(\bigcup_{m\ge n}A_m)\ge\limsup_nP(A_n)>0$. The pairs $\xi_i=(X_i,Y_i)$ are i.i.d. Let $\pi$ permute $(\xi_1,\dots,\xi_k)$. For every $n\ge k$, the product over $X_1,\dots,X_n$ and the symmetric estimator $t_n(Y_1,\dots,Y_n)$ are unchanged, so $A_n$ is unchanged. Hence $\limsup A_n$ is invariant under finite permutations. By the Hewitt--Savage zero--one law (Hewitt and Savage, 1955) its probability is $0$ or $1$, and so it is $1$. $\square$

The hypothesis is mild. In the Gaussian model $Q_n=\log E_n(\theta)$ has the same law for every $n$ (Theorem 1(i)), so $P(A_n)=P(Z_1Z_2-Z_2^2/2\ge\log(1/\varepsilon))>0$ for independent standard normals $Z_1,Z_2$; at $\varepsilon=0.20$ this is $1.13\%$. Under the conditions of Theorem 2 the limiting non-coverage is also positive; in dimension $d$ under (R4) the limit law of $Q_n$ is $\lambda_1\chi^2_d+\lambda_2\tilde\chi^2_d$ with independent chi-squares (for $d=2$ the non-coverage at $\varepsilon=0.20$ is about $2.0\%$).

Proposition 0 settles the conjecture but says nothing about *how often* or *how fast* the exclusions recur. Those rates determine both the repair and the practical severity, and the rest of the paper is about them. (Proposition 0 was suggested during our internal review. We include it because it is the shortest complete refutation.)

**Lemma 1 (vector LIL for continuous 2-homogeneous functionals).** Let $Z_1,Z_2,\dots$ be i.i.d. in $\mathbb R^m$ with mean zero and nonsingular covariance $\Sigma$, and let $S_n=\sum_{i\le n}Z_i$. Let $\mathcal E=\{\Sigma^{1/2}y:\|y\|\le1\}$, and let $f:\mathbb R^m\to\mathbb R$ be continuous with $f(tx)=t^2f(x)$ for $t>0$. Then almost surely
$$
\limsup_{n\to\infty}\frac{f(S_n)}{n\log\log n}\ \le\ 2\max_{x\in\mathcal E}f(x),
$$
with equality whenever $\max_{\mathcal E}f>0$. In particular, for symmetric $A$ with $\lambda_{\max}(\Sigma^{1/2}A\Sigma^{1/2})>0$, $\limsup_n S_n^\top AS_n/(n\log\log n)=2\lambda_{\max}(\Sigma^{1/2}A\Sigma^{1/2})$.

*Proof.* Put $x_n=S_n/\sqrt{2n\log\log n}$ for $n\ge3$. By homogeneity $f(S_n)/(n\log\log n)=2f(x_n)$. Fix a vector $u$. The scalar Hartman--Wintner LIL (Hartman and Wintner, 1941) applied to $\langle u,Z_i\rangle$ gives $\limsup_n\langle u,x_n\rangle=\sqrt{u^\top\Sigma u}=h(u)$ almost surely, where $h$ is the support function of $\mathcal E$.

(a) *Eventual closeness to $\mathcal E$.* Applying this to the coordinate vectors gives $\sup_n\|x_n\|<\infty$ a.s. Applying it to a countable dense set $\mathcal U$ of unit vectors gives, on one event of probability one, $\limsup\langle u,x_n\rangle\le h(u)$ for every $u\in\mathcal U$. Suppose that along a subsequence $\operatorname{dist}(x_{n_j},\mathcal E)\ge\eta>0$. By boundedness we may assume $x_{n_j}\to\bar x$, so that $\operatorname{dist}(\bar x,\mathcal E)\ge\eta$. Since $\mathcal E$ is compact and convex, some unit vector $u$ has $\langle u,\bar x\rangle\ge h(u)+\eta$. By continuity, some $u'\in\mathcal U$ has $\langle u',\bar x\rangle>h(u')+\eta/2$, a contradiction. Hence $\operatorname{dist}(x_n,\mathcal E)\to0$ a.s. Uniform continuity of $f$ on a compact neighbourhood of $\mathcal E$ then gives $\limsup f(x_n)\le\max_{\mathcal E}f$.

(b) *Reaching a maximiser.* Let $\max_{\mathcal E}f>0$ be attained at $x^*$. Then $x^*\in\partial\mathcal E$: otherwise $tx^*\in\mathcal E$ and $f(tx^*)=t^2f(x^*)>f(x^*)$ for some $t>1$. Write $x^*=\Sigma^{1/2}y^*$ with $\|y^*\|=1$, and put $u=\Sigma^{-1/2}y^*$. Then $h(u)=1$, and $x^*$ is the unique maximiser of $\langle u,\cdot\rangle$ on $\mathcal E$, because $\langle u,\Sigma^{1/2}y\rangle=\langle y^*,y\rangle$ attains its maximum over the unit ball only at $y=y^*$. By Hartman--Wintner, choose a subsequence with $\langle u,x_{n_j}\rangle\to1$, and then a further subsequence converging to some $\bar x$. By (a), $\bar x\in\mathcal E$ and $\langle u,\bar x\rangle=1$, so $\bar x=x^*$ and $f(x_{n_j})\to f(x^*)$. For the quadratic case, $\max_{\mathcal E}x^\top Ax=\max_{\|y\|\le1}y^\top\Sigma^{1/2}A\Sigma^{1/2}y$. $\square$

*Remark.* Lemma 1 is a standard consequence of the multivariate cluster-set form of the LIL (Strassen, 1964; for $\mathbb R^m$ see, for example, Berkes and Philipp, 1979; Kuelbs, 1976). We derive it from the scalar theorem to keep the paper self-contained. Applied to $-f$, it gives $\liminf f(S_n)/(n\log\log n)=2\min_{\mathcal E}f$ when $\min_{\mathcal E}f<0$.

# 3. The Gaussian case: the exact rate

**Theorem 1.** Let $X_i,Y_i$ be independent $N(\theta,\sigma^2)$ with $\sigma$ known, and let $\hat\theta_n=\bar Y_n$. Put $S_n=\sum_{i\le n}(X_i-\theta)/\sigma$ and $T_n=\sum_{i\le n}(Y_i-\theta)/\sigma$.

(i) For every $n$, $Q_n:=\log E_n(\theta)=(S_nT_n-T_n^2/2)/n=(S_n,T_n)A(S_n,T_n)^\top/n$, where $A=\begin{pmatrix}0&1/2\\1/2&-1/2\end{pmatrix}$ has eigenvalues $\lambda_1=(\sqrt5-1)/4$ and $\lambda_2=-(\sqrt5+1)/4$.

(ii) Almost surely, $\limsup_nQ_n/\log\log n=2\lambda_1=(\sqrt5-1)/2\approx0.618$ and $\liminf_nQ_n/\log\log n=2\lambda_2=-(\sqrt5+1)/2$.

(iii) For every $\varepsilon\in(0,1)$ and every $n_0$, $P_\theta(\theta\in C_n\ \forall n\ge n_0)=0$.

(iv) Almost surely $\bigcap_{n\ge n_0}C_n=\varnothing$ for every $n_0$, and $P_\theta\big(\bigcap_{n_0\le n\le N}C_n=\varnothing\big)\to1$ as $N\to\infty$.

*Proof.* (i) With $d_n=\bar Y_n-\theta=\sigma T_n/n$,
$$\log E_n(\theta)=\sum_{i\le n}\frac{(X_i-\theta)^2-(X_i-\theta-d_n)^2}{2\sigma^2}=\frac{\sigma S_nd_n-nd_n^2/2}{\sigma^2}=\frac{S_nT_n-T_n^2/2}{n}.$$
The characteristic polynomial of $A$ is $\lambda^2+\lambda/2-1/4$.

(ii) The vectors $((X_i-\theta)/\sigma,(Y_i-\theta)/\sigma)$ are i.i.d. $N(0,I_2)$. Apply Lemma 1 to $f(x)=x^\top Ax$ and to $-f$.

(iii) By (ii), $Q_n>\log(1/\varepsilon)$ for infinitely many $n$ a.s., and $\theta\notin C_n$ exactly when $Q_n\ge\log(1/\varepsilon)$.

(iv) Completing the square gives
$C_n=\big\{\vartheta:|\vartheta-\bar X_n|<\sqrt{(\bar X_n-\bar Y_n)^2+2\sigma^2\log(1/\varepsilon)/n}\big\}$, which is Pace et al.'s equation (9).
The centre tends to $\theta$ and the half-width to $0$ almost surely. By (ii) there is an $m\ge n_0$ with $Q_m>\log(1/\varepsilon)+1$. Then $\theta$ is at some positive distance $\delta$ from $\overline{C_m}$, since $\vartheta\mapsto\log E_m(\vartheta)$ is a continuous convex quadratic. For all large $n$, $C_n\subset(\theta-\delta/2,\theta+\delta/2)$, which is disjoint from $C_m$. So a finite pair $(m,n)$ already has empty intersection. The events $\{\bigcap_{n_0\le n\le N}C_n=\varnothing\}$ increase in $N$, and their union has probability one. $\square$

*Interpretation.* At each fixed $n$, $E_n=e^{Q_n}$ is an e-value. As a process, though, $Q_n$ is a quadratic form in a two-dimensional random walk scaled by $1/n$, and it returns to heights of order $\log\log n$ forever. A fixed threshold cannot contain it.

# 4. The rate in every regular model

Let $\ell(\vartheta;x)=\log p_\vartheta(x)$ be a working log-likelihood on an open set $\Theta\subseteq\mathbb R^d$. Let the data be i.i.d. from $P$, and let $\theta_0$ be a point at which the following hold:

* **(R1)** On a closed ball $B$ around $\theta_0$, $\ell(\cdot;x)$ is $C^3$ for every $x$, with all third partial derivatives bounded by $M(x)$, where $E_PM<\infty$.
* **(R2)** $E_P\nabla\ell(\theta_0;Y)=0$; $J=\operatorname{Cov}_P\nabla\ell(\theta_0;Y)$ and $H=-E_P\nabla^2\ell(\theta_0;Y)$ are positive definite; and $E_P\|\nabla^2\ell(\theta_0;Y)\|^2<\infty$.
* **(R3)** The training estimator satisfies $\hat\theta_n\to\theta_0$ a.s. and, for all large $n$, $\sum_{i\le n}\nabla\ell(\hat\theta_n;Y_i)=0$.
* **(R4)** (correct specification) $P=p_{\theta_0}$ and $J=H$, the Fisher information.

(R4) is needed for the universal constant. (R1)--(R3) do not imply the information identity $J=H$. Our AI-assisted adversarial reviewer (Astra) supplied an explicit counterexample: the mixture $(1-\theta^2)\phi(x-\theta)+\theta^2\phi(x-1/\theta)$ at $\theta_0=0$, which has $J=1$ and $H=3$. In standard regular models, $J=H$ follows from differentiating $\int p_\vartheta=1$ twice under the integral sign.

For a full regular exponential family with $\theta_0$ interior, (R1)--(R4) hold with $\hat\theta_n$ the MLE. The MLE exists for all large $n$ a.s., because the mean sufficient statistic converges to an interior point. This covers Pace et al.'s Poisson example with $\lambda>0$; their replacement of $\hat\lambda$ by $0.5/n$ when the training sum is zero happens only finitely often. It also covers the Bernoulli model with $0<p<1$.

**Theorem 2.** Under (R1)--(R3), with $U_n=\sum_{i\le n}\nabla\ell(\theta_0;X_i)$ and $V_n=\sum_{i\le n}\nabla\ell(\theta_0;Y_i)$, almost surely
$$
\log E_n(\theta_0)=\frac{U_n^\top H^{-1}V_n-\tfrac12V_n^\top H^{-1}V_n}{n}+O\!\left(\frac{(\log\log n)^{3/2}}{\sqrt n}\right),
\qquad
\limsup_{n\to\infty}\frac{\log E_n(\theta_0)}{\log\log n}=\frac{\sqrt5-1}{2}\,\mu_{\max}(H^{-1}J).
$$
Under (R4) the constant is $(\sqrt5-1)/2$ for every fixed dimension $d$, and the split-naive sequence has asymptotic persistence $0$ for every $\varepsilon\in(0,1)$.

*Proof.* All $O(\cdot)$ and $o(\cdot)$ statements below hold almost surely.

*Preliminaries.* By Hartman--Wintner coordinatewise, $\|U_n\|,\|V_n\|=O(\sqrt{n\log\log n})$. Put $G_n=\sum_{i\le n}\nabla^2\ell(\theta_0;Y_i)$ and $G'_n=\sum_{i\le n}\nabla^2\ell(\theta_0;X_i)$. By (R2) and Hartman--Wintner entrywise, $G_n,G'_n=-nH+O(\sqrt{n\log\log n})$. By the strong law, $\sum_{i\le n}M(Y_i)$ and $\sum_iM(X_i)$ are $O(n)$.

*Step 1 (the estimator).* Let $h_n=\hat\theta_n-\theta_0\to0$; for large $n$, $\hat\theta_n\in B$. Apply Taylor's theorem with integral remainder to each component of the training score. This gives $0=V_n+G_nh_n+R_n$ with $\|R_n\|\le c_d\big(\sum_iM(Y_i)\big)\|h_n\|^2$, where $c_d$ depends only on $d$. Since $G_n/n\to-H$ and $h_n\to0$, we get $(H+o(1))h_n=V_n/n$, so $\|h_n\|=O(\sqrt{\log\log n/n})$. Writing $nHh_n=V_n+R_n+(G_n+nH)h_n$,
$$
h_n=H^{-1}V_n/n+e_n,\qquad\|e_n\|\le\|H^{-1}\|\big(\|R_n\|+\|G_n+nH\|\|h_n\|\big)/n=O(\log\log n/n).
$$

*Step 2 (the validation log-likelihood).* Expanding around $\theta_0$,
$$\log E_n(\theta_0)=h_n^\top U_n+\tfrac12h_n^\top G'_nh_n+R'_n,\qquad |R'_n|\le c_d\textstyle\sum_iM(X_i)\|h_n\|^3=O\big(n(\log\log n/n)^{3/2}\big).$$
By Step 1:

- $h_n^\top U_n=V_n^\top H^{-1}U_n/n+O\big((\log\log n/n)\sqrt{n\log\log n}\big)$;
- $\tfrac12h_n^\top G'_nh_n=-\tfrac n2h_n^\top Hh_n+O(\sqrt{n\log\log n}\,\log\log n/n)$;
- $-\tfrac n2h_n^\top Hh_n=-V_n^\top H^{-1}V_n/(2n)+O\big(\|V_n\|\|e_n\|+n\|e_n\|^2\big)$.

Every error term is $O((\log\log n)^{3/2}/\sqrt n)$.

*Step 3 (the constant).* The vectors $(J^{-1/2}\nabla\ell(\theta_0;X_i),J^{-1/2}\nabla\ell(\theta_0;Y_i))$ are i.i.d. with mean $0$ and covariance $I_{2d}$. In these coordinates the quadratic form has matrix $A\otimes M$, where $M=J^{1/2}H^{-1}J^{1/2}$ is positive definite. The eigenvalues of $A\otimes M$ are the products $\lambda_i(A)\mu_j(M)$, so the largest is $\lambda_1\mu_{\max}(M)=\lambda_1\mu_{\max}(H^{-1}J)$. Lemma 1 gives the constant; the vanishing remainder does not affect a limit normalised by $\log\log n$. Under (R4), $M=I_d$ and $E_n(\theta_0)$ is a valid e-value, so persistence $0$ follows as in Theorem 1(iii). $\square$

*Remarks.*

(a) Composite likelihoods, the main subject of Pace et al., are covered by the working-likelihood form. With a composite likelihood in place of $\ell$, $J\ne H$ in general, and the fixed-$n$ e-value property itself may fail (Nguyen, 2020). The constant $\tfrac{\sqrt5-1}{2}\mu_{\max}(H^{-1}J)$ then describes the growth of a statistic that need not be an e-value.

(b) Strieder and Drton (2022, Theorem 3.1) give the fixed-$n$ limiting law of the split likelihood ratio. Theorem 2 is its pathwise, sequential counterpart.

(c) A superefficient training estimator can make $E_n(\theta_0)=1$ eventually at a single point $\theta_0$. Proposition 0 still applies at every other parameter value where the fixed-$n$ non-coverage does not vanish.

(d) Boundary parameters, growing dimension and estimators without the stated expansion are excluded from Theorem 2, though not from Proposition 0.

# 5. Split fractions and aggregation

**Theorem 3 (unequal splits).** Suppose each update adds $a$ validation and $b$ training observations, $a,b\ge1$, and let $r=a/b$. In the Gaussian model, and in every model satisfying (R1)--(R4), almost surely
$$
\limsup_{n\to\infty}\frac{\log E_n(\theta_0)}{\log\log n}=c(r):=\frac{\sqrt{r^2+4r}-r}{2}\in(0,1).
$$
The function $c$ is increasing, with $c(1)=(\sqrt5-1)/2$, $c(r)\sim\sqrt r$ as $r\to0$, and $c(r)\to1$ as $r\to\infty$. No fixed split fraction restores (1.1).

*Proof.* Gaussian case. Let $S_n,T_n$ be the standardised sums of the $na$ validation and $nb$ training observations, and put $s_n=S_n/\sqrt a$ and $t_n=T_n/\sqrt b$. Then $(s_n,t_n)$ is a sum of $n$ i.i.d. $N(0,I_2)$ vectors (the whitened block sums). The computation of Theorem 1(i), with $d_n=\sigma T_n/(nb)$, gives
$$\log E_n=\frac{S_nT_n}{nb}-\frac{na\,T_n^2}{2(nb)^2}=\frac{\sqrt r\,s_nt_n-\tfrac r2t_n^2}{n}.$$
This is the quadratic form of $A_r=\begin{pmatrix}0&\sqrt r/2\\\sqrt r/2&-r/2\end{pmatrix}$, divided by $n$. The characteristic polynomial of $A_r$ is $\lambda^2+\tfrac r2\lambda-\tfrac r4$, so $\lambda_{\max}=(\sqrt{r^2+4r}-r)/4$, and Lemma 1 applies. For regular models, Steps 1--3 of Theorem 2 applied to the block sums give the whitened matrix $A_r\otimes I_d$. The properties of $c$ follow from $c(r)=2r/(\sqrt{r^2+4r}+r)$. $\square$

For example, $c(1/9)=0.282$, $c(1/3)=0.434$, $c(3)=0.791$ and $c(9)=0.908$. Spending less of the data on validation slows the growth but never stops it.

**$K$-fold cross-fitting.** Take $K\ge2$ independent streams of $n$ observations each. Fold $j$ validates stream $j$ at the mean of the other $K-1$ streams, giving fold e-values $E^{(1)}_n,\dots,E^{(K)}_n$. An *e-merging function* (Vovk and Wang, 2021) is a map $F$ such that $F(E^{(1)},\dots,E^{(K)})$ is an e-value whenever the inputs are. For positive weights $w_j$ summing to one, the weighted power mean $M_{p,w}=(\sum_jw_jE_j^p)^{1/p}$ is an e-merging function for every $p\le1$, because it is at most the arithmetic mean. Here $p=0$ means the weighted geometric mean and $p=-\infty$ the minimum.

**Theorem 4 (aggregation).** In the Gaussian $K$-stream scheme, let
$$\Gamma(F)=\limsup_{n\to\infty}\frac{\log F(E^{(1)}_n,\dots,E^{(K)}_n)}{\log\log n},$$
and let $A_j$ be the matrix of fold $j$ defined in the proof. Almost surely:

(i) for $0<p\le1$, $\Gamma(M_{p,w})=c(1/(K-1))$, the constant of a single fold;

(ii) for $p=0$, $\Gamma(M_{0,w})=2\lambda_{\max}(\sum_jw_jA_j)\in[1/K,\,c(1/(K-1))]$, with $\Gamma=1/K$ for equal weights;

(iii) for $-\infty\le p<0$ and any positive weights, $\Gamma(M_{p,w})=1/K$;

(iv) every symmetric e-merging function $F$ has $\Gamma(F)\le c(1/(K-1))$, and every e-merging function with $F\ge\gamma\min_jE_j$ for some constant $\gamma>0$ has $\Gamma(F)\ge1/K$.

Hence every weighted power mean with $p\le1$, and every e-merging function that dominates a positive multiple of the minimum, yields asymptotic persistence $0$. Symmetry alone is not enough: the trivial e-merging function $F\equiv1$ has persistence $1$ and no power. For $K=2$ the constants are $c(1)=0.618$ (arithmetic) and $1/2$ (geometric); for $K=5$ they are $0.390$ and $0.2$.

*Proof.* Let $\mathbf S_n\in\mathbb R^K$ collect the standardised stream sums; it is a sum of $n$ i.i.d. $N(0,I_K)$ vectors. As in Theorem 1(i), $\log E^{(j)}_n=\mathbf S_n^\top A_j\mathbf S_n/n$, where $A_j$ is the matrix of the form $s_j\bar s_{-j}-\bar s_{-j}^2/2$ with $\bar s_{-j}=(K-1)^{-1}\sum_{i\ne j}s_i$. In the orthonormal coordinates $s_j$ and $t=(K-1)^{-1/2}\sum_{i\ne j}s_i$, this form equals $\sqrt r\,s_jt-\tfrac r2t^2$ with $r=1/(K-1)$, and it vanishes on the orthogonal complement. Hence $\max_{\|x\|\le1}x^\top A_jx=c(r)/2$.

(i) For $p>0$, $w_{\min}^{1/p}\max_jE_j\le M_{p,w}\le\max_jE_j$, so $\log M_{p,w}=\max_j\log E^{(j)}_n+O(1)$. Apply Lemma 1 to $f=\max_jx^\top A_jx$, whose maximum over the unit ball is $c(r)/2$.

(ii) $\log M_{0,w}=\mathbf S_n^\top(\sum_jw_jA_j)\mathbf S_n/n$, and Lemma 1 gives the formula. The bounds follow from $\min_jE_j\le M_{0,w}\le\max_jE_j$, together with (iii) and (i). For equal weights a direct sum gives
$$\frac1K\sum_jA_j=A_K:=\frac{\mathbf 1\mathbf 1^\top}{2(K-1)^2}-\frac{(2K-1)I_K}{2K(K-1)^2},$$
whose eigenvalue on $\mathbf 1$ is $1/(2K)$ and on $\mathbf 1^\perp$ is $-(2K-1)/(2K(K-1)^2)$.

(iii) For $p<0$, $\min_jE_j\le M_{p,w}\le w_{\min}^{1/p}\min_jE_j$, so it suffices to treat $f=\min_jx^\top A_jx$. On the unit vector $\mathbf 1/\sqrt K$ the validation and training means coincide, so every fold equals $\tfrac1K(1-\tfrac12)=1/(2K)$. Also $f(x)\le x^\top A_Kx\le\|x\|^2/(2K)$. Hence $\max_{\|x\|\le1}f=1/(2K)$, and Lemma 1 applies.

(iv) By Vovk and Wang (2021, Theorem 3.2), every symmetric e-merging function satisfies $F\le\lambda+(1-\lambda)\bar E$ for some $\lambda\in[0,1]$, where $\bar E$ is the arithmetic mean. Hence $\log F\le\log2+\max(0,\log\bar E)$, and (i) gives the upper bound. The lower bound follows from (iii), because $\log F\ge\log\gamma+\log\min_jE_j$. $\square$

*Remark.* Theorem 4 assumes a fixed number of folds with a fixed assignment. It does not cover growing $K=K_n$, splits re-randomised as $n$ grows (common in practice), or leave-one-out averaging. Since $c(1/(K-1))\sim(K-1)^{-1/2}$ and $1/K$ both tend to zero, such schemes might behave differently; we claim nothing about them. Proposition 0 applies to any fixed assignment, and more generally whenever the exclusion events are invariant under finite permutations of the i.i.d. inputs; re-randomised splits with fresh randomisation at each $n$ are not covered.

# 6. Repairs

The failure lies in the *threshold*, not the statistic. The split statistic grows like $2\lambda_1\log\log n$, so the threshold must grow too. The first repair is sharp for the Gaussian scheme; the second is valid in every model.

**Theorem 5 (reduction).** In the Gaussian two-stream scheme, let $u$ be a unit eigenvector of $A$ for $\lambda_1$ and $v\perp u$ a unit vector. Set $W_n=u^\top(S_n,T_n)^\top$ and $\widetilde W_n=v^\top(S_n,T_n)^\top$. Then $W$ and $\widetilde W$ are independent standard Gaussian random walks, and pathwise
$$
Q_n=\frac{\lambda_1W_n^2+\lambda_2\widetilde W_n^2}{n}\ \le\ \frac{\lambda_1W_n^2}{n}.
$$
Consequently, if $g:\mathbb N\to(0,\infty)$ satisfies $P(\exists n\ge1:|W_n|\ge g(n))\le\alpha$, then with $b(n)=\lambda_1g(n)^2/n$ the sets
$$
C^{b}_n=\Big\{\vartheta:\ |\vartheta-\bar X_n|<\sqrt{(\bar X_n-\bar Y_n)^2+2\sigma^2b(n)/n}\Big\}
$$
satisfy $P_\theta(\exists n\ge1:\theta\notin C^b_n)\le\alpha$, so they form a $(1-\alpha)$ confidence sequence. The same holds, with the Gaussian working likelihood, for non-Gaussian streams whose standardised increment vector is jointly conditionally sub-Gaussian,
$$E[\exp(a\Delta S_n+b\Delta T_n)\mid\mathcal F_{n-1}]\le\exp\{(a^2+b^2)/2\}$$
(for example, independent streams with $1$-sub-Gaussian standardised increments), provided $g$ is valid for $1$-sub-Gaussian martingales.

*Proof.* $(u,v)$ is an orthonormal eigenbasis of $A$, so $x^\top Ax=\lambda_1(u^\top x)^2+\lambda_2(v^\top x)^2$. Orthogonal images of i.i.d. $N(0,I_2)$ increments are again i.i.d. $N(0,I_2)$, so $W$ and $\widetilde W$ are independent standard walks. Since $\lambda_2<0$ the inequality follows, and $\{\theta\notin C^b_n\}=\{Q_n\ge b(n)\}\subseteq\{|W_n|\ge g(n)\}$. The interval comes from completing the square, as in Theorem 1(iv). In the sub-Gaussian case, the stated condition with $(a,b)=su$ makes $W$ a $1$-sub-Gaussian martingale, and the identity for $Q_n$ is algebraic. $\square$

Two classical choices of $g$ are valid for every $n\ge1$:

* **Normal mixture** (Robbins, 1970). $(1+n/\rho)^{-1/2}\exp\{W_n^2/(2(n+\rho))\}$ is a nonnegative martingale with initial value $1$. Ville's inequality gives $g(n)^2=(n+\rho)\{2\log(1/\alpha)+\log(1+n/\rho)\}$, so $b(n)=\lambda_1(1+\rho/n)\{2\log(1/\alpha)+\log(1+n/\rho)\}$, which grows like $\lambda_1\log n$.
* **Stitching** (Howard et al., 2021, Theorem 1). Fix $\eta>1$ and $\delta>0$. For $\eta^k\le n<\eta^{k+1}$ put
  - $c_k=1+\log(k+1)$;
  - $w_k=\log(k+e)^{-\delta}-\log(k+1+e)^{-\delta}$, which is positive and telescopes to at most $1$;
  - $\rho_k=\eta^k/c_k$ and $L_k=\log\{1/(\alpha w_k)\}+\tfrac12\log(1+\eta c_k)$;
  - $g(n)^2=2n(1+1/c_k)L_k$.

  *Validity:* on epoch $k$, $n/\rho_k<\eta c_k$ and $\rho_k/n\le1/c_k$. So $|W_n|\ge g(n)$ forces the normal-mixture martingale with $\rho=\rho_k$ above $1/(\alpha w_k)$, and Ville's inequality charges at most $\alpha w_k$ to epoch $k$. *Growth:* as $n\to\infty$, $b(n)=\lambda_1\{2\log\log n+(3+2\delta)\log\log\log n\}+O(1)$. This follows from $\log(1/w_k)=\log k+(1+\delta)\log\log k+O(1)$, $\log k=\log\log n+O(1)$ and $(1+1/c_k)L_k=L_k+1+o(1)$. This is the iterated-logarithm form of Darling and Robbins (1967) and Robbins (1970); compare equation (35) of the arXiv version of Howard et al.

**Theorem 6 (sharpness).** In the Gaussian two-stream scheme let $b_{\kappa,C}(n)=\lambda_1\{2\log\log n+\kappa\log\log\log n+C\}$.

(i) If $\kappa\le3$, then for every $C\in\mathbb R$, almost surely $Q_n\ge b_{\kappa,C}(n)$ for infinitely many $n$. Therefore no threshold $b$ with $b(n)\le b_{3,C}(n)$ for all large $n$ gives a confidence sequence at any level $\alpha<1$.

(ii) If $\kappa>3$, then for every $\alpha\in(0,1)$ there is a valid threshold $b$ with $b(n)-\lambda_1\{2\log\log n+\kappa\log\log\log n\}$ bounded: the stitched threshold with $\delta=(\kappa-3)/2$.

So $2\lambda_1$ is the sharp first-order coefficient of the threshold and $3\lambda_1$ the sharp second-order one. Part (ii) is classical stitching. Part (i) is the Kolmogorov--Erdős--Feller test carried through the reduction of Theorem 5. What is specific here is the proof that the negative-eigenvalue walk $\widetilde W$ cannot spoil the lower-class times of $W$.

*Proof.* (ii) is the stitched construction above. For (i) it suffices to take $\kappa=3$, because $b_{\kappa,C}\le b_{3,C}$ eventually when $\kappa<3$. Fix $C$.

*Lower-class times for $W$.* Represent $W_n=B(n)$ for a standard Brownian motion $B$. Let $\psi(t)^2=2\log\log t+3\log\log\log t+C+2$, which is nondecreasing for large $t$. Then $\int^\infty\psi(t)e^{-\psi(t)^2/2}\,dt/t$ diverges, since the integrand is of exact order $(t\log t\log\log t)^{-1}$. By the Kolmogorov--Erdős--Feller integral test for Brownian motion (Erdős, 1942; Feller, 1943; Motoo, 1959), almost surely $B(t)\ge\sqrt t\,\psi(t)$ for arbitrarily large $t$.

*Discretisation.* $P(\sup_{s\in[0,1]}|B(n)-B(n-s)|>2\sqrt{\log n})\le4n^{-2}$. By Borel--Cantelli, almost surely for all large $n$ the maximum deviation of $B$ from $B(n)$ on $[n-1,n]$ is at most $2\sqrt{\log n}$. At a large $t$ with $B(t)\ge\sqrt t\psi(t)$, put $n=\lceil t\rceil$. Then
$$W_n/\sqrt n\ge\psi(n-1)\sqrt{(n-1)/n}-2\sqrt{\log n/n}=\psi(n)-o(1/\psi(n)),$$
so $W_n^2/n\ge\psi(n)^2-1$ for all large such $n$. Hence the random set $\mathcal N=\{n:W_n^2/n\ge2\log\log n+3\log\log\log n+C+1\}$ is infinite a.s., and it is a function of $W$ alone.

*Controlling $\widetilde W$.* Put $\varrho^2=\lambda_1/|\lambda_2|$. On $\{n\in\mathcal N,\ \widetilde W_n^2/n\le\varrho^2\}$,
$$Q_n\ge\lambda_1(2\log\log n+3\log\log\log n+C+1)-\lambda_1=b_{3,C}(n).$$
Condition on $W$. Because $\widetilde W$ is independent of $W$, $\mathcal N$ is then a fixed infinite set. Choose $m_1<m_2<\cdots$ in $\mathcal N$ with $m_{j+1}\ge m_j^3$, and let $A_j=\{|\widetilde W_{m_j}|\le\varrho\sqrt{m_j}\}$ and $\mathcal F_j=\sigma(\widetilde W_1,\dots,\widetilde W_{m_j})$.

- By the LIL, $|\widetilde W_{m_j}|\le\sqrt{3m_j\log\log m_j}$ for all large $j$, a.s.
- Given $\mathcal F_j$, $\widetilde W_{m_{j+1}}$ is $N(\widetilde W_{m_j},m_{j+1}-m_j)$, and $|\widetilde W_{m_j}|/\sqrt{m_{j+1}}\le\sqrt{3\log\log m_j}/m_j\to0$.
- Hence $P(A_{j+1}\mid\mathcal F_j)\to P(|N(0,1)|\le\varrho)>0$ a.s., so $\sum_jP(A_{j+1}\mid\mathcal F_j)=\infty$ a.s.

By Lévy's conditional Borel--Cantelli lemma, $A_j$ occurs infinitely often a.s. Integrating over $W$ shows that, almost surely, $Q_n\ge b_{3,C}(n)$ infinitely often. $\square$

**Proposition 5A (a model-free checkpoint repair).** Consider any model in which each $E_n(\theta)$ is a fixed-$n$ e-value. Fix checkpoints $n_1<n_2<\cdots$ and weights $w_k\ge0$ with $\sum_kw_k\le1$. For $n_k\le n<n_{k+1}$ let $C^{\mathrm{cp}}_n=\{\vartheta:E_{n_k}(\vartheta)<1/(\alpha w_k)\}$, and let $C^{\mathrm{cp}}_n=\Theta$ for $n<n_1$. Then $P_\theta(\exists n:\theta\notin C^{\mathrm{cp}}_n)\le\alpha$.

*Proof.* The event $\{\exists n:\theta\notin C^{\mathrm{cp}}_n\}$ equals $\bigcup_k\{E_{n_k}(\theta)\ge1/(\alpha w_k)\}$. By Markov's inequality and the union bound its probability is at most $\sum_k\alpha w_k\le\alpha$. $\square$

With $n_k=2^k$ and the stitching weights above, the threshold is $\log(1/\alpha)+\log\log n+(1+\delta)\log\log\log n+O(1)$. Its leading coefficient is $1$, against the sharp $2\lambda_1=0.618$ of Theorem 6, which we prove only in the setting of Theorem 5. Proposition 5A is dominated by the time-mixture e-process of Ramdas and Wang (2025, Proposition 7.23) with weights on the checkpoints: $M_n=\sum_{k:n_k\le n}w_kE_{n_k}\ge w_kE_{n_k}$, so the sets $\{M_n<1/\alpha\}$ are contained in ours. Unlike Theorem 5, it needs no regularity, so it is the repair to use for irregular models.

*Remark.* Theorem 6 is about the iterated-logarithm *scale*. There is no single smallest valid boundary, since validity is a statement at level $\alpha$. For example, weights $w_k\propto1/(k\log k(\log\log k)^2)$ give valid thresholds $\lambda_1\{2L_2+3L_3+4\log L_3\}+O(1)$, which lie between every $\kappa>3$ and the impossible $\kappa=3$. At practical sample sizes the normal-mixture choice is tighter than any stitched choice we tried (Section 8). For horizons fixed in advance, Section 7.3 gives calibrated thresholds.

# 7. How bad is it at finite horizons?

Theorem 1 is a statement about monitoring forever. Pace et al.'s simulations, which looked nearly fine, ran over finite horizons. This section reconciles the two.

## 7.1 An exact Ornstein--Uhlenbeck representation

**Proposition 7.** Couple $S_n=B_1(n)$ and $T_n=B_2(n)$ with independent standard Brownian motions, and let $\mathbf X(t)=e^{-t/2}(B_1(e^t),B_2(e^t))$. Then $\mathbf X$ is a stationary two-dimensional Ornstein--Uhlenbeck process with covariance $e^{-|t-s|/2}I_2$, and $Q_n=\mathbf X(\log n)^\top A\,\mathbf X(\log n)$ exactly for every $n$. For fixed $R>1$, as $n_0\to\infty$,
$$
P\Big(\max_{n_0\le n\le Rn_0}Q_n\ge c\Big)\ \longrightarrow\ \pi(R,c):=P\Big(\sup_{0\le t\le\log R}\mathbf X(t)^\top A\mathbf X(t)\ge c\Big)
$$
at every continuity point $c$ of the limit law, that is, for all but countably many $c$.

*Proof.* The covariance follows from $\operatorname{Cov}(B(e^t),B(e^s))=e^{\min(t,s)}$, and stationarity from Brownian scaling. The identity is Theorem 1(i) with $S_n/\sqrt n=e^{-\log n/2}B_1(e^{\log n})$. By stationarity, $\max_{n_0\le n\le Rn_0}Q_n$ has the law of $M_{n_0}=\max\{F(\mathbf X(s)):s\in G_{n_0}\}$, where $F(x)=x^\top Ax$ and $G_{n_0}=\{\log(n/n_0):n_0\le n\le Rn_0\}\subset[0,\log R]$ has mesh at most $1/n_0$. By path continuity, $M_{n_0}\to\sup_{[0,\log R]}F(\mathbf X)$ almost surely, which gives convergence in distribution. $\square$

So, asymptotically, the uncoverage over $n_{\min}\le n\le n_{\max}$ depends only on the horizon ratio $R=n_{\max}/n_{\min}$. Increasing the burn-in at a fixed ratio does *not* make the problem disappear.

## 7.2 Failure exponents and the horizons at which error exceeds nominal

Let $D_c=\{x:x^\top Ax<c\}$, the unbounded region between the two branches of a hyperbola. Let $\mu(c)>0$ be the principal Dirichlet eigenvalue of $-\mathcal L$ on $D_c$, where $\mathcal Lf=\tfrac12\Delta f-\tfrac12x\cdot\nabla f$ is the OU generator. The OU operator has compact resolvent in $L^2$ of the standard Gaussian measure, and the Dirichlet restriction keeps this property. Hence $1-\pi(R,c)\sim\mathrm{const}\times R^{-\mu(c)}$ as $R\to\infty$: persistence decays like a *power* of the horizon ratio. We use this asymptotic only to interpret the results and to extrapolate beyond $R=10^6$; the tabulated values are computed directly.

*Numerical method* (script `ou_spectral.py`).

- *Discretisation:* second-order finite differences in the eigen-coordinates of $A$, with a Shortley--Weller treatment of the curved boundary and truncation at $|v|\le8$ along the negative-eigenvalue axis.
- *Eigenvalues:* Richardson-extrapolated from three meshes.
- *Survival curves:* Crank--Nicolson time stepping (step $0.01$, Rannacher start), Richardson-extrapolated from meshes $0.08$ and $0.04$.
- *Checks:* in the hardest case ($\varepsilon=0.01$, $R=10^6$), a finer-mesh extrapolation agrees to within 1% relative, and halving the time step changes the result by under 0.5%. The long-horizon decay rate matches the eigenvalue to within 3%, and the initial value matches the exact Gaussian mass of $D_c$.

{{OU_TABLE}}

Three features stand out.

1. The exponents are small, and they fall steeply with $\varepsilon$. So the excess error grows very slowly with the horizon, which is why simulations did not reveal the problem.
2. The single-time coverage $P(\mathbf X(0)\in D_c)$ is far above $1-\varepsilon$, because Markov's inequality is loose for this e-value. A large part of the error budget is unused at any one time, and it takes a long monitoring period to use it up.
3. Uncoverage first exceeds $\varepsilon$ at a horizon ratio of about {{R_eps20}} for $\varepsilon=0.20$, about $10^{{{L10_eps10}}}$ for $\varepsilon=0.10$, about $10^{{{L10_eps05}}}$ for $\varepsilon=0.05$ and about $10^{{{L10_eps01}}}$ for $\varepsilon=0.01$ (the last by power-law extrapolation).

At conventional levels, then, the split-naive sequence satisfies Pace et al.'s finite-horizon property up to substantial horizons: about $660$ times the first look at $90\%$ and about $8{,}600$ times at $95\%$. Beyond those horizons, and at looser levels, it does not.

## 7.3 Reproducing Pace et al.'s Gaussian tables, and calibrated thresholds

The Gaussian scheme can be simulated exactly (script `finite_horizon_mc.py`):

{{MC_TABLE}}

All 16 published Gaussian percentages agree with ours to within Monte Carlo error; the largest discrepancy is $2.65$ standard errors, and most are under one. The continuum values $\pi(R,c)$ are slightly higher than the discrete-time ones, because discrete monitoring after a finite burn-in misses a few crossings.

If the horizon ratio $R$ is fixed in advance, the threshold can be calibrated: replace $\log(1/\varepsilon)$ by the $c^*(\varepsilon,R)$ that solves $\pi(R,c^*)=\varepsilon$. These numerically calibrated thresholds target the continuum limit, and in our discrete-time simulations the continuum crossing probability was always the larger; their numerical accuracy is stated in the caption:

{{CAL_TABLE}}

## 7.4 A finite number of looks

Group-sequential designs look at the data only a few times. The next table takes $L$ looks spaced geometrically between $n_{\mathrm{first}}$ and $Rn_{\mathrm{first}}$ (exact simulation with $10^6$ replications, script `looks_and_mixture_controls.py`):

{{LOOKS_TABLE}}

As in classical repeated significance testing (Armitage, McPherson and Rowe, 1969; Pocock, 1977), the error grows with the number of looks. Because Markov's inequality is loose here, however, for the two schedules and levels shown, even 50 looks at the nominal threshold stay within $\varepsilon$ for $\varepsilon\le0.10$.

## 7.5 The Poisson example

Pace et al.'s Table 5 reports split-naive uncoverages of $43.00$, $20.34$, $8.76$ and $1.34\%$ for Poisson($1$) data. Their text gives $n_{\min}=100$ and $n_{\max}=10{,}000$, while the caption says total sizes 100 to 10,000. We simulated their equation (13) exactly as printed, including the $0.5/n$ replacement, under both readings with $100{,}000$ replications, and obtained {{POIS}}. These values are close to the Gaussian values at the same horizon ratio, as Theorem 2 and Proposition 7 predict, and less than half the published figures.

To locate the difference we also simulated the same paper's *mixture* intervals. Our Gaussian normal-mixture uncoverages reproduce their Table 2 ({{GMIX}}). Our Poisson Gamma-mixture uncoverages are {{PMIX}}, against $9.20$, $5.44$, $3.34$ and $1.00\%$ published. The Poisson discrepancy therefore also affects a method that does not involve splitting; its source remains unresolved. Our conclusions do not depend on it: Proposition 0 and Theorem 2 cover the Poisson model whatever the finite-horizon numbers.

# 8. Efficiency: repaired intervals against existing anytime-valid intervals

We compared seven intervals for a Gaussian mean ($\sigma=1$, $\theta=0$, two halves of $n$ each, $2n$ in total) over $1\le n\le10^5$, with $4{,}000$ replications (script `repair_efficiency.py`):

* split-naive, which is not time-uniform;
* the *repaired split* interval (Theorem 5), with the normal-mixture $g$ ($\rho=1$, matching Pace et al.'s prior) and with the stitched $g$ ($\eta=2$, $\delta=1/2$);
* the *checkpoint* repair (Proposition 5A, $n_k=2^k$);
* the *running-MLE* ("split-exact") confidence sequence of Wasserman et al. (2020) on the pooled stream, with the mean of the previous observations as plug-in;
* the *Robbins normal-mixture* and *stitched* confidence sequences on the pooled $2n$ observations.

{{EFF_TABLE}}

Every valid method kept its time-uniform non-coverage below $\alpha$, and the split-naive interval exceeded it at $\alpha=0.20$. At $\alpha=0.05$ the split-naive interval stayed below $\alpha$ over this horizon (a ratio of $10^5$, $4.2\%$), consistent with Section 7.2, and it was the narrowest of the displayed intervals at $n=10^4$ and $10^5$. This is the practical sense in which the conjecture "almost" holds at conventional levels.

At the tabulated sizes $n\ge10^3$ the repaired split interval with the mixture boundary is about 13--14% wider than the pooled normal-mixture interval and about 10--11% wider than the running-MLE interval. For the stitched repair against the pooled stitched interval, the ratio of half-widths tends in probability to $\sqrt{4\lambda_1}=\sqrt{\sqrt5-1}\approx1.11$. Its almost-sure $\limsup$ is $\sqrt{4\lambda_1+4}$, because of the $(\bar X_n-\bar Y_n)^2$ term. This is the price of validating on half of the data. At practical sizes the stitched boundaries, despite their better growth rate, are wider than the mixture ones, and the model-free checkpoint repair is wider again.

In practice:

- for Gaussian data, or data satisfying Theorem 5's joint sub-Gaussian condition with the Gaussian working likelihood, if one wants to keep the split statistic, use the repaired threshold of Theorem 5;
- for irregular models with genuine fixed-$n$ split e-values, use the checkpoint or time-mixture repair (Proposition 5A);
- if efficiency matters, use a pooled anytime-valid method.

# 9. Discussion

*What changes.* A fixed-$n$ e-value recomputed as data accumulate need not be an e-process (Ramdas and Wang, 2025, Section 7.9). For the split likelihood ratio the consequence is now exact: zero asymptotic persistence (Proposition 0), with explicit rate constants in regular models. The conjecture of Pace et al. (2023, Section 3) is false. Their simulations correctly showed small finite-horizon excesses at nominal $80\%$. Proposition 7 and Section 7.2 explain why the excesses were small, and why the error probability tends to one as the horizon grows.

*What does not change.* Fixed-$n$ universal inference is valid. The running-MLE construction gives exact confidence sequences (Wasserman et al., 2020, Theorem 11). Pace et al.'s mixture sequences, and their separate Section 6 conjecture, are untouched.

*Limitations.*

- Theorem 6 is specific to the Gaussian scheme. For other regular models only the rate (Theorem 2) and the model-free repair (Proposition 5A) are proved.
- The continuum calculations in Section 7 are numerical, not closed form; their accuracy is stated with each table.
- The calibrated thresholds assume the Gaussian scheme and a large first look.
- Growing-$K$ and re-randomised cross-fitting, rate constants for irregular models and a plug-in repair for composite likelihoods remain open.
- Our prior-art search was extensive but bounded.

# 10. Reproducibility

All scripts are in the release archive (MIT licence) and use only NumPy, SciPy and SymPy. Checks raise errors rather than using `assert`, so they stay active under `python -O`. Each script contains negative controls that it must detect: a wrong constant, a non-dominating direction, the split-naive excess and a wrong Fisher information.

* `algebra_checks.py`: exact symbolic identities and eigenvalues (Theorems 1, 3 and 4, and the sandwich constant).
* `universality_checks.py`: the quadratic expansion of Theorem 2 for Poisson and Bernoulli streams, up to $n=10^8$.
* `finite_horizon_mc.py`: exact Monte Carlo of Pace et al.'s Gaussian and Poisson experiments.
* `ou_spectral.py`: eigenvalues, survival curves and calibrated thresholds for the OU continuum, with convergence checks.
* `looks_and_mixture_controls.py`: finite numbers of looks, and the mixture-row controls of Section 7.5.
* `repair_efficiency.py`: validity and widths of the repaired and comparison intervals.
* `fill_tables.py`: inserts every tabulated number into the manuscript from the result files.

**Assurance.** This is an unrefereed research release. The proofs were checked internally, by AI-assisted review and by the scripts above. That checking is not external peer review or independent reproduction.

# References

{{REFS}}
