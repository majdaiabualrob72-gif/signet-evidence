"""
decode_core.py
Reference implementation of the linear-Gaussian core of DECODE
(Diagnostic Error-Correcting Decoding of Evidence).

Model
-----
    y = Z theta + S b + u + e
    e_i ~ N(0, s_i^2)                  sampling error (reported standard errors)
    u_i ~ N(0, tau^2)                  idiosyncratic between-study deviation
    b_k = gamma_k * beta_k             violation of assumption k
    gamma_k ~ Bernoulli(pi_k)          prior probability that assumption k fails
    beta_k ~ N(0, sigma_k^2)           size of the violation if it fails
    theta ~ flat                       target estimands

    Z : n x m   which estimand each evidence item measures
    S : n x K   assumption signature matrix (exposure of item i to assumption k)

Everything below is exact for this model: the posterior over the 2^K
assumption configurations is enumerated, with closed-form Gaussian
marginal likelihoods. No sampling is involved, so results are reproducible
to floating-point precision.
"""
import itertools
import numpy as np
from scipy import linalg as sla
from scipy import optimize as sopt
from scipy import stats as sst
from scipy.special import logsumexp

DEFAULT_TAU_GRID = np.array([0.0, 0.02, 0.04, 0.06, 0.08, 0.12, 0.16, 0.24])


# ----------------------------------------------------------------------
# Ledger
# ----------------------------------------------------------------------
class Ledger:
    """The evidence ledger: items (rows), estimands, assumptions, priors."""

    def __init__(self, y, s, Z, S, pi, sigma, item_labels=None,
                 estimand_names=None, assumption_names=None, het=None):
        self.y = np.asarray(y, float)
        self.s = np.asarray(s, float)
        self.Z = np.atleast_2d(np.asarray(Z, float))
        if self.Z.shape[0] != len(self.y):
            self.Z = self.Z.T
        self.S = np.asarray(S, float)
        self.pi = np.asarray(pi, float)
        self.sigma = np.asarray(sigma, float)
        self.n, self.m = self.Z.shape
        self.K = self.S.shape[1]
        self.item_labels = item_labels or [f"item{i}" for i in range(self.n)]
        self.estimand_names = estimand_names or [f"theta{j}" for j in range(self.m)]
        self.assumption_names = assumption_names or [f"A{k+1}" for k in range(self.K)]
        # het[i] = 1 if the idiosyncratic heterogeneity term applies to item i
        self.het = np.ones(self.n) if het is None else np.asarray(het, float)

    def append(self, y, s, Z, S, labels=None, het=None):
        """Return a new ledger with extra evidence items appended."""
        y = np.atleast_1d(y); s = np.atleast_1d(s)
        Z = np.atleast_2d(Z); S = np.atleast_2d(S)
        het_new = np.ones(len(y)) if het is None else np.asarray(het, float)
        return Ledger(np.r_[self.y, y], np.r_[self.s, s],
                      np.vstack([self.Z, Z]), np.vstack([self.S, S]),
                      self.pi, self.sigma,
                      self.item_labels + (labels or [f"new{i}" for i in range(len(y))]),
                      self.estimand_names, self.assumption_names,
                      np.r_[self.het, het_new])


# ----------------------------------------------------------------------
# Mixture of Gaussians (posterior of any linear functional)
# ----------------------------------------------------------------------
class Mixture:
    def __init__(self, w, mu, sd):
        self.w = np.asarray(w, float)
        self.mu = np.asarray(mu, float)
        self.sd = np.maximum(np.asarray(sd, float), 1e-10)

    def mean(self):
        return float(np.sum(self.w * self.mu))

    def var(self):
        return float(np.sum(self.w * (self.sd ** 2 + self.mu ** 2)) - self.mean() ** 2)

    def sd_(self):
        return float(np.sqrt(max(self.var(), 0.0)))

    def cdf(self, x):
        return float(np.sum(self.w * sst.norm.cdf((x - self.mu) / self.sd)))

    def quantile(self, p):
        lo = float(np.min(self.mu - 12 * self.sd)) - 1e-6
        hi = float(np.max(self.mu + 12 * self.sd)) + 1e-6
        return float(sopt.brentq(lambda x: self.cdf(x) - p, lo, hi, xtol=1e-10))

    def interval(self, level=0.95):
        a = (1 - level) / 2
        return self.quantile(a), self.quantile(1 - a)

    def logpdf(self, x, extra_var=0.0):
        sd = np.sqrt(self.sd ** 2 + extra_var)
        return float(logsumexp(np.log(self.w + 1e-300) + sst.norm.logpdf(x, self.mu, sd)))


# ----------------------------------------------------------------------
# Exact decoder
# ----------------------------------------------------------------------
class Posterior:
    def __init__(self, ledger, comps, logw, tau_grid):
        self.L = ledger
        self.comps = comps            # list of (tau_index, gamma tuple, mu, cov)
        self.logw = logw
        self.w = np.exp(logw)
        self.tau_grid = tau_grid

    # --- assumption-level summaries -----------------------------------
    def pip(self):
        K = self.L.K
        out = np.zeros(K)
        for w, (_, g, _, _) in zip(self.w, self.comps):
            for k in g:
                out[k] += w
        return out

    def config_posterior(self):
        """Posterior over assumption configurations, tau marginalised."""
        d = {}
        for w, (_, g, _, _) in zip(self.w, self.comps):
            d[g] = d.get(g, 0.0) + w
        return sorted(d.items(), key=lambda kv: -kv[1])

    def tau_posterior(self):
        out = np.zeros(len(self.tau_grid))
        for w, (ti, _, _, _) in zip(self.w, self.comps):
            out[ti] += w
        return out

    def joint_failure(self, ks):
        """Posterior probability that at least one assumption in ks fails."""
        ks = set(ks)
        return float(sum(w for w, (_, g, _, _) in zip(self.w, self.comps) if ks & set(g)))

    # --- linear functionals f = c' theta + a' b ------------------------
    def functional(self, c=None, a=None):
        m, K = self.L.m, self.L.K
        c = np.zeros(m) if c is None else np.asarray(c, float)
        a = np.zeros(K) if a is None else np.asarray(a, float)
        mus, sds = [], []
        for (_, g, mu, cov) in self.comps:
            v = np.r_[c, a[list(g)]]
            mus.append(v @ mu)
            sds.append(np.sqrt(max(v @ cov @ v, 0.0)))
        return Mixture(self.w, mus, sds)

    def theta(self, j=0):
        c = np.zeros(self.L.m); c[j] = 1.0
        return self.functional(c=c)

    def b(self, k):
        a = np.zeros(self.L.K); a[k] = 1.0
        return self.functional(a=a)

    # --- prediction of a new evidence item ------------------------------
    def predict(self, z_new, s_new_row, var_new):
        """Posterior predictive mixture for a new item's reported estimate."""
        mix = self.functional(c=z_new, a=s_new_row)
        return Mixture(mix.w, mix.mu, np.sqrt(mix.sd ** 2 + var_new))


def decode(L, tau_grid=None, tau_weights=None, max_size=None):
    """Exact posterior over assumption configurations and targets."""
    tau_grid = DEFAULT_TAU_GRID if tau_grid is None else np.asarray(tau_grid, float)
    if tau_weights is None:
        # half-normal(0.10) prior on tau, discretised on the grid
        tau_weights = sst.norm.pdf(tau_grid, 0, 0.10)
        if len(tau_grid) > 1:
            # trapezoid-type cell widths so unequal spacing is respected
            widths = np.gradient(tau_grid)
            tau_weights = tau_weights * widths
    tau_weights = np.asarray(tau_weights, float) / np.sum(tau_weights)
    n, m, K = L.n, L.m, L.K
    max_size = K if max_size is None else max_size
    X = np.hstack([L.Z, L.S])
    configs = [g for r in range(max_size + 1) for g in itertools.combinations(range(K), r)]
    logpi, log1mpi = np.log(L.pi), np.log1p(-L.pi)
    base_prior = np.sum(log1mpi)
    comps, logw = [], []
    for ti, (tau, wt) in enumerate(zip(tau_grid, tau_weights)):
        if wt <= 0:
            continue
        w = 1.0 / (L.s ** 2 + (tau * L.het) ** 2)
        G = X.T @ (w[:, None] * X)
        h = X.T @ (w * L.y)
        yWy = L.y @ (w * L.y)
        ld = -np.sum(np.log(w))
        for g in configs:
            idx = list(range(m)) + [m + k for k in g]
            Lam = G[np.ix_(idx, idx)].copy()
            for j, k in enumerate(g):
                Lam[m + j, m + j] += 1.0 / L.sigma[k] ** 2
            cf = sla.cho_factor(Lam, lower=True)
            mu = sla.cho_solve(cf, h[idx])
            logdet = 2.0 * np.sum(np.log(np.diag(cf[0])))
            gl = list(g)
            logml = -0.5 * (ld + np.sum(np.log(L.sigma[gl] ** 2)) + logdet + yWy - h[idx] @ mu)
            lp = base_prior + np.sum(logpi[gl] - log1mpi[gl]) + np.log(wt)
            cov = sla.cho_solve(cf, np.eye(len(idx)))
            comps.append((ti, g, mu, cov))
            logw.append(logml + lp)
    logw = np.asarray(logw)
    logw -= logsumexp(logw)
    return Posterior(L, comps, logw, tau_grid)


# ----------------------------------------------------------------------
# Structural (data-free) analysis of a ledger
# ----------------------------------------------------------------------
def parity_basis(Z, tol=1e-10):
    """Rows span the left null space of Z: P Z = 0 (the parity checks)."""
    U, sv, _ = np.linalg.svd(Z, full_matrices=True)
    r = int(np.sum(sv > tol * max(Z.shape) * (sv[0] if len(sv) else 1.0)))
    return U[:, r:].T


def detectability(Z, S):
    """Norm of each assumption's parity image; 0 means undetectable."""
    P = parity_basis(Z)
    return np.linalg.norm(P @ S, axis=0)


def alias_classes(Z, S, tol=1e-8):
    """Group assumptions whose single violations leave collinear syndromes."""
    P = parity_basis(Z)
    PS = P @ S
    norms = np.linalg.norm(PS, axis=0)
    K = S.shape[1]
    classes, assigned = [], [False] * K
    undetectable = [k for k in range(K) if norms[k] < tol]
    for k in undetectable:
        assigned[k] = True
    for k in range(K):
        if assigned[k]:
            continue
        cls = [k]; assigned[k] = True
        for l in range(k + 1, K):
            if assigned[l]:
                continue
            cosang = abs(PS[:, k] @ PS[:, l]) / (norms[k] * norms[l])
            if cosang > 1 - tol:
                cls.append(l); assigned[l] = True
        classes.append(cls)
    return classes, undetectable


def _null(A, tol=1e-9):
    """Null space with an absolute tolerance (inputs are O(1) by construction)."""
    A = np.atleast_2d(A)
    if A.shape[0] == 0:
        return np.eye(A.shape[1])
    _, sv, Vt = np.linalg.svd(A, full_matrices=True)
    rank = int(np.sum(sv > tol))
    return Vt[rank:].T


def blind_dimension(Z, S, tol=1e-9):
    P = parity_basis(Z)
    PS = P @ S
    return S.shape[1] - np.linalg.matrix_rank(PS, tol=tol)


def refutation_distance(Z, S, c=None, a=None, tmax=7, tol=1e-8, return_all=False):
    """Smallest number of simultaneous assumption violations that shifts
    the functional f = c' theta + a' b while leaving no syndrome.
    Returns (distance, supporting set) or (inf, None)."""
    n, m = Z.shape
    K = S.shape[1]
    c = np.zeros(m) if c is None else np.asarray(c, float)
    a = np.zeros(K) if a is None else np.asarray(a, float)
    P = parity_basis(Z)
    Zp = np.linalg.pinv(Z)
    found = []
    for t in range(1, min(tmax, K) + 1):
        for T in itertools.combinations(range(K), t):
            T = list(T)
            PS_T = P @ S[:, T]
            N = _null(PS_T)
            if N.shape[1] == 0:
                continue
            shift = (c @ Zp @ S[:, T] - a[T]) @ N
            if np.any(np.abs(shift) > tol):
                if not return_all:
                    return t, T
                found.append(T)
        if found:
            return t, found
    return np.inf, None


# ----------------------------------------------------------------------
# Bias budgets (a joint, multi-assumption breakdown analysis)
# ----------------------------------------------------------------------
def bias_budget(L, c=None, a=None, direction=None, thresh=0.0, tau=0.0,
                level=0.95):
    """Minimum total assumption violation, in prior standard deviations,

        minimise   sum_k |b_k| / sigma_k

    subject to
      (coherence)  the bias-corrected data y - S b show no residual
                   contradiction: Q(b) <= chi2_{level, n-m}, where Q is the
                   weighted residual sum of squares after fitting theta;
      (tenability) optional. With direction = +1 the statement
                   f = c' theta + a' b >= thresh must be tenable, meaning the
                   upper confidence limit of f under bias pattern b reaches
                   thresh; direction = -1 is the mirror image.

    With direction=None the result is the coherence budget B0: the least
    assumption failure that must be conceded for the literature to be
    mutually consistent. Convex (linear objective, one linear and one
    convex quadratic constraint); solved by SLSQP from several starts.
    """
    n, m, K = L.n, L.m, L.K
    c = np.zeros(m) if c is None else np.asarray(c, float)
    a = np.zeros(K) if a is None else np.asarray(a, float)
    W = 1.0 / (L.s ** 2 + (tau * L.het) ** 2)
    ZtWZ = L.Z.T @ (W[:, None] * L.Z)
    V = np.linalg.inv(ZtWZ)
    M = V @ (L.Z.T * W)                       # theta_hat = M (y - S b)
    Pres = np.eye(n) - L.Z @ M
    r0 = Pres @ L.y
    Rm = Pres @ L.S
    f0 = c @ M @ L.y
    gvec = a - c @ M @ L.S                    # f_hat(b) = f0 + gvec @ b
    se_f = float(np.sqrt(max(c @ V @ c, 0.0)))
    zq = sst.norm.ppf(0.5 + level / 2)
    chi = sst.chi2.ppf(level, n - m)
    sw = np.sqrt(W)
    A_q = sw[:, None] * Rm; r_q = sw * r0     # Q(b) = || r_q - A_q b ||^2

    def unpack(x):
        return x[:K] - x[K:]

    cost = np.r_[1 / L.sigma, 1 / L.sigma]
    cons = [dict(type="ineq",
                 fun=lambda x: chi - np.sum((r_q - A_q @ unpack(x)) ** 2),
                 jac=lambda x: np.r_[2 * A_q.T @ (r_q - A_q @ unpack(x)),
                                     -2 * A_q.T @ (r_q - A_q @ unpack(x))])]
    if direction is not None:
        cons.append(dict(type="ineq",
                         fun=lambda x: direction * (f0 + gvec @ unpack(x) - thresh) + zq * se_f,
                         jac=lambda x: direction * np.r_[gvec, -gvec]))
    q_min, _ = min_incoherence(L, tau)
    if q_min > chi:
        return dict(budget=np.inf, b=None,
                    status="incoherent at this tau: no listed bias restores coherence "
                           "(minimum Q %.2f exceeds %.2f); raise tau or extend the ledger" % (q_min, chi))
    best = None
    rng = np.random.default_rng(0)
    # starts: least-squares bias fit, zero, and random points
    b_ls = np.linalg.lstsq(A_q, r_q, rcond=None)[0]
    starts = [np.r_[np.maximum(b_ls, 0), np.maximum(-b_ls, 0)], np.zeros(2 * K)]
    starts += [np.abs(rng.normal(0, 0.3, 2 * K)) for _ in range(6)]
    for x0 in starts:
        res = sopt.minimize(lambda x: cost @ x, x0, jac=lambda x: cost,
                            bounds=[(0, None)] * (2 * K), constraints=cons,
                            method="SLSQP", options=dict(maxiter=500, ftol=1e-12))
        feas = all(cn["fun"](res.x) >= -1e-6 for cn in cons)
        if feas and (best is None or res.fun < best.fun - 1e-10):
            best = res
    if best is None:
        return dict(budget=np.inf, b=None, status="infeasible")
    bvec = unpack(best.x)
    bvec[np.abs(bvec) < 1e-6] = 0.0
    return dict(budget=float(np.sum(np.abs(bvec) / L.sigma)), b=bvec,
                theta=M @ (L.y - L.S @ bvec), f=float(f0 + gvec @ bvec),
                se_f=se_f, status="optimal",
                n_engaged=int(np.sum(bvec != 0)))


def min_incoherence(L, tau=0.0):
    """Smallest weighted residual sum of squares attainable by any bias
    pattern (no penalty), and the coherence threshold it must not exceed."""
    n, m = L.n, L.m
    W = 1.0 / (L.s ** 2 + (tau * L.het) ** 2)
    X = np.hstack([L.Z, L.S])
    sw = np.sqrt(W)
    coef = np.linalg.lstsq(sw[:, None] * X, sw * L.y, rcond=None)[0]
    res = sw * (L.y - X @ coef)
    return float(res @ res), float(sst.chi2.ppf(0.95, n - m))


def coherent_tau(L, tau0=0.0, tau_max=2.0):
    """Smallest heterogeneity SD, not below tau0, at which some bias pattern
    makes the literature coherent. Equal to tau0 whenever tau0 suffices."""
    q, chi = min_incoherence(L, tau0)
    if q <= chi:
        return tau0
    lo, hi = tau0, max(2 * tau0, 0.01)
    while min_incoherence(L, hi)[0] > min_incoherence(L, hi)[1] and hi < tau_max:
        lo, hi = hi, 2 * hi
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        q, chi = min_incoherence(L, mid)
        lo, hi = (lo, mid) if q <= chi else (mid, hi)
    return hi


def refutation_premium(L, c=None, a=None, claim_sign=-1.0, thresh=0.0, tau=None):
    """Claim: sign * (f - thresh) > 0, for example claim_sign=-1 states f < thresh.
    Returns the coherence budget B0, the least budget under which the claim is
    tenable, the least budget under which its negation is tenable, and the
    refutation premium (negation budget minus B0).

    tau: heterogeneity SD used in the weights. If None, the posterior mean from
    the decoder is used. In either case tau is raised to the smallest value at
    which coherence is attainable, and the value used is returned."""
    if tau is None:
        post = decode(L)
        tau = float(post.tau_posterior() @ post.tau_grid)
    tau_used = coherent_tau(L, tau)
    B0 = bias_budget(L, tau=tau_used)
    Bc = bias_budget(L, c, a, direction=claim_sign, thresh=thresh, tau=tau_used)
    Bn = bias_budget(L, c, a, direction=-claim_sign, thresh=thresh, tau=tau_used)
    return dict(B0=B0, claim=Bc, negation=Bn,
                premium=Bn["budget"] - B0["budget"],
                support_premium=Bc["budget"] - B0["budget"],
                tau_used=tau_used)


# ----------------------------------------------------------------------
# Value of a candidate new study (expected information gain, in bits)
# ----------------------------------------------------------------------
def expected_information_gain(post, z_new, s_new_row, var_new, feature,
                              prune=1e-9, grid=3001):
    """Mutual information (bits) between a discrete feature of the
    assumption configuration and the result of a candidate study.
    feature: function gamma_tuple -> hashable label."""
    m = post.L.m
    keep = post.w > prune
    w = post.w[keep] / post.w[keep].sum()
    comps = [cp for cp, kp in zip(post.comps, keep) if kp]
    labels, mus, vs = [], [], []
    z_new = np.asarray(z_new, float); s_new_row = np.asarray(s_new_row, float)
    for (_, g, mu, cov) in comps:
        v = np.r_[z_new, s_new_row[list(g)]]
        mus.append(v @ mu); vs.append(v @ cov @ v + var_new)
        labels.append(feature(g))
    mus, vs = np.array(mus), np.array(vs)
    uniq = sorted(set(labels), key=str)
    lab_idx = np.array([uniq.index(l) for l in labels])
    lo = np.min(mus - 8 * np.sqrt(vs)); hi = np.max(mus + 8 * np.sqrt(vs))
    yg = np.linspace(lo, hi, grid)
    dens = w[:, None] * sst.norm.pdf(yg[None, :], mus[:, None], np.sqrt(vs)[:, None])
    joint = np.vstack([dens[lab_idx == j].sum(axis=0) for j in range(len(uniq))])
    py = joint.sum(axis=0)
    pl = np.array([w[lab_idx == j].sum() for j in range(len(uniq))])
    H_prior = -np.sum(pl[pl > 0] * np.log2(pl[pl > 0]))
    with np.errstate(divide="ignore", invalid="ignore"):
        post_l = np.where(py > 0, joint / py, 0.0)
        ent = -np.sum(np.where(post_l > 0, post_l * np.log2(post_l), 0.0), axis=0)
    H_post = np.trapezoid(ent * py, yg)
    return float(H_prior - H_post), float(H_prior)


# ----------------------------------------------------------------------
# Conventional comparators
# ----------------------------------------------------------------------
def re_meta(y, s, hksj=True):
    """DerSimonian-Laird random-effects meta-analysis, optional
    Hartung-Knapp-Sidik-Jonkman interval (modified, never narrower)."""
    y = np.asarray(y, float); s = np.asarray(s, float); k = len(y)
    w = 1 / s ** 2
    mu_fe = np.sum(w * y) / np.sum(w)
    Q = np.sum(w * (y - mu_fe) ** 2)
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - (k - 1)) / C) if k > 1 else 0.0
    ws = 1 / (s ** 2 + tau2)
    mu = np.sum(ws * y) / np.sum(ws)
    se = np.sqrt(1 / np.sum(ws))
    if hksj and k > 1:
        q = np.sum(ws * (y - mu) ** 2) / (k - 1)
        se = se * np.sqrt(max(1.0, q))
        crit = sst.t.ppf(0.975, k - 1)
        p = 2 * (1 - sst.t.cdf(abs(mu / se), k - 1))
    else:
        crit = sst.norm.ppf(0.975)
        p = 2 * (1 - sst.norm.cdf(abs(mu / se)))
    I2 = max(0.0, (Q - (k - 1)) / Q) if (k > 1 and Q > 0) else 0.0
    return dict(est=mu, se=se, lo=mu - crit * se, hi=mu + crit * se, p=p,
                tau2=tau2, I2=I2, k=k, Q=Q)


def dense_bias_model(L, tau_grid=None):
    """Closest established comparator: every listed bias is always active
    with a Gaussian prior matched to the prior second moment
    (variance pi_k * sigma_k^2). No sparsity, no diagnosis."""
    L2 = Ledger(L.y, L.s, L.Z, L.S, np.full(L.K, 1 - 1e-12),
                np.sqrt(L.pi) * L.sigma, L.item_labels, L.estimand_names,
                L.assumption_names, L.het)
    tau_grid = DEFAULT_TAU_GRID if tau_grid is None else tau_grid
    # only the all-on configuration carries weight
    n, m, K = L2.n, L2.m, L2.K
    X = np.hstack([L2.Z, L2.S])
    tw = sst.norm.pdf(tau_grid, 0, 0.10) * (np.gradient(tau_grid) if len(tau_grid) > 1 else 1.0)
    tw = tw / tw.sum()
    comps, logw = [], []
    g = tuple(range(K))
    for ti, (tau, wt) in enumerate(zip(tau_grid, tw)):
        w = 1.0 / (L2.s ** 2 + (tau * L2.het) ** 2)
        Lam = X.T @ (w[:, None] * X)
        for k in range(K):
            Lam[m + k, m + k] += 1.0 / L2.sigma[k] ** 2
        h = X.T @ (w * L2.y)
        cf = sla.cho_factor(Lam, lower=True)
        mu = sla.cho_solve(cf, h)
        logdet = 2.0 * np.sum(np.log(np.diag(cf[0])))
        logml = -0.5 * (-np.sum(np.log(w)) + logdet + L2.y @ (w * L2.y) - h @ mu)
        comps.append((ti, g, mu, sla.cho_solve(cf, np.eye(m + K))))
        logw.append(logml + np.log(wt))
    logw = np.asarray(logw); logw -= logsumexp(logw)
    return Posterior(L2, comps, logw, tau_grid)


def auroc(score, truth):
    score = np.asarray(score, float); truth = np.asarray(truth, bool)
    pos, neg = score[truth], score[~truth]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    ranks = sst.rankdata(np.r_[pos, neg])
    return float((ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))
