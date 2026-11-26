"""Independent checks of signet.core: (1) DECODE's refutation distance as a
special case; (2) the test statistic, certificates and identified sets
recomputed with a different modelling layer (CVXPY with SCS).
Run from the repository root: python3 tests/verify_crosssolver.py"""
import itertools, os, sys
import numpy as np
import cvxpy as cp
_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.join(_ROOT, "experiments", "simulation"))
from decode_core import refutation_distance          # noqa: E402
from signet.core import *                            # noqa: E402,F401,F403
import worked_example as W                           # noqa: E402
assert hasattr(W, "simulate")

rng = np.random.default_rng(11)
# ---------------------------------------------------------------- (1)
agree, total = 0, 0
for trial in range(60):
    n, K = rng.integers(5, 10), rng.integers(3, 7)
    S = (rng.random((n, K)) < 0.35).astype(float)
    S = S[:, S.sum(axis=0) > 0]
    K = S.shape[1]
    Z = np.ones((n, 1))
    d_dec = refutation_distance(Z, S, c=[1.0], tmax=K)[0]
    m = Model(1)
    for i in range(n):
        m.add_row(f"r{i}", [1.0], [f"a{k}" for k in range(K) if S[i, k]])
    shares, zeros = [], []
    for k in range(K):
        if len(m.rows_with(f"a{k}")) > 1:
            m.share(f"sh{k}", f"a{k}"); shares.append(f"sh{k}")
        m.zero(f"v{k}", f"a{k}", rung=3); zeros.append(f"v{k}")
    d_sig = np.inf
    for kk in range(0, K + 1):
        hit = False
        for R in itertools.combinations(zeros, kk):
            inc = shares + [z for z in zeros if z not in R]
            cap = capacity(m, inc, [1.0])
            if not (cap["above"] and cap["below"]):
                hit = True; break
        if hit:
            d_sig = kk; break
    total += 1; agree += (d_sig == d_dec)
print(f"(1) refutation distance = smallest tolerance that unidentifies the claim: {agree}/{total} ledgers agree")

# ---------------------------------------------------------------- (2)
rng2 = np.random.default_rng(W.SEED)
y, s, studies, info = W.simulate(rng2)
m = W.build(); m.set_data(y, s)
cpd = m.compile(W.RUNG3)
A, D, E = cpd["A"], cpd["D"], cpd["E"]
Winv = np.diag(1 / s ** 2)
c = m.cvec(W.FUNCS["Early-stage population"])
F = Fitter(m, W.RUNG3)
maxdiff = 0.0
for t0 in (-0.5, -0.35, -0.25, -0.1, 0.05):
    x = cp.Variable(m.p)
    cons = [D @ x <= 0, E @ x == 0, c @ x == t0]
    pr = cp.Problem(cp.Minimize(cp.quad_form(y - A @ x, Winv)), cons)
    pr.solve(solver="SCS", eps=1e-9, max_iters=200000)
    T_ref = pr.value
    T_sig = F.fit(c, t0)["T"]
    maxdiff = max(maxdiff, abs(T_ref - T_sig))
print(f"(2a) T(theta0) by two solvers: largest absolute difference {maxdiff:.2e}")
# certificate: lexicographic rule re-solved with cvxpy
n_, mD, mE = A.shape[0], D.shape[0], E.shape[0]
lam, nu, kap = cp.Variable(n_), cp.Variable(mD), cp.Variable(mE)
base = [A.T @ lam + D.T @ nu + E.T @ kap == c, nu >= 0]
p1 = cp.Problem(cp.Minimize(cp.sum(nu)), base); p1.solve(solver="CLARABEL")
p2 = cp.Problem(cp.Minimize(cp.quad_form(lam, np.diag(s ** 2))), base + [cp.sum(nu) <= p1.value + 1e-9])
p2.solve(solver="CLARABEL")
ct = certificate(m, W.RUNG3, W.FUNCS["Early-stage population"], "hi")
print(f"(2b) certificate sd {ct['sd']:.6f} vs {np.sqrt(p2.value):.6f}; estimate {ct['lam'] @ y:+.6f} vs {lam.value @ y:+.6f}")
# identified sets at the true means
mu_true = np.array([d["mu_true"] for d in info])
worst = 0.0
for fn, cf in W.FUNCS.items():
    cc = m.cvec(cf)
    for inc in (W.RUNG2, W.RUNG3):
        cq = m.compile(inc)
        x = cp.Variable(m.p)
        cons = [cq["A"] @ x == mu_true, cq["D"] @ x <= 0, cq["E"] @ x == 0]
        hi = cp.Problem(cp.Maximize(cc @ x), cons); hi.solve(solver="CLARABEL")
        lo = cp.Problem(cp.Minimize(cc @ x), cons); lo.solve(solver="CLARABEL")
        idf = identified(m, inc, cf, mu_true)
        for a, b in ((hi.value, idf["hi"]), (lo.value, idf["lo"])):
            if np.isfinite(a) and np.isfinite(b):
                worst = max(worst, abs(a - b))
print(f"(2c) identified sets at the true means, two solvers: largest difference {worst:.2e}")
# exact and fast ranks along the early-population profile
Fe, Ff = Fitter(m, W.RUNG3, exact_rank=True), Fitter(m, W.RUNG3, exact_rank=False)
same = 0; tot = 0
for t0 in np.linspace(-0.6, 0.1, 15):
    fe, ff = Fe.fit(c, t0), Ff.fit(c, t0)
    if fe["ok"] and fe["T"] > 1e-9:
        tot += 1; same += Fe.rank(fe) == Ff.rank(ff)
print(f"(2d) exact and fast rank agree at {same}/{tot} points")
