"""GENESIS VI independent computational audit of signet.core.
Each check uses a different formulation or solver from the package, or a Monte Carlo experiment.
(A) capacity: package primal LP (HiGHS) vs Farkas dual feasibility (CVXPY + Clarabel), random assumption subsets
(B) T(theta0): package QP (Clarabel, null-space form) vs CVXPY + OSQP on the original constrained problem
(C) size of the conditional chi-squared test under exactly Gaussian cells with known variances
(D) confidence-set inversion: coverage of theta under the true assumption set, and cases where theta is
    accepted by the test but lies outside the reported set (missed acceptance regions); dense-grid scan
(E) certificate validity: residual of c = A'lam + D'nu + E'kap, sign of nu
Run from the repository root: python3 audit/independent_checks.py"""
import os, sys, time, itertools
import numpy as np, cvxpy as cp, scipy.stats as sst
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "experiments", "simulation"))
from signet.core import *            # noqa
import worked_example as W           # noqa
rng = np.random.default_rng(6060)
m = W.build()
y0, s0, _, info = W.simulate(np.random.default_rng(W.SEED))
V0 = np.diag(s0 ** 2); m.set_data(y0, s0)
groups = W.RUNG3
TRUE = [g for g in W.RUNG3 if g != "Surrogate valid"]      # the worked example breaks surrogate validity
out = []
def say(s): print(s); out.append(s)
# ---------------- (A)
t = time.time(); agree = tot = 0; fails = []
for r in range(250):
    inc = [g for g in groups if rng.random() < 0.6]
    cpd = m.compile(inc); A, D, E = cpd["A"], cpd["D"], cpd["E"]
    for fn, cf in W.FUNCS.items():
        c = m.cvec(cf); cap = capacity(m, inc, cf)
        for side, sg in (("above", 1.0), ("below", -1.0)):
            lam, nu = cp.Variable(A.shape[0]), cp.Variable(D.shape[0])
            terms = A.T @ lam + D.T @ nu
            if E.shape[0]:
                kap = cp.Variable(E.shape[0]); terms = terms + E.T @ kap
            pr = cp.Problem(cp.Minimize(cp.norm1(terms - sg * c)), [nu >= 0])
            try:
                pr.solve(solver="CLARABEL")
            except cp.error.SolverError:
                pr.solve(solver="SCS", eps=1e-10, max_iters=200000)
            dual_ok = pr.value is not None and pr.value < 1e-6 * (1 + np.abs(c).sum())
            tot += 1; agree += (dual_ok == cap[side])
            if dual_ok != cap[side]: fails.append((inc, fn, side))
say(f"(A) capacity, primal LP vs Farkas dual: {agree}/{tot} agree ({time.time()-t:.0f} s)")
# ---------------- (B)
t = time.time(); worst = 0.0; worst_rel = 0.0; nb = 0
for r in range(40):
    y = np.array([d["mu_true"] for d in info]) + s0 * rng.standard_normal(len(s0))
    F = Fitter(m, TRUE, y, V0); cpd = m.compile(TRUE); A, D, E = cpd["A"], cpd["D"], cpd["E"]
    for fn in ("Early-stage population", "0-3 y", "Late minus early population"):
        c = m.cvec(W.FUNCS[fn])
        for t0 in (-0.6, -0.3, 0.0, 0.2):
            x = cp.Variable(m.p)
            cons = [D @ x <= 0, c @ x == t0] + ([E @ x == 0] if E.shape[0] else [])
            pr = cp.Problem(cp.Minimize(cp.sum_squares((y - A @ x) / s0)), cons)
            pr.solve(solver="OSQP", eps_abs=1e-11, eps_rel=1e-11, max_iter=400000, polish=True)
            fr = F.fit(c, t0)
            if pr.status == "optimal" and fr["ok"]:
                d_ = abs(pr.value - fr["T"]); worst = max(worst, d_); worst_rel = max(worst_rel, d_ / max(1.0, fr["T"])); nb += 1
say(f"(B) T(theta0), package vs CVXPY+OSQP on {nb} problems: max abs diff {worst:.2e}, max rel diff {worst_rel:.2e} ({time.time()-t:.0f} s)")
# ---------------- (C) and (D): exactly Gaussian cells with known variances
mu_sim = np.array([d["mu_true"] for d in info])
def mu_with(Bover):
    return np.array([W.true_mean(K, procs, {**over, **Bover}) for (_, K, procs, _, _, over) in W.CELLS])
configs = {"simulation truth": mu_sim, "least favourable (H=P=R=0)": mu_with(dict(H=0.0, P=0.0, R=0.0))}
for fn in ("Early-stage population", "Late minus early population"):
    cf = W.FUNCS[fn]; th = float(np.asarray(cf) @ W.F_TRUE)
    for lab, mu in configs.items():
        t = time.time(); R = 2000; rej05 = rej10 = 0
        for r in range(R):
            y = mu + s0 * rng.standard_normal(len(s0))
            p, T, rk = cc_pvalue(m, TRUE, cf, th, y, V0)
            rej05 += p < 0.05; rej10 += p < 0.10
        se = np.sqrt(0.05 * 0.95 / R)
        say(f"(C) size, {fn}, {lab}: rejection {rej05/R:.4f} at 0.05 (MC SE {se:.4f}), {rej10/R:.4f} at 0.10 ({time.time()-t:.0f} s)")
t = time.time(); R = 300; cov = cov_cc = miss = 0; islands = 0; dense_n = 0
cf = W.FUNCS["Early-stage population"]; th = float(np.asarray(cf) @ W.F_TRUE); c = m.cvec(cf)
for r in range(R):
    y = mu_sim + s0 * rng.standard_normal(len(s0))
    st = signet_set(m, TRUE, cf, 0.05, y, V0)
    lo, hi = st.get("lo"), st.get("hi")
    cov += (lo is not None and lo - 1e-9 <= th <= hi + 1e-9)
    cs = conf_set(m, TRUE, cf, 0.025, y, V0)
    inside = (cs.get("lo") is not None and cs["lo"] - 1e-9 <= th <= cs["hi"] + 1e-9)
    cov_cc += inside
    F = Fitter(m, TRUE, y, V0)
    if not inside and F.accept(c, th, 0.025)[0]:
        miss += 1
    if r < 25 and cs.get("lo") is not None and np.isfinite(cs["lo"]) and np.isfinite(cs["hi"]):
        dense_n += 1
        w = cs["hi"] - cs["lo"]
        grid = np.linspace(cs["lo"] - 3 * w, cs["hi"] + 3 * w, 300)
        acc_out = [g for g in grid if (g < cs["lo"] - 1e-6 or g > cs["hi"] + 1e-6) and F.accept(c, g, 0.025)[0]]
        islands += bool(acc_out)
say(f"(D) {R} literatures, true assumption set, Gaussian cells: SIGNET set coverage {cov/R:.3f}; "
    f"CC set (alpha/2) coverage {cov_cc/R:.3f}; theta accepted but outside reported CC set: {miss}; "
    f"dense scan (300 points, 3 widths each side) found accepted values outside the set in {islands}/{dense_n} literatures ({time.time()-t:.0f} s)")
# ---------------- (E)
t = time.time(); worst = 0.0; negnu = 0; ncert = 0
for r in range(150):
    inc = [g for g in groups if rng.random() < 0.6]
    cpd = m.compile(inc); A, D, E = cpd["A"], cpd["D"], cpd["E"]
    for fn, cf in W.FUNCS.items():
        for side in ("hi", "lo"):
            ct = certificate(m, inc, cf, side, V0)
            if ct is None: continue
            c = m.cvec(cf) * (1 if side == "hi" else -1)
            # recover nu, kappa for the returned lam by least squares on the remaining terms
            resid0 = c - A.T @ ct["lam"]
            M = np.hstack([D.T, E.T]) if E.shape[0] else D.T
            sol = np.linalg.lstsq(M, resid0, rcond=None)[0] if M.shape[1] else np.zeros(0)
            nu_ = sol[:D.shape[0]]
            # feasibility of nu >= 0 with exact equality (LP) as the strict check
            nu = cp.Variable(D.shape[0]); terms = D.T @ nu
            if E.shape[0]:
                kap = cp.Variable(E.shape[0]); terms = terms + E.T @ kap
            pr = cp.Problem(cp.Minimize(cp.norm1(terms - resid0)), [nu >= 0]); pr.solve(solver="CLARABEL")
            worst = max(worst, pr.value); ncert += 1
say(f"(E) {ncert} certificates: largest L1 residual of c - A'lam over D'nu + E'kap with nu >= 0: {worst:.2e} ({time.time()-t:.0f} s)")
open(os.path.join(ROOT, "audit", "independent_checks_output.txt"), "w").write("\n".join(out) + "\n")
