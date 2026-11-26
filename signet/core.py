"""
signet_core.py
Reference implementation of SIGNET (SIGNed Evidence Triangulation).

Model
-----
    y_i = K_i' f + sum_{p in P(i)} b[p,i] + e_i,     e ~ N(0, V), V known
    f       effect field over G contexts (G = 1 for a single effect)
    b[p,i]  contribution of bias process p to evidence cell i ("incidence")
    x = (f, b) is the full parameter; mu = A x is the vector of cell means.

Knowledge enters only as named assumptions, each a group of linear
constraints on x (homogeneous unless a magnitude bound is declared):
    sign      s * b[p,i] >= 0                    direction of a bias
    zero      b[p,i] = 0                         absence of a bias in a cell
    share     b[p,j] = r_j * b[p,i0]             equal (or proportional) size
    monotone  f_g <= f_{g+1}                     shape of the effect field
    fieldsign s * f_g >= 0                       mechanistic sign of the effect
    bound     lo <= b[p,i] <= hi                 optional magnitude rung

Inference on theta = c'x inverts the conditional chi-squared test of
Cox and Shi (2023) for inequalities in which nuisance parameters enter
linearly: T(theta0) = min (y - Ax)' V^-1 (y - Ax) over feasible x with
c'x = theta0, compared with a chi-squared quantile whose degrees of
freedom equal the rank of the active constraints in the space of cell
means. The rank is computed as n - rank(A N), where N spans the affine
hull of the minimal face containing the optimal set (its implicit
equalities are verified by linear programmes). Over-counting the rank only
makes the test conservative.
"""
import itertools
import numpy as np
import scipy.sparse as sp
from scipy import linalg as sla
from scipy import optimize as sopt
from scipy import stats as sst
import clarabel

INF = float("inf")

# Numerical thresholds (v1.2.0: named and adjustable; the defaults are the values hard-coded in v1.0.0
# and v1.1.0, so results are unchanged). They are absolute. T is dimensionless. The two slack
# thresholds are in the units of the parameters, so the defaults assume estimates on the log-ratio
# scale with standard errors of ordinary size (roughly 0.001 to 10); the sensitivity of every reported
# result to these values is in experiments/audit_v12 (tolerance_sensitivity.md).
TOL = dict(T_zero=1e-9,     # T(theta0) at or below this is treated as zero (p = 1, no active constraint)
           slack=1e-4,      # slack below slack * (1 + |d_j|) makes constraint j a candidate active constraint
           implicit=1e-7)   # a candidate whose largest attainable slack is below this is an implicit equality


# ----------------------------------------------------------------------
# Model specification
# ----------------------------------------------------------------------
class Model:
    def __init__(self, G=1, field_names=None):
        self.G = int(G)
        self.field_names = field_names or [f"f{g + 1}" for g in range(self.G)]
        self.rows = []          # dict(name, K, procs)
        self.inc = []           # (process, row index) pairs
        self.groups = []        # dict(name, rung, kind, cons=[(coef{j:v}, 'le'|'eq', rhs)])
        self.y = None
        self.V = None
        # v1.1.0: constraint rows are scaled to unit length when compiled, so that certificate
        # selection and the active-set tolerances do not depend on how a constraint is written
        # (for example b1 <= b2 versus 2 b1 <= 2 b2). Set to False to reproduce v1.0.0 exactly.
        self.normalise_rows = True

    # --- structure ----------------------------------------------------
    def add_row(self, name, K, procs=()):
        K = np.atleast_1d(np.asarray(K, float))
        assert K.size == self.G
        i = len(self.rows)
        self.rows.append(dict(name=name, K=K, procs=list(procs)))
        for p in procs:
            self.inc.append((p, i))
        return i

    def row(self, name):
        for i, r in enumerate(self.rows):
            if r["name"] == name:
                return i
        raise KeyError(name)

    @property
    def n(self):
        return len(self.rows)

    @property
    def p(self):
        return self.G + len(self.inc)

    def pidx(self, proc, row):
        i = self.row(row) if isinstance(row, str) else row
        return self.G + self.inc.index((proc, i))

    def rows_with(self, proc):
        return [i for (q, i) in self.inc if q == proc]

    def param_names(self):
        return self.field_names + [f"{q}@{self.rows[i]['name']}" for (q, i) in self.inc]

    def A(self):
        A = np.zeros((self.n, self.p))
        for i, r in enumerate(self.rows):
            A[i, :self.G] = r["K"]
        for m, (q, i) in enumerate(self.inc):
            A[i, self.G + m] = 1.0
        return A

    def set_data(self, y, V):
        self.y = np.asarray(y, float)
        V = np.asarray(V, float)
        self.V = np.diag(V ** 2) if V.ndim == 1 else V

    # --- assumptions ----------------------------------------------------
    def _add(self, name, rung, kind, cons):
        assert name not in [g["name"] for g in self.groups], name
        self.groups.append(dict(name=name, rung=rung, kind=kind, cons=cons))

    def sign(self, name, proc, s, rows=None, rung=1):
        rows = self.rows_with(proc) if rows is None else [self.row(r) if isinstance(r, str) else r for r in rows]
        cons = [({self.pidx(proc, i): -float(s)}, "le", 0.0) for i in rows]
        self._add(name, rung, "sign", cons)

    def zero(self, name, proc, rows=None, rung=1):
        rows = self.rows_with(proc) if rows is None else [self.row(r) if isinstance(r, str) else r for r in rows]
        cons = [({self.pidx(proc, i): 1.0}, "eq", 0.0) for i in rows]
        self._add(name, rung, "zero", cons)

    def share(self, name, proc, rows=None, ratios=None, rung=2):
        rows = self.rows_with(proc) if rows is None else [self.row(r) if isinstance(r, str) else r for r in rows]
        ratios = [1.0] * len(rows) if ratios is None else list(ratios)
        j0 = self.pidx(proc, rows[0])
        cons = []
        for i, r in zip(rows[1:], ratios[1:]):
            cons.append(({self.pidx(proc, i): 1.0, j0: -float(r) / ratios[0]}, "eq", 0.0))
        self._add(name, rung, "share", cons)

    def monotone(self, name, direction=1, rung=2):
        cons = []
        for g in range(self.G - 1):
            # direction +1: f_g <= f_{g+1}
            cons.append(({g: float(direction), g + 1: -float(direction)}, "le", 0.0))
        self._add(name, rung, "monotone", cons)

    def fieldsign(self, name, g, s, rung=2):
        self._add(name, rung, "fieldsign", [({g: -float(s)}, "le", 0.0)])

    def bound(self, name, proc, lo, hi, rows=None, rung=4):
        rows = self.rows_with(proc) if rows is None else [self.row(r) if isinstance(r, str) else r for r in rows]
        cons = []
        for i in rows:
            j = self.pidx(proc, i)
            cons.append(({j: 1.0}, "le", float(hi)))
            cons.append(({j: -1.0}, "le", -float(lo)))
        self._add(name, rung, "bound", cons)

    def order(self, name, proc, row_lo, row_hi, rung=2):
        """b[proc, row_lo] <= b[proc, row_hi] (for example: a negative-control row whose
        bias is at least as negative as the outcome row's)."""
        i = self.row(row_lo) if isinstance(row_lo, str) else row_lo
        j = self.row(row_hi) if isinstance(row_hi, str) else row_hi
        self._add(name, rung, "order", [({self.pidx(proc, i): 1.0, self.pidx(proc, j): -1.0}, "le", 0.0)])

    def linear(self, name, coefs, sense="le", rhs=0.0, rung=2):
        """Generic linear assumption sum_k coef_k * param_k (sense) rhs, with parameters named
        as in param_names(): field names, or 'proc@row' for bias contributions."""
        names = self.param_names()
        coef = {names.index(k): float(v) for k, v in coefs.items()}
        if sense == "ge":
            coef = {k: -v for k, v in coef.items()}; rhs = -float(rhs); sense = "le"
        assert sense in ("le", "eq")
        self._add(name, rung, "linear", [(coef, sense, float(rhs))])

    def group_names(self, max_rung=None):
        return [g["name"] for g in self.groups if max_rung is None or g["rung"] <= max_rung]

    def compile(self, include):
        include = set(include)
        Dl, dl, gl, El, el, ge = [], [], [], [], [], []
        for gi, g in enumerate(self.groups):
            if g["name"] not in include:
                continue
            for coef, typ, rhs in g["cons"]:
                row = np.zeros(self.p)
                for j, v in coef.items():
                    row[j] += v
                if typ == "le":
                    Dl.append(row); dl.append(rhs); gl.append(gi)
                else:
                    El.append(row); el.append(rhs); ge.append(gi)
        p = self.p
        D, d = np.array(Dl).reshape(-1, p), np.array(dl, float)
        E, e = np.array(El).reshape(-1, p), np.array(el, float)
        if getattr(self, "normalise_rows", True):
            for M_, v_ in ((D, d), (E, e)):
                if M_.shape[0]:
                    s = np.linalg.norm(M_, axis=1)
                    s[s == 0] = 1.0
                    M_ /= s[:, None]
                    v_ /= s
        return dict(A=self.A(), D=D, d=d, gD=np.array(gl, int), E=E, e=e, gE=np.array(ge, int))

    def cvec(self, cf):
        """Functional on the field (length G) or on all parameters (length p)."""
        cf = np.atleast_1d(np.asarray(cf, float))
        if cf.size == self.p:
            return cf
        assert cf.size == self.G
        return np.r_[cf, np.zeros(self.p - self.G)]


# ----------------------------------------------------------------------
# Convex programmes
# ----------------------------------------------------------------------
def _settings(tol=1e-9, equilibrate=False):
    st = clarabel.DefaultSettings()
    st.verbose = False
    st.equilibrate_enable = equilibrate
    st.tol_gap_abs = tol
    st.tol_gap_rel = tol
    st.tol_feas = tol
    st.max_iter = 400
    return st


_ST = _settings()
_ST_EQ = _settings(1e-9, True)
_ST_LOOSE = _settings(1e-7, True)


def _solve_ineq_qp(Pm, q, Dc, dc):
    """min 0.5 w'Pw + q'w  s.t.  Dc w <= dc  (no equalities: they are eliminated)."""
    k = Pm.shape[0]
    Ptri = sp.triu(sp.csc_matrix(Pm), format="csc")
    if Dc.shape[0]:
        Acon = sp.csc_matrix(Dc); b = np.asarray(dc, float)
        cones = [clarabel.NonnegativeConeT(Dc.shape[0])]
    else:
        Acon = sp.csc_matrix((0, k)); b = np.zeros(0); cones = []
    for st in (_ST, _ST_EQ, _ST_LOOSE):
        sol = clarabel.DefaultSolver(Ptri, q, Acon, b, cones, st).solve()
        status = str(sol.status)
        if status in ("Solved", "AlmostSolved"):
            return True, status, np.array(sol.x)
    # fallback: sequential quadratic programming from the conic solver's last iterate
    w0 = np.array(sol.x) if sol.x is not None and np.all(np.isfinite(sol.x)) else np.zeros(k)
    cons = [dict(type="ineq", fun=lambda w: b - Dc @ w, jac=lambda w: -Dc)] if Dc.shape[0] else []
    res = sopt.minimize(lambda w: 0.5 * w @ Pm @ w + q @ w, w0, jac=lambda w: Pm @ w + q,
                        constraints=cons, method="SLSQP", options=dict(maxiter=500, ftol=1e-14))
    if res.success and (not Dc.shape[0] or np.all(b - Dc @ res.x > -1e-7)):
        return True, "SLSQP", res.x
    return False, status, None


def _linprog(c, A_eq=None, b_eq=None, A_ub=None, b_ub=None):
    res = sopt.linprog(c, A_ub=A_ub if A_ub is not None and len(A_ub) else None,
                       b_ub=b_ub if b_ub is not None and len(b_ub) else None,
                       A_eq=A_eq if A_eq is not None and len(A_eq) else None,
                       b_eq=b_eq if b_eq is not None and len(b_eq) else None,
                       bounds=[(None, None)] * len(c), method="highs")
    return res


def _affine(E, e, p):
    """x = x0 + N z parametrises {x : E x = e}."""
    if E.shape[0] == 0:
        return np.zeros(p), np.eye(p)
    x0 = np.linalg.lstsq(E, e, rcond=None)[0]
    N = sla.null_space(E, rcond=1e-11)
    return x0, N


class Fitter:
    """Holds a compiled assumption set and data; evaluates T(theta0).
    Equality constraints are eliminated by a null-space parametrisation so
    that the conic solver sees only inequalities."""

    def __init__(self, model, include, y=None, V=None, exact_rank=True):
        self.m = model
        self.cp = model.compile(include)
        self.y = model.y if y is None else np.asarray(y, float)
        V = model.V if V is None else V
        self.W = np.linalg.inv(V)
        self.A = self.cp["A"]
        self.n, self.p = self.A.shape
        self.x0, self.N = _affine(self.cp["E"], self.cp["e"], self.p)
        self.exact_rank = exact_rank
        self._cache = {}

    def _space(self, c):
        """Precompute the parametrisation with the extra row c'x = theta0."""
        key = None if c is None else c.tobytes()
        if key in self._cache:
            return self._cache[key]
        if c is None:
            sp_ = dict(B=self.N, base=self.x0, dirn=np.zeros(self.p), degenerate=False)
        else:
            g = c @ self.N
            if np.linalg.norm(g) < 1e-10:
                sp_ = dict(B=self.N, base=self.x0, dirn=np.zeros(self.p), degenerate=True,
                           fixed=float(c @ self.x0))
            else:
                z0 = g / (g @ g)             # c'N z0 = 1
                N2 = sla.null_space(g[None, :], rcond=1e-11)
                sp_ = dict(B=self.N @ N2, base=self.x0 - self.N @ z0 * float(c @ self.x0),
                           dirn=self.N @ z0, degenerate=False)
        AB = self.A @ sp_["B"]
        P0 = 2 * AB.T @ self.W @ AB
        sp_["P"] = P0 + 1e-9 * max(1.0, np.trace(P0) / max(1, P0.shape[0])) * np.eye(AB.shape[1])
        sp_["AB"] = AB
        sp_["DB"] = self.cp["D"] @ sp_["B"]
        self._cache[key] = sp_
        return sp_

    def fit(self, c=None, theta0=None):
        sp_ = self._space(c)
        E, e = self.cp["E"], self.cp["e"]
        if c is not None:
            E = np.vstack([E, c[None, :]]); e = np.r_[e, theta0]
            if sp_["degenerate"]:
                if abs(theta0 - sp_["fixed"]) > 1e-9:
                    return dict(ok=False, T=INF, status="theta fixed by assumptions")
                xp = sp_["base"]
            else:
                xp = sp_["base"] + sp_["dirn"] * theta0
        else:
            xp = sp_["base"]
        r0 = self.y - self.A @ xp
        q = -2 * sp_["AB"].T @ self.W @ r0
        D, d = self.cp["D"], self.cp["d"]
        ok, status, w = _solve_ineq_qp(sp_["P"], q, sp_["DB"], d - D @ xp)
        if not ok:
            return dict(ok=False, T=INF, status=status)
        x = xp + sp_["B"] @ w
        r = self.y - self.A @ x
        T = float(r @ self.W @ r)
        slack = (d - D @ x) if D.shape[0] else np.zeros(0)
        return dict(ok=True, T=max(T, 0.0), x=x, slack=slack, mu=self.A @ x, E=E, e=e, status=status)

    def rank(self, fr):
        """Rank of the active constraints in the space of cell means."""
        D, d = self.cp["D"], self.cp["d"]
        E, e = fr["E"], fr["e"]
        cand = np.where(fr["slack"] < TOL["slack"] * (1 + np.abs(d)))[0] if D.shape[0] else np.zeros(0, int)
        implicit = []
        if len(cand):
            Aeq = np.vstack([self.A, E]) if E.shape[0] else self.A
            beq = np.r_[fr["mu"], e] if E.shape[0] else fr["mu"]
            if self.exact_rank:
                # constraint j is an implicit equality of the optimal face iff it
                # cannot be made slack while A x = mu_hat and all constraints hold
                for j in cand:
                    cc = np.r_[np.zeros(self.p), -1.0]
                    Aub = np.vstack([np.c_[D, np.zeros(len(d))], np.r_[D[j], 1.0][None, :],
                                     np.r_[np.zeros(self.p), 1.0][None, :]])
                    bub = np.r_[d, d[j], 1.0]
                    res = _linprog(cc, np.c_[Aeq, np.zeros(Aeq.shape[0])], beq, Aub, bub)
                    if res.status != 0 or -res.fun < TOL["implicit"]:
                        implicit.append(j)
            else:
                # one programme maximising the sum of capped slacks (conservative)
                k = len(cand)
                cc = np.r_[np.zeros(self.p), -np.ones(k)]
                T1 = np.eye(k)
                Aub = np.vstack([np.c_[D, np.zeros((len(d), k))],
                                 np.c_[D[cand], T1],
                                 np.c_[np.zeros((k, self.p)), T1]])
                bub = np.r_[d, d[cand], np.ones(k)]
                res = _linprog(cc, np.c_[Aeq, np.zeros((Aeq.shape[0], k))], beq, Aub, bub)
                if res.status != 0:
                    implicit = list(cand)
                else:
                    t = res.x[self.p:]
                    implicit = [j for j, tj in zip(cand, t) if tj < TOL["implicit"]]
        M = E if not implicit else (np.vstack([E, D[implicit]]) if E.shape[0] else D[implicit])
        N = sla.null_space(M, rcond=1e-10) if M.shape[0] else np.eye(self.p)
        if N.shape[1] == 0:
            return self.n
        sv = np.linalg.svd(self.A @ N, compute_uv=False)
        rk = int(np.sum(sv > 1e-9 * max(1.0, sv.max() if sv.size else 1.0)))
        return self.n - rk

    def accept(self, c, theta0, alpha=0.05):
        fr = self.fit(c, theta0)
        if not fr["ok"]:
            return False, INF, None
        T = fr["T"]
        if T <= TOL["T_zero"]:
            return True, T, 0
        if T <= sst.chi2.ppf(1 - alpha, 1):
            return True, T, None          # rank >= 1 whenever T > 0
        if T > sst.chi2.ppf(1 - alpha, self.n):
            return False, T, None
        r = self.rank(fr)
        if r == 0:
            return False, T, 0
        return bool(T <= sst.chi2.ppf(1 - alpha, r)), T, r


# ----------------------------------------------------------------------
# Capacity, identified sets and certificates (noise-free geometry)
# ----------------------------------------------------------------------
def capacity(model, include, cf):
    """Is theta = c'x bounded above / below on {x feasible: A x = mu}?
    By Farkas' lemma this depends only on the design and the assumptions:
    theta is bounded above iff no direction d with A d = 0 that satisfies the
    (homogenised) constraints raises c'd."""
    cp = model.compile(include)
    c = model.cvec(cf)
    out = {}
    for side, sgn in (("above", 1.0), ("below", -1.0)):
        Aeq = np.vstack([cp["A"], cp["E"]])
        beq = np.zeros(Aeq.shape[0])
        Aub = np.vstack([cp["D"], sgn * c[None, :]])
        bub = np.r_[np.zeros(cp["D"].shape[0]), 1.0]
        res = _linprog(-sgn * c, Aeq, beq, Aub, bub)
        out[side] = bool(res.status == 0 and -res.fun < 1e-7)
    return out


def identified(model, include, cf, mu):
    """Sharp bounds on c'x given exact cell means mu, with certificates:
    U = lam' mu + nu' d + kap' e, where lam weights the cells and nu, kap
    are the multipliers of the assumptions that the bound uses."""
    cp = model.compile(include)
    c = model.cvec(cf)
    Aeq = np.vstack([cp["A"], cp["E"]]) if cp["E"].shape[0] else cp["A"]
    beq = np.r_[mu, cp["e"]] if cp["E"].shape[0] else np.asarray(mu, float)
    out = dict(feasible=True)
    for side, sgn in (("hi", 1.0), ("lo", -1.0)):
        res = _linprog(-sgn * c, Aeq, beq, cp["D"], cp["d"])
        if res.status == 2:
            out["feasible"] = False
            out[side] = None
            continue
        if res.status == 3:
            out[side] = sgn * INF
            out["cert_" + side] = None
            continue
        val = sgn * (-res.fun)
        out[side] = float(val)
        lam_all = -sgn * np.asarray(res.eqlin.marginals)   # d(bound)/d(rhs)
        lam = lam_all[:model.n]
        kap = lam_all[model.n:]
        nu = -sgn * np.asarray(res.ineqlin.marginals) if cp["D"].shape[0] else np.zeros(0)
        used = {}
        for j, gi in enumerate(cp["gD"]):
            if abs(nu[j]) > 1e-9:
                used[model.groups[gi]["name"]] = used.get(model.groups[gi]["name"], 0) + abs(nu[j])
        for j, gi in enumerate(cp["gE"]):
            if abs(kap[j]) > 1e-9:
                used[model.groups[gi]["name"]] = used.get(model.groups[gi]["name"], 0) + abs(kap[j])
        out["cert_" + side] = dict(
            weights={model.rows[i]["name"]: float(lam[i]) for i in range(model.n) if abs(lam[i]) > 1e-9},
            assumptions=sorted(used))
    return out


# ----------------------------------------------------------------------
# Confidence sets by test inversion
# ----------------------------------------------------------------------
def conf_set(model, include, cf, alpha=0.05, y=None, V=None, ngrid=40,
             exact_rank=True, span=6.0, tol=1e-4, outer_scan=16, inversion="hull", max_span=1e4):
    """Confidence set for theta = c'x by inverting the conditional chi-squared test.

    T(theta0) is convex in theta0 (partial minimisation of a convex quadratic over a convex set),
    and every value the test accepts satisfies T <= chi2_n(1 - alpha), because the rank-adjusted
    critical value never exceeds the n-degree one. The outer set {theta0 : T <= chi2_n(1 - alpha)}
    is therefore an interval that contains every accepted value; its end points are found by
    bracketing and bisection. Under the model it is itself a valid set, because
    T(theta) <= (y - mu)' V^-1 (y - mu), which is chi-squared with n degrees of freedom.

    inversion="hull" (default) reports the convex hull of the accepted values: a grid of ngrid
    points spans the outer set, the hull is refined by bisection, and outer_scan extra points on
    each side between the hull and the outer limits are checked, the hull being extended to any
    accepted value found. inversion="outer" reports the outer set (guaranteed conservative,
    wider). v1.0.0 behaviour is reproduced by outer_scan=0, except that v1.0.0 stopped the
    bracketing search at span; v1.1.0 continues to max_span and reports an infinite end point,
    with a warning, if the outer limit is not bracketed."""
    c = model.cvec(cf)
    F = Fitter(model, include, y, V, exact_rank)
    cmax = sst.chi2.ppf(1 - alpha, F.n)
    base = F.fit()
    if not base["ok"]:
        return dict(lo=None, hi=None, empty=True, reason="infeasible")
    th0 = float(c @ base["x"])
    cap = capacity(model, include, c)
    warnings = []

    def acc(t):
        return F.accept(c, t, alpha)[0]

    def T_of(t):
        return F.fit(c, t)["T"]

    # first doubling step at or beyond span (where v1.0.0 stopped searching)
    s_legacy = 0.05
    while s_legacy < span:
        s_legacy *= 2
    lims, bracketed, legacy_trunc = {}, {}, {}
    for side, sgn in (("hi", 1.0), ("lo", -1.0)):
        if not cap["above" if sgn > 0 else "below"]:
            lims[side] = th0 + sgn * span
            bracketed[side] = None              # unbounded by capacity: no outer limit exists
            legacy_trunc[side] = False
            continue
        step = 0.05
        t_in = th0
        t_out = th0 + sgn * step
        while T_of(t_out) <= cmax and abs(t_out - th0) < max_span:
            t_in, step = t_out, step * 2
            t_out = th0 + sgn * step
        bracketed[side] = bool(T_of(t_out) > cmax)
        if not bracketed[side]:
            warnings.append(f"outer limit ({side}) not bracketed within {max_span:g} of the estimate; "
                            "the set is reported as unbounded on that side")
        for _ in range(60):
            mid = 0.5 * (t_in + t_out)
            if T_of(mid) <= cmax:
                t_in = mid
            else:
                t_out = mid
            if abs(t_out - t_in) < tol:
                break
        lims[side] = t_out
        legacy_trunc[side] = bool(abs(t_out - th0) > s_legacy)
    outer = dict(lo=-INF if (not cap["below"] or bracketed["lo"] is False) else lims["lo"],
                 hi=INF if (not cap["above"] or bracketed["hi"] is False) else lims["hi"])
    info = dict(th0=th0, T0=base["T"], capacity=cap, outer=outer, bracketed=bracketed,
                legacy_truncation=legacy_trunc, warnings=warnings, inversion=inversion)
    if inversion == "outer":
        if base["T"] > cmax:
            return dict(lo=None, hi=None, empty=True, reason="rejected", **info)
        return dict(lo=float(outer["lo"]), hi=float(outer["hi"]), empty=False, hull_extended=False, **info)
    if inversion != "hull":
        raise ValueError("inversion must be 'hull' or 'outer'")
    grid = np.unique(np.r_[np.linspace(lims["lo"], lims["hi"], ngrid), th0])
    accs = np.array([acc(t) for t in grid])
    if not accs.any():
        return dict(lo=None, hi=None, empty=True, reason="rejected", **info)
    ia = np.where(accs)[0]

    def refine(t_acc, t_rej):
        for _ in range(40):
            mid = 0.5 * (t_acc + t_rej)
            if acc(mid):
                t_acc = mid
            else:
                t_rej = mid
            if abs(t_rej - t_acc) < tol:
                break
        return t_acc

    i_lo, i_hi = ia[0], ia[-1]
    if (i_lo == 0 and not cap["below"]) or bracketed["lo"] is False:
        lo = -INF
    else:
        lo = refine(grid[i_lo], grid[i_lo - 1]) if i_lo > 0 else grid[0]
    if (i_hi == len(grid) - 1 and not cap["above"]) or bracketed["hi"] is False:
        hi = INF
    else:
        hi = refine(grid[i_hi], grid[i_hi + 1]) if i_hi < len(grid) - 1 else grid[-1]
    # dense check (v1.1.0): accepted values strictly between the refined hull and the outer limits
    extended = False
    for side in (("lo", "hi") if outer_scan else ()):
        e = lo if side == "lo" else hi
        if not np.isfinite(e):
            continue
        pts = np.linspace(lims[side], e, outer_scan + 2)[1:-1]
        hit = [t for t in pts if acc(t)]
        if not hit:
            continue
        extended = True
        if side == "lo":
            t_new = min(hit)
            beyond = [t for t in pts if t < t_new]
            lo = refine(t_new, max(beyond) if beyond else lims["lo"])
        else:
            t_new = max(hit)
            beyond = [t for t in pts if t > t_new]
            hi = refine(t_new, min(beyond) if beyond else lims["hi"])
    return dict(lo=float(lo), hi=float(hi), empty=False, hull_extended=extended, **info)


def certificate(model, include, cf, side="hi", V=None):
    """Certificate for one side of theta = c'x, chosen without looking at the
    estimates y (the second stage uses V; when V is itself estimated from the same studies the
    choice depends on the data through V unless the certificate is unique, see certificate_range,
    or the weights are chosen with a fixed covariance, see cert_bounds(cert_V=...)).
    A certificate is a weight vector lam over cells with
    c'x <= lam'mu (side 'hi') for every x allowed by the assumptions; by
    Farkas' lemma it exists iff c = A'lam + D'nu + E'kap with nu >= 0.
    Choice rule (lexicographic): first minimise the total weight sum(nu)
    placed on inequality (sign, shape) assumptions, so that bounds resting
    on exact equalities are preferred because they carry no slack; then,
    among those, minimise the variance lam'V lam. Returns None when no
    certificate exists, i.e. the side is unbounded (cone criterion)."""
    cp = model.compile(include)
    if np.any(np.abs(cp["d"]) > 0) or np.any(np.abs(cp["e"]) > 0):
        raise ValueError("certificates are defined here for homogeneous assumptions only")
    c = model.cvec(cf) * (1.0 if side == "hi" else -1.0)
    V = model.V if V is None else V
    A, D, E = cp["A"], cp["D"], cp["E"]
    n, mD, mE = A.shape[0], D.shape[0], E.shape[0]
    k = n + mD + mE
    Aeq = np.hstack([A.T, D.T, E.T])                 # A'lam + D'nu + E'kap = c
    # step 1: least reliance on inequalities
    cost = np.r_[np.zeros(n), np.ones(mD), np.zeros(mE)]
    bounds = [(None, None)] * n + [(0, None)] * mD + [(None, None)] * mE
    res = sopt.linprog(cost, A_eq=Aeq, b_eq=c, bounds=bounds, method="highs")
    if res.status != 0:
        return None
    s1 = float(res.fun)
    # step 2: least variance among those
    Pm = np.zeros((k, k)); Pm[:n, :n] = 2 * V
    Pm += 1e-12 * np.eye(k)
    blocks = [Aeq]; rhs = [c]; cones = [clarabel.ZeroConeT(Aeq.shape[0])]
    if mD:
        G = np.zeros((mD + 1, k)); G[:mD, n:n + mD] = -np.eye(mD); G[mD, n:n + mD] = 1.0
        blocks.append(G); rhs.append(np.r_[np.zeros(mD), s1 + 1e-9 * (1 + s1)])
        cones.append(clarabel.NonnegativeConeT(mD + 1))
    Acon = sp.csc_matrix(np.vstack(blocks)); b = np.concatenate(rhs)
    sol = None
    for st in (_ST, _ST_EQ, _ST_LOOSE):
        sol = clarabel.DefaultSolver(sp.triu(sp.csc_matrix(Pm), format="csc"), np.zeros(k), Acon, b, cones, st).solve()
        if str(sol.status) in ("Solved", "AlmostSolved"):
            break
    if str(sol.status) not in ("Solved", "AlmostSolved"):
        z = res.x                                # fall back to the step-1 solution
    else:
        z = np.array(sol.x)
    lam, nu, kap = z[:n], z[n:n + mD], z[n + mD:]
    used = {}
    for j, gi in enumerate(cp["gD"]):
        if nu[j] > 1e-7:
            used[model.groups[gi]["name"]] = True
    for j, gi in enumerate(cp["gE"]):
        if abs(kap[j]) > 1e-7:
            used[model.groups[gi]["name"]] = True
    return dict(lam=lam, sd=float(np.sqrt(max(lam @ V @ lam, 0.0))), side=side, slack_weight=s1,
                weights={model.rows[i]["name"]: float(lam[i] * (1 if side == "hi" else -1))
                         for i in range(n) if abs(lam[i]) > 1e-6},
                assumptions=sorted(used))


def check_certificate(model, include, cf, lam, side="hi"):
    """Dual feasibility of a given cell-weight vector lam (for example one fixed in a protocol).
    lam certifies side 'hi' when c = A'lam + D'nu + E'kap for some nu >= 0. With bounds or other
    non-homogeneous assumptions (D x <= d, E x = e) the implied bound is
        c'x <= lam'mu + nu'd + kap'e,
    and the smallest offset nu'd + kap'e over the feasible multipliers is returned."""
    cp = model.compile(include)
    c = model.cvec(cf) * (1.0 if side == "hi" else -1.0)
    lam = np.asarray(lam, float) * (1.0 if side == "hi" else -1.0)
    A, D, E = cp["A"], cp["D"], cp["E"]
    mD, mE = D.shape[0], E.shape[0]
    rhs = c - A.T @ lam
    if mD + mE == 0:
        ok = bool(np.linalg.norm(rhs) < 1e-8)
        return dict(valid=ok, offset=0.0 if ok else None, residual=float(np.linalg.norm(rhs)))
    Aeq = np.hstack([D.T, E.T])
    cost = np.r_[cp["d"], cp["e"]]
    bounds = [(0, None)] * mD + [(None, None)] * mE
    res = sopt.linprog(cost, A_eq=Aeq, b_eq=rhs, bounds=bounds, method="highs")
    if res.status == 0:
        return dict(valid=True, offset=float(res.fun), residual=float(np.linalg.norm(Aeq @ res.x - rhs)))
    if res.status == 3:                          # feasible, offset unbounded below (cannot occur for a sound model)
        return dict(valid=True, offset=-INF, residual=0.0)
    return dict(valid=False, offset=None, residual=None)


def certificate_range(model, include, cf, side="hi", tol=1e-7):
    """Range of each cell weight over all certificates that place the least total weight on
    inequality assumptions (the first, data-free stage of the choice rule). When every range has
    zero width the certificate is unique: its choice cannot depend on V or on the estimates, so
    the certificate bound is an exact z-bound whenever V is known. Returns None when no
    certificate exists (v1.1.0)."""
    cp = model.compile(include)
    c = model.cvec(cf) * (1.0 if side == "hi" else -1.0)
    A, D, E = cp["A"], cp["D"], cp["E"]
    n, mD, mE = A.shape[0], D.shape[0], E.shape[0]
    Aeq = np.hstack([A.T, D.T, E.T])
    cost = np.r_[np.zeros(n), np.ones(mD), np.zeros(mE)]
    bounds = [(None, None)] * n + [(0, None)] * mD + [(None, None)] * mE
    res = sopt.linprog(cost, A_eq=Aeq, b_eq=c, bounds=bounds, method="highs")
    if res.status != 0:
        return None
    s1 = float(res.fun)
    Aub = cost[None, :] if mD else None
    bub = np.array([s1 + 1e-9 * (1 + s1)]) if mD else None
    lo, hi = np.zeros(n), np.zeros(n)
    for i in range(n):
        for sgn in (1.0, -1.0):
            obj = np.zeros(n + mD + mE); obj[i] = sgn
            r = sopt.linprog(obj, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=c, bounds=bounds, method="highs")
            v = sgn * INF * -1 if r.status == 3 else (sgn * r.fun if r.status == 0 else np.nan)
            if sgn > 0:
                lo[i] = v
            else:
                hi[i] = v
    width = hi - lo
    unique = bool(np.all(np.isfinite(width)) and np.all(width < tol * (1 + np.abs(lo))))
    return dict(unique=unique, lo=lo, hi=hi, slack_weight=s1)


def cert_bounds(model, include, cf, a=0.0125, y=None, V=None, cert_V=None):
    """One-sided bounds from the two minimum-variance certificates, each at level a.
    With homogeneous assumptions they are exact z-bounds when the certificate weights do not
    depend on y. cert_V (v1.1.0), when given, is the covariance used to choose the weights,
    for example one computed from the within-study standard errors only, or fixed in a protocol;
    the standard deviation of the certified combination is always computed from V."""
    y = model.y if y is None else np.asarray(y, float)
    V = model.V if V is None else V
    z = sst.norm.ppf(1 - a)
    out = {}
    for side in ("hi", "lo"):
        ct = certificate(model, include, cf, side, V if cert_V is None else cert_V)
        if ct is None:
            out[side] = INF if side == "hi" else -INF
            out["cert_" + side] = None
        else:
            if cert_V is not None:
                ct["sd"] = float(np.sqrt(max(ct["lam"] @ V @ ct["lam"], 0.0)))
                ct["selected_with"] = "cert_V"
            est = float(ct["lam"] @ y)
            out[side] = est + z * ct["sd"] if side == "hi" else -(est + z * ct["sd"])
            out["cert_" + side] = ct
            out["est_" + side] = est if side == "hi" else -est
    return out


def signet_set(model, include, cf, alpha=0.05, y=None, V=None, cert_V=None, **kw):
    """SIGNET confidence set: the minimum-variance certificate bounds (level
    alpha/4 per side) intersected with the conditional chi-squared set
    (level alpha/2). Valid at level alpha by Bonferroni when V is known and the
    certificate weights do not depend on y (see cert_bounds for cert_V)."""
    cb = cert_bounds(model, include, cf, alpha / 4, y, V, cert_V)
    cs = conf_set(model, include, cf, alpha / 2, y, V, **kw)
    if cs["empty"]:
        return dict(lo=None, hi=None, empty=True, reason=cs.get("reason"), cert=cb, cc=cs)
    lo, hi = max(cb["lo"], cs["lo"]), min(cb["hi"], cs["hi"])
    if lo > hi:
        return dict(lo=None, hi=None, empty=True, reason="certificate and test disagree", cert=cb, cc=cs)
    return dict(lo=lo, hi=hi, empty=False, cert=cb, cc=cs)


def box_set(model, include, cf, alpha=0.05, y=None, V=None):
    """Version-1 inference: LP bounds over a Bonferroni box for the cell means
    (valid but conservative)."""
    cp = model.compile(include)
    c = model.cvec(cf)
    y = model.y if y is None else np.asarray(y, float)
    V = model.V if V is None else V
    s = np.sqrt(np.diag(V))
    z = sst.norm.ppf(1 - alpha / (2 * len(y)))
    A = cp["A"]
    Aub = np.vstack([cp["D"], A, -A])
    bub = np.r_[cp["d"], y + z * s, -(y - z * s)]
    out = {}
    for side, sgn in (("hi", 1.0), ("lo", -1.0)):
        res = _linprog(-sgn * c, cp["E"], cp["e"], Aub, bub)
        if res.status == 2:
            return dict(lo=None, hi=None, empty=True)
        out[side] = sgn * INF if res.status == 3 else float(sgn * (-res.fun))
    out["empty"] = False
    return out


def model_test(model, include, alpha=0.05, y=None, V=None, exact_rank=True):
    """Conditional chi-squared test of the assumption set itself."""
    F = Fitter(model, include, y, V, exact_rank)
    fr = F.fit()
    if not fr["ok"]:
        return dict(reject=True, T=INF, r=None, crit=None)
    T = fr["T"]
    r = 0 if T <= TOL["T_zero"] else F.rank(fr)
    crit = sst.chi2.ppf(1 - alpha, r) if r > 0 else 0.0
    pval = float(sst.chi2.sf(T, r)) if r > 0 else 1.0
    return dict(reject=bool(T > crit + 1e-9), T=float(T), r=int(r), crit=float(crit), p=pval)


def repairs(model, include, alpha=0.05, kmax=2, droppable=None, y=None, V=None):
    """Minimal sets of assumption groups whose removal stops rejection
    (the minimal hitting sets of the rejected conflicts)."""
    include = list(include)
    droppable = include if droppable is None else [g for g in droppable if g in include]
    found = []
    for k in range(0, kmax + 1):
        for S in itertools.combinations(droppable, k):
            if any(set(f) <= set(S) for f in found):
                continue
            keep = [g for g in include if g not in S]
            mt = model_test(model, keep, alpha, y, V)
            if not mt["reject"]:
                found.append(S)
        if found and k == 0:
            return [()]
    return found


def tolerant_set(model, include, cf, k=1, alpha=0.05, droppable=None, y=None, V=None,
                 setfun=None, **kw):
    """Union of confidence sets over all ways of dropping up to k groups."""
    setfun = signet_set if setfun is None else setfun
    include = list(include)
    droppable = include if droppable is None else [g for g in droppable if g in include]
    lo, hi, any_ok = INF, -INF, False
    for kk in range(k + 1):
        for S in itertools.combinations(droppable, kk):
            keep = [g for g in include if g not in S]
            cs = setfun(model, keep, cf, alpha, y, V, **kw)
            if cs["empty"]:
                continue
            any_ok = True
            lo, hi = min(lo, cs["lo"]), max(hi, cs["hi"])
    if not any_ok:
        return dict(lo=None, hi=None, empty=True)
    return dict(lo=lo, hi=hi, empty=False)


class _OmitRows:
    """View of a Model with some evidence cells omitted: their estimates leave the analysis, while their
    bias parameters stay in the model, constrained only by the assumptions (added in the GENESIS VI audit)."""

    def __init__(self, base, omit):
        self._b = base
        omit = set(base.row(o) if isinstance(o, str) else int(o) for o in omit)
        self._keep = [i for i in range(base.n) if i not in omit]
        self.rows = [base.rows[i] for i in self._keep]
        self.y = None if base.y is None else np.asarray(base.y, float)[self._keep]
        self.V = None if base.V is None else np.asarray(base.V, float)[np.ix_(self._keep, self._keep)]

    def __getattr__(self, k):
        return getattr(self._b, k)

    @property
    def n(self):
        return len(self._keep)

    def A(self):
        return self._b.A()[self._keep]

    def compile(self, include):
        out = dict(self._b.compile(include))
        out["A"] = out["A"][self._keep]
        return out


def _union(sets):
    ok = [s for s in sets if not s.get("empty")]
    if not ok:
        return dict(lo=None, hi=None, empty=True)
    return dict(lo=min(s["lo"] for s in ok), hi=max(s["hi"] for s in ok), empty=False)


def cell_tolerant_set(model, include, cf, alpha=0.05, y=None, V=None, setfun=None, **kw):
    """Union of the full analysis and the analyses that each omit one evidence cell (Proposition S4).
    It covers theta with probability at least 1 - alpha when at most one cell is broken and the named
    assumptions hold for the others. When omitting a cell leaves nothing that bounds theta, that
    analysis is unbounded and so is the union, as the proposition requires."""
    setfun = signet_set if setfun is None else setfun
    y = model.y if y is None else np.asarray(y, float)
    V = model.V if V is None else np.asarray(V, float)
    sets = [setfun(model, include, cf, alpha, y, V, **kw)]
    for i in range(model.n):
        keep = [j for j in range(model.n) if j != i]
        if not keep:
            return dict(lo=-INF, hi=INF, empty=False)
        kw_i = dict(kw)
        if kw.get("cert_V") is not None:
            kw_i["cert_V"] = np.asarray(kw["cert_V"], float)[np.ix_(keep, keep)]
        sets.append(setfun(_OmitRows(model, [i]), include, cf, alpha, y[keep], V[np.ix_(keep, keep)], **kw_i))
    return _union(sets)


def default_set(model, include, cf, alpha=0.05, droppable=None, y=None, V=None, **kw):
    """The reporting default: tolerant to one failed assumption group or one broken evidence cell."""
    g = tolerant_set(model, include, cf, k=1, alpha=alpha, droppable=droppable, y=y, V=V, **kw)
    c = cell_tolerant_set(model, include, cf, alpha=alpha, y=y, V=V, **kw)
    return _union([g, c])


def default_capacity(model, include, cf, droppable=None):
    """Capacity of the reporting default (v1.2.0). The default is a union of analyses, so it can be
    bounded on a side only if every analysis in the union is: the full set of groups, each set that omits
    one droppable group, and each analysis that omits one evidence cell. Like capacity(), this depends on
    the design, the assumptions and the target, not on the estimates, so it can be computed and registered
    before extraction. When one cell alone measures the target, omitting it leaves nothing to bound the
    target and the default is unbounded on both sides whatever the data.
    Returns dict(above, below, limiting), where limiting lists what removes each bound."""
    include = list(include)
    droppable = include if droppable is None else [g for g in droppable if g in include]
    out = dict(above=True, below=True, limiting=dict(above=[], below=[]))

    def note(cap, label):
        for side in ("above", "below"):
            if not cap[side]:
                out[side] = False
                out["limiting"][side].append(label)

    note(capacity(model, include, cf), "all groups in force")
    for g in droppable:
        note(capacity(model, [h for h in include if h != g], cf), f"without group {g}")
    for i in range(model.n):
        if model.n == 1:
            note(dict(above=False, below=False), f"without cell {model.rows[i]['name']}")
        else:
            note(capacity(_OmitRows(model, [i]), include, cf), f"without cell {model.rows[i]['name']}")
    return out


def pool_cell(y, s):
    """Random-effects (DerSimonian-Laird) pooling within one evidence cell."""
    y = np.asarray(y, float); s = np.asarray(s, float)
    w = 1 / s ** 2
    mu = np.sum(w * y) / np.sum(w)
    k = len(y)
    Q = float(np.sum(w * (y - mu) ** 2))
    tau2 = 0.0
    if k > 1:
        C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
        tau2 = max(0.0, (Q - (k - 1)) / C)
    ws = 1 / (s ** 2 + tau2)
    m = float(np.sum(ws * y) / np.sum(ws))
    return m, float(np.sqrt(1 / np.sum(ws))), dict(k=k, Q=Q, tau2=tau2)


# ----------------------------------------------------------------------
# p-value for a directional claim, consistent with signet_set
# ----------------------------------------------------------------------
def cc_pvalue(model, include, cf, theta0, y=None, V=None, exact_rank=True):
    """Conditional chi-squared p-value for H0: c'x = theta0 (Cox and Shi, 2023)."""
    c = model.cvec(cf)
    F = Fitter(model, include, y, V, exact_rank)
    fr = F.fit(c, theta0)
    if not fr["ok"]:
        return 0.0, INF, None
    T = fr["T"]
    if T <= TOL["T_zero"]:
        return 1.0, T, 0
    r = F.rank(fr)
    if r == 0:
        return 0.0, T, 0
    return float(sst.chi2.sf(T, r)), T, int(r)


def claim_pvalue(model, include, cf, claim="benefit", y=None, V=None, span=3.0, ngrid=61, cert_V=None,
                 alpha=0.05, empty=None):
    """p-value for a directional claim at the level used by signet_set.
    claim='benefit' tests H0: theta >= 0; claim='harm' tests H0: theta <= 0.
    A non-empty SIGNET set certifies the claim at level alpha exactly when
    min(4 p_cert, 2 p_cc) < alpha, because it intersects a certificate bound
    at level alpha/4 per side with the conditional chi-squared set at alpha/2.

    Empty sets (v1.2.0). When the SIGNET set at level alpha is empty, the data reject the assumption
    set (every value of theta is rejected, or the certificate bounds and the test disagree). The
    formula would then be small for benefit and for harm at once, so no directional p-value is
    reported: 'p' and 'p_upper' are None, 'set_empty' is True and 'reason' says why; the components
    p_cert, p_cc and p_cc_upper are still returned for the record. Pass empty=True or False when the
    caller has already computed the set at level alpha with its own inversion settings; with
    empty=None the set is computed here with the default settings of signet_set.

    p_cc is the supremum of the conditional chi-squared p-value over the null. The grid value
    (ngrid points on [0, span], plus the null value closest to the estimate) is a lower
    estimate of that supremum. v1.1.0 adds a validated upper bound: the p-value at theta0 never
    exceeds P(chi2_n > T(theta0)), and T is convex, so the supremum over the null is at most
    P(chi2_n > T(t*)), where t* is the null value nearest the unrestricted minimiser.
    p_upper = min(1, 4 p_cert, 2 p_cc_upper) is therefore a valid directional p-value."""
    y = model.y if y is None else np.asarray(y, float)
    V = model.V if V is None else V
    side = "hi" if claim == "benefit" else "lo"
    ct = certificate(model, include, cf, side, V if cert_V is None else cert_V)
    if ct is not None and cert_V is not None:
        ct["sd"] = float(np.sqrt(max(ct["lam"] @ V @ ct["lam"], 0.0)))
    if ct is None or ct["sd"] <= 0:
        p_cert = 1.0
    else:
        p_cert = float(sst.norm.cdf(float(ct["lam"] @ y) / ct["sd"]))
    c = model.cvec(cf)
    F = Fitter(model, include, y, V, True)
    base = F.fit()
    if not base["ok"]:
        t_star, p_cc_upper = 0.0, 0.0
    else:
        th0 = float(c @ base["x"])
        t_star = max(0.0, th0) if claim == "benefit" else min(0.0, th0)
        T_star = F.fit(c, t_star)["T"]
        p_cc_upper = 1.0 if T_star <= TOL["T_zero"] else float(sst.chi2.sf(T_star, F.n))
    grid = np.linspace(0.0, span, ngrid) if claim == "benefit" else np.linspace(-span, 0.0, ngrid)
    grid = np.unique(np.r_[grid, t_star])
    p_cc = max(cc_pvalue(model, include, cf, t, y, V)[0] for t in grid)
    p_cc_upper = max(p_cc_upper, p_cc)
    reason = None
    if empty is None:
        s_ = signet_set(model, include, cf, alpha, y, V, cert_V=cert_V)
        empty, reason = bool(s_["empty"]), s_.get("reason")
    cert_out = None if ct is None else dict(weights=ct["weights"], assumptions=ct["assumptions"], sd=ct["sd"])
    if empty:
        return dict(p_cert=p_cert, p_cc=p_cc, p=None, p_cc_upper=p_cc_upper, p_upper=None, set_empty=True,
                    reason=reason or "the SIGNET set is empty",
                    note="no directional p-value is reported for an empty set: the data reject the assumption set",
                    certificate=cert_out)
    return dict(p_cert=p_cert, p_cc=p_cc, p=float(min(1.0, 4 * p_cert, 2 * p_cc)),
                p_cc_upper=p_cc_upper, p_upper=float(min(1.0, 4 * p_cert, 2 * p_cc_upper)), set_empty=False,
                certificate=cert_out)
