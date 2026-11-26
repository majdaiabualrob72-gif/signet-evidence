"""
Specification layer: a JSON document describing evidence cells, assumptions and targets is
validated and compiled into a signet.core.Model.  Running the specification produces a
results dictionary that records capacity, SIGNET sets per rung, tolerance sets,
certificates, the model test, claim p-values and an input hash.

Specification fields
--------------------
name         free text
field        list of context names (length G); ["theta"] for a single effect
cells        list of {id, estimate, lower, upper, scale ("ratio" | "log"), se (optional),
                      kernel (length G), processes (list of names), design, source}
assumptions  list of {name, type, rung, source, ...type-specific fields}
             sign:      process, direction ("negative" | "positive"), rows (optional)
             zero:      process, rows (optional)
             share:     process, rows (optional), ratios (optional)
             order:     process, row_lo, row_hi            (b[row_lo] <= b[row_hi])
             fieldsign: field, direction
             monotone:  direction ("increasing" | "decreasing")
             bound:     process, lo, hi, rows (optional)
             linear:    coefs {parameter name: coefficient}, sense ("le" | "ge" | "eq"), rhs
certificate_se {cell id: standard error on the log scale}  (optional, v1.1.0; expected standard
             errors recorded before extraction; certificate weights are then chosen with these, so the
             choice cannot depend on the estimates; every cell must be listed)
correlations list of {cells: [id, id], rho}  (optional, v1.1.0; correlation between the
             sampling errors of two cells, for example estimates from the same trials; cells are
             independent unless listed)
targets      list of {name, coefficients {field name: weight}}
analysis     {alpha, rungs, tolerance, claims, inversion ("hull" | "outer")}
"""
import hashlib
import json

import numpy as np
from scipy import stats as sst

from .core import (Model, capacity, signet_set, conf_set, tolerant_set, default_set, model_test, repairs,
                   claim_pvalue, identified, certificate_range, default_capacity, INF)

Z975 = float(sst.norm.ppf(0.975))
TYPES = {"sign", "zero", "share", "order", "fieldsign", "monotone", "bound", "linear"}


class SpecError(ValueError):
    pass


def canonical_hash(spec):
    return hashlib.sha256(json.dumps(spec, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def validate(spec):
    """Quality control of the input; raises SpecError with every problem found."""
    errs = []
    field = spec.get("field") or ["theta"]
    G = len(field)
    ids = [c.get("id") for c in spec.get("cells", [])]
    if not ids:
        errs.append("no cells")
    if len(set(ids)) != len(ids):
        errs.append("cell ids are not unique")
    procs = set()
    for c in spec.get("cells", []):
        cid = c.get("id")
        k = c.get("kernel", [1.0] * G)
        if len(k) != G:
            errs.append(f"cell {cid}: kernel has length {len(k)}, field has {G}")
        scale = c.get("scale", "ratio")
        est = c.get("estimate")
        if est is None:
            errs.append(f"cell {cid}: missing estimate")
            continue
        if scale == "ratio":
            if est <= 0:
                errs.append(f"cell {cid}: ratio estimate must be positive")
            if "se" not in c:
                lo, hi = c.get("lower"), c.get("upper")
                if lo is None or hi is None:
                    errs.append(f"cell {cid}: give se or both confidence limits")
                elif not (0 < lo < est < hi):
                    errs.append(f"cell {cid}: limits must satisfy 0 < lower < estimate < upper")
        elif scale != "log":
            errs.append(f"cell {cid}: scale must be 'ratio' or 'log'")
        if "se" in c and not (c["se"] > 0):
            errs.append(f"cell {cid}: se must be positive")
        procs |= set(c.get("processes", []))
    names = [a.get("name") for a in spec.get("assumptions", [])]
    if len(set(names)) != len(names):
        errs.append("assumption names are not unique")
    for a in spec.get("assumptions", []):
        t = a.get("type")
        if t not in TYPES:
            errs.append(f"assumption {a.get('name')}: unknown type {t}")
            continue
        if t in ("sign", "zero", "share", "order", "bound") and a.get("process") not in procs:
            errs.append(f"assumption {a.get('name')}: process {a.get('process')} is not attached to any cell")
        for r in list(a.get("rows", []) or []) + [a.get("row_lo"), a.get("row_hi")]:
            if r is not None and r not in ids:
                errs.append(f"assumption {a.get('name')}: unknown cell {r}")
        if t == "fieldsign" and a.get("field") not in field:
            errs.append(f"assumption {a.get('name')}: unknown field {a.get('field')}")
        if not a.get("source"):
            errs.append(f"assumption {a.get('name')}: a documented source is required")
    for tg in spec.get("targets", []):
        for f in tg.get("coefficients", {}):
            if f not in field:
                errs.append(f"target {tg.get('name')}: unknown field {f}")
    for k, cr in enumerate(spec.get("correlations", []) or []):
        cc = cr.get("cells", [])
        if len(cc) != 2 or cc[0] == cc[1] or any(x not in ids for x in cc):
            errs.append(f"correlation {k + 1}: give two different known cell ids")
        rho = cr.get("rho")
        if rho is None or not (-1 < float(rho) < 1):
            errs.append(f"correlation {k + 1}: rho must lie strictly between -1 and 1")
    cse = spec.get("certificate_se")
    if cse is not None:
        if not isinstance(cse, dict):
            errs.append("certificate_se: give a mapping from cell id to a positive standard error")
        else:
            miss = [x for x in ids if x not in cse]
            if miss:
                errs.append(f"certificate_se: missing cells {miss}")
            for x, v in cse.items():
                if x not in ids:
                    errs.append(f"certificate_se: unknown cell {x}")
                elif not (isinstance(v, (int, float)) and v > 0):
                    errs.append(f"certificate_se: {x} must be a positive number")
    if not errs and spec.get("correlations"):
        R = correlation_matrix(spec)
        if np.linalg.eigvalsh(R).min() < -1e-10:
            errs.append("correlations: the implied correlation matrix is not positive semidefinite")
    if errs:
        raise SpecError("; ".join(errs))
    return True


def cell_log_se(c):
    if c.get("scale", "ratio") == "log":
        return float(c["estimate"]), float(c["se"])
    y = float(np.log(c["estimate"]))
    if "se" in c:
        return y, float(c["se"])
    return y, float((np.log(c["upper"]) - np.log(c["lower"])) / (2 * Z975))


def scale_warnings(spec):
    """v1.2.0: the numerical settings are absolute and assume the log-ratio scale. Set limits are found
    by bisection to an absolute tolerance of 1e-4 and the active-set thresholds (signet.core.TOL) are in
    the units of the estimates, so very small or very large standard errors are flagged. The analysis
    still runs; see experiments/audit_v12/tolerance_sensitivity.md."""
    out = []
    for c in spec.get("cells", []):
        s = cell_log_se(c)[1]
        if s < 0.01:
            out.append(f"cell {c.get('id')}: standard error {s:.3g} on the log scale is below 0.01; set limits are computed "
                       f"to an absolute tolerance of 1e-4 ({100 * 1e-4 / s:.1f}% of this standard error)")
        elif s > 100:
            out.append(f"cell {c.get('id')}: standard error {s:.3g} on the log scale is above 100; the numerical settings "
                       "assume estimates on the log-ratio scale")
    return out


def correlation_matrix(spec):
    ids = [c["id"] for c in spec["cells"]]
    R = np.eye(len(ids))
    for cr in spec.get("correlations", []) or []:
        i, j = ids.index(cr["cells"][0]), ids.index(cr["cells"][1])
        R[i, j] = R[j, i] = float(cr["rho"])
    return R


def build(spec):
    validate(spec)
    field = spec.get("field") or ["theta"]
    M = Model(G=len(field), field_names=list(field))
    ys, ses = [], []
    for c in spec["cells"]:
        M.add_row(c["id"], c.get("kernel", [1.0] * len(field)), c.get("processes", []))
        y, s = cell_log_se(c)
        ys.append(y); ses.append(s)
    if spec.get("correlations"):
        S = np.diag(ses)
        M.set_data(ys, S @ correlation_matrix(spec) @ S)
    else:
        M.set_data(ys, ses)
    M.cert_V = None
    if spec.get("certificate_se"):
        Sc = np.diag([float(spec["certificate_se"][c["id"]]) for c in spec["cells"]])
        M.cert_V = Sc @ correlation_matrix(spec) @ Sc
    rungs = {}
    for a in spec.get("assumptions", []):
        t, nm, rg = a["type"], a["name"], int(a.get("rung", 1))
        if t == "sign":
            s = -1 if a.get("direction", "negative") == "negative" else +1
            M.sign(nm, a["process"], s, rows=a.get("rows"), rung=rg)
        elif t == "zero":
            M.zero(nm, a["process"], rows=a.get("rows"), rung=rg)
        elif t == "share":
            M.share(nm, a["process"], rows=a.get("rows"), ratios=a.get("ratios"), rung=rg)
        elif t == "order":
            M.order(nm, a["process"], a["row_lo"], a["row_hi"], rung=rg)
        elif t == "fieldsign":
            s = +1 if a.get("direction", "positive") == "positive" else -1
            M.fieldsign(nm, field.index(a["field"]), s, rung=rg)
        elif t == "monotone":
            M.monotone(nm, direction=+1 if a.get("direction", "increasing") == "increasing" else -1, rung=rg)
        elif t == "bound":
            M.bound(nm, a["process"], a["lo"], a["hi"], rows=a.get("rows"), rung=rg)
        elif t == "linear":
            M.linear(nm, a["coefs"], a.get("sense", "le"), a.get("rhs", 0.0), rung=rg)
        rungs[nm] = rg
    return M, rungs


def _rr(lo, hi):
    f = lambda v: ("0" if v < 0 else "inf") if not np.isfinite(v) else f"{np.exp(v):.3f}"
    return f"{f(lo)} to {f(hi)}"


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (bool, np.bool_)):
        return bool(o)
    if isinstance(o, (float, np.floating)):
        return float(o) if np.isfinite(o) else ("inf" if o > 0 else "-inf")
    if isinstance(o, (int, np.integer)):
        return int(o)
    if isinstance(o, np.ndarray):
        return _jsonable(o.tolist())
    return o


def run(spec):
    M, rungs = build(spec)
    an = spec.get("analysis", {})
    alpha = float(an.get("alpha", 0.05))
    use_rungs = an.get("rungs", sorted(set(rungs.values())))
    tol = an.get("tolerance", [0, 1])
    claims = an.get("claims", ["benefit", "harm"])
    field = spec.get("field") or ["theta"]
    inv = an.get("inversion", "hull")
    from . import __version__
    CV = getattr(M, "cert_V", None)
    out = dict(name=spec.get("name"), input_sha256=canonical_hash(spec), alpha=alpha, signet_version=__version__,
               inversion=inv, certificate_selection="frozen at certificate_se" if CV is not None else "covariance of the estimates",
               warnings=scale_warnings(spec), targets={})
    for tg in spec.get("targets", [{"name": field[0], "coefficients": {field[0]: 1.0}}]):
        cf = [float(tg["coefficients"].get(f, 0.0)) for f in field]
        res = {}
        for r in use_rungs:
            include = [g for g, rg in rungs.items() if rg <= r]
            homogeneous = all(rungs[g] < 4 for g in include)
            rec = dict(include=include, capacity=capacity(M, include, cf))
            kwh = dict(inversion=inv) if not homogeneous else dict(inversion=inv, cert_V=CV)
            st = signet_set(M, include, cf, alpha, **kwh) if homogeneous else conf_set(M, include, cf, alpha, **kwh)
            rec["set"] = None if st["empty"] else dict(lo=st["lo"], hi=st["hi"], ratio_scale=_rr(st["lo"], st["hi"]))
            cc_ = st.get("cc", st)
            rec["inversion_checks"] = dict(outer_set=cc_.get("outer"), bracketed=cc_.get("bracketed"),
                                           hull_extended=cc_.get("hull_extended"), warnings=cc_.get("warnings", []))
            if homogeneous:
                rec["certificate_unique"] = {s: (None if (cr := certificate_range(M, include, cf, s)) is None else cr["unique"])
                                             for s in ("hi", "lo")}
            idf = identified(M, include, cf, M.y)
            rec["identified_at_estimates"] = dict(lo=idf.get("lo"), hi=idf.get("hi"),
                                                  certificate_hi=idf.get("cert_hi"), certificate_lo=idf.get("cert_lo"))
            mt = model_test(M, include, alpha)
            rec["model_test"] = mt
            if mt["reject"]:
                rec["minimal_repairs"] = repairs(M, include, alpha, kmax=2)
            if homogeneous:
                # v1.2.0: no directional p-value is reported when the set at this rung is empty
                rec["claims"] = {c: claim_pvalue(M, include, cf, c, cert_V=CV, alpha=alpha, empty=bool(st["empty"]))
                                 for c in claims}
            for k in tol:
                if k == 0:
                    continue
                t = tolerant_set(M, include, cf, k=k, alpha=alpha,
                                 setfun=None if homogeneous else conf_set, **kwh)
                rec[f"tolerance_k{k}"] = None if t["empty"] else dict(lo=t["lo"], hi=t["hi"], ratio_scale=_rr(t["lo"], t["hi"]))
            if r == max(use_rungs):
                # primary result (manuscript Section 4, step 11): one failed group or one broken cell
                d = default_set(M, include, cf, alpha=alpha, setfun=None if homogeneous else conf_set, **kwh)
                rec["default"] = None if d["empty"] else dict(lo=d["lo"], hi=d["hi"], ratio_scale=_rr(d["lo"], d["hi"]))
                # v1.2.0: what the default can certify whatever the data, and the set without tolerance beside it
                rec["default_capacity"] = default_capacity(M, include, cf)
                rec["n_cells_measuring_target"] = int(sum(1 for row in M.rows if abs(float(np.dot(row["K"], cf))) > 0))
            res[f"rung{r}"] = rec
        out["targets"][tg["name"]] = res
    return _jsonable(out)


def report_md(result):
    L = [f"# SIGNET report: {result.get('name')}", "", f"Input SHA-256: {result['input_sha256']}", ""]
    for w in result.get("warnings", []) or []:
        L += [f"Warning: {w}", ""]
    if result.get("certificate_selection", "").startswith("frozen"):
        L += ["Certificate weights were chosen with the expected standard errors recorded in the specification.", ""]
    for tname, res in result["targets"].items():
        L.append(f"## Target: {tname}")
        L.append("")
        L.append("| Rung | Capacity (upper, lower) | 95% SIGNET set (ratio scale) | k = 1 set | p benefit | p harm | Model test p |")
        L.append("| --- | --- | --- | --- | --- | --- | --- |")
        for rung, rec in res.items():
            cap = rec["capacity"]
            s = rec["set"]["ratio_scale"] if rec["set"] else "empty"
            k1 = rec.get("tolerance_k1")
            k1 = k1["ratio_scale"] if k1 else "n/a"
            cl = rec.get("claims", {})
            _pf = lambda c_: "n/a" if c_ not in cl else ("not reported (empty set)" if cl[c_]["p"] is None else f"{cl[c_]['p']:.4f}")
            pb, ph = _pf("benefit"), _pf("harm")
            mp = rec["model_test"].get("p")
            L.append(f"| {rung} | {'bounded' if cap['above'] else 'unbounded'}, {'bounded' if cap['below'] else 'unbounded'} "
                     f"| {s} | {k1} | {pb} | {ph} | {mp if mp is None else round(mp, 4)} |")
        top = [(rung, rec) for rung, rec in res.items() if "default" in rec]
        if top:
            rung, rec = top[-1]
            dflt = rec["default"]["ratio_scale"] if rec["default"] else "empty"
            L.append("")
            L.append(f"Default at {rung}, tolerant to one failed assumption group or one broken evidence cell: {dflt}")
            L.append(f"Beside it, {rung} without tolerance (conditional on every named assumption): "
                     f"{rec['set']['ratio_scale'] if rec['set'] else 'empty'}")
            dc = rec.get("default_capacity")
            if dc is not None:
                f = lambda b: "can be bounded" if b else "unbounded whatever the data"
                L.append(f"Capacity of the default: upper limit {f(dc['above'])}; lower limit {f(dc['below'])}.")
                if rec.get("n_cells_measuring_target", 2) <= 1:
                    L.append("One cell measures this target, so the default cannot support a directional claim here; "
                             "its silence carries no information about the effect.")
        L.append("")
    return "\n".join(L) + "\n"
