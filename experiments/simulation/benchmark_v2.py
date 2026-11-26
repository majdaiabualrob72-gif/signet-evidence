"""
benchmark_v2.py  (Parts VI and VII of the evaluation)

1. Adds established comparators to the ten original scenarios of benchmark.py, on the
   identical simulated literatures (same seeds):
     FE        inverse-variance fixed effect over all designs that measure the effect
     BA-U      bias-adjusted synthesis, unsigned additive bias priors (Turner-type)
     BA-S      bias-adjusted synthesis, elicited directions and magnitudes (Shapland-type)
     NCO-cal   negative-control calibration of the cells that have a control, then RE pooling
2. Runs five stress scenarios with all methods (original and new):
     extreme outlier, publication bias, heterogeneous populations,
     missing controls and variances, contradictory evidence.
3. Records computing time per method and certificate sparsity.

Usage:
  python3 benchmark_v2.py comparators NREP        # cheap comparators on the 10 original scenarios
  python3 benchmark_v2.py stress NREP [idx ...]   # full method set on the stress scenarios
"""
import json
import os
import sys
import time

import numpy as np
from multiprocessing import Pool
from scipy import stats as sst

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
sys.path.insert(0, HERE)
import benchmark as B                                   # noqa: E402
from decode_core import Ledger, decode, re_meta          # noqa: E402
from signet.core import (signet_set, tolerant_set, model_test, repairs,  # noqa: E402
                         cert_bounds, pool_cell, INF)

Z = float(sst.norm.ppf(0.975))
STRESS = ["extreme outlier", "publication bias", "heterogeneous populations",
          "missing controls and variances", "contradictory evidence"]

# elicited bias priors for the bias-adjusted comparators (generic values, not the simulation truth)
SIGN = dict(H=-1.0, P=-1.0, R=+1.0, C=0.0, A=0.0, S=0.0, Pl=0.0)
SD_SIGNED, SD_UNSIGNED, SD_TRIAL = 0.15, 0.20, 0.05


# ----------------------------------------------------------------------
# stress-scenario literatures
# ----------------------------------------------------------------------
def literature_v2(sc, rng):
    if sc in B.ALL_SCEN:
        return B.literature(sc, rng)
    theta = B.TH
    b = B.draw_biases("reference", rng)
    present = []
    while not any(K > 0 for (_n, K, *_r) in present):
        present = [d for d in B.DESIGNS if rng.random() < d[3]]
        if sc == "missing controls and variances":
            present = [d for d in present if d[0] not in ("NU-NCO", "AC-NCO")]
        if sc == "contradictory evidence":
            names = [d[0] for d in present]
            for must in ("NU-L", "RCT"):
                if must not in names:
                    present.append([d for d in B.DESIGNS if d[0] == must][0])
    studies = []
    for (name, K, procs, _pr, (lo, hi)) in present:
        k = 1 + rng.poisson(2)
        th_d = theta + (rng.normal(0, 0.10) if (sc == "heterogeneous populations" and K > 0) else 0.0)
        mu = K * th_d + sum(b[p] for p in procs if not (p == "P" and name == "NU-L"))
        if sc == "contradictory evidence" and name == "NU-L":
            mu += 0.30          # a design-specific error that contradicts the signed assumptions
        selected = sc == "publication bias" and K > 0 and name != "RCT"
        got, tries = 0, 0
        while got < k and tries < 6 * k:
            tries += 1
            se = rng.uniform(lo, hi)
            yv = mu + rng.normal(0, se)
            if selected and not (yv / se < -1.96 or rng.random() < 0.35):
                continue
            studies.append(dict(design=name, K=K, procs=procs, y=yv, s=se))
            got += 1
        if got == 0:            # keep at least one study per design
            se = rng.uniform(lo, hi)
            studies.append(dict(design=name, K=K, procs=procs, y=mu + rng.normal(0, se), s=se))
    if sc == "extreme outlier":
        j = int(rng.integers(len(studies)))
        studies[j]["y"] += (1 if rng.random() < 0.5 else -1) * 5.0 * studies[j]["s"]
    if sc == "missing controls and variances":
        # v1.1.0: sorted, so that the random draws do not depend on the per-process order of a set of strings
        for name in sorted(set(d["design"] for d in studies)):
            st = [d for d in studies if d["design"] == name]
            med = float(np.median([d["s"] for d in st]))
            for d in st:
                if rng.random() < 0.30:
                    d["s_analysis"] = med          # variance unreported, imputed by the design median
    for d in studies:
        d.setdefault("s_analysis", d["s"])
        d["s_true"] = d["s"]
        d["s"] = d["s_analysis"]                    # analysts only ever see the analysis SE
    return theta, b, present, studies


# ----------------------------------------------------------------------
# established comparators on cell-level data
# ----------------------------------------------------------------------
def cells(present, studies):
    out = []
    for (name, K, procs, *_r) in present:
        st = [d for d in studies if d["design"] == name]
        y, s, info = pool_cell([d["y"] for d in st], [d["s"] for d in st])
        eff = [p for p in procs if not (p == "P" and name == "NU-L")]
        out.append(dict(name=name, K=K, procs=eff, y=y, s=s))
    return out


def bias_prior(cell, signed):
    if cell["name"] == "RCT":
        return 0.0, SD_TRIAL
    if not signed:
        return 0.0, SD_UNSIGNED
    m, v = 0.0, 0.0
    for p in cell["procs"]:
        if SIGN[p] != 0.0:
            m += SIGN[p] * SD_SIGNED
            v += SD_SIGNED ** 2
        else:
            v += SD_UNSIGNED ** 2
    return m, float(np.sqrt(max(v, SD_TRIAL ** 2)))


def bias_adjusted(cl, signed):
    ys, ws = [], []
    for c in cl:
        if c["K"] == 0:
            continue
        m, s = bias_prior(c, signed)
        ys.append(c["y"] - m); ws.append(1 / (c["s"] ** 2 + s ** 2))
    ys, ws = np.array(ys), np.array(ws)
    mu, se = float(np.sum(ws * ys) / np.sum(ws)), float(np.sqrt(1 / np.sum(ws)))
    return mu - Z * se, mu + Z * se


def nco_calibrated(cl):
    by = {c["name"]: c for c in cl}
    ys, ss = [], []
    for c in cl:
        if c["K"] == 0:
            continue
        partner = {"NU-U": "NU-NCO", "NU-L": "NU-NCO", "AC": "AC-NCO"}.get(c["name"])
        if partner in by:
            ys.append(c["y"] - by[partner]["y"]); ss.append(np.sqrt(c["s"] ** 2 + by[partner]["s"] ** 2))
        else:
            ys.append(c["y"]); ss.append(c["s"])
    r = re_meta(ys, ss)
    return r["lo"], r["hi"]


def fixed_effect(studies):
    main = [d for d in studies if d["K"] > 0]
    w = np.array([1 / d["s"] ** 2 for d in main]); y = np.array([d["y"] for d in main])
    mu, se = float(np.sum(w * y) / np.sum(w)), float(np.sqrt(1 / np.sum(w)))
    return mu - Z * se, mu + Z * se


def comparators(theta, present, studies):
    cl = cells(present, studies)
    out = {}
    out["FE (all designs)"] = B.summarise_set(*fixed_effect(studies), theta)
    out["BA-U (unsigned bias priors)"] = B.summarise_set(*bias_adjusted(cl, False), theta)
    out["BA-S (elicited directions and magnitudes)"] = B.summarise_set(*bias_adjusted(cl, True), theta)
    out["NCO-cal (calibrated, RE)"] = B.summarise_set(*nco_calibrated(cl), theta)
    return out


def one_comparators(args):
    sc, seed = args
    rng = np.random.default_rng(seed)
    theta, b, present, studies = literature_v2(sc, rng)
    out = dict(theta=theta)
    t0 = time.perf_counter()
    out.update(comparators(theta, present, studies))
    out["t_comparators"] = time.perf_counter() - t0
    return out


# ----------------------------------------------------------------------
# full method set (stress scenarios)
# ----------------------------------------------------------------------
def one_full(args):
    sc, seed = args
    rng = np.random.default_rng(seed)
    theta, b, present, studies = literature_v2(sc, rng)
    out = dict(theta=theta)
    t = {}
    t0 = time.perf_counter()
    main = [d for d in studies if d["K"] > 0]
    ma = re_meta([d["y"] for d in main], [d["s"] for d in main])
    out["RE (all designs)"] = B.summarise_set(ma["lo"], ma["hi"], theta)
    tr = [d for d in studies if d["design"] == "RCT"]
    if tr:
        mt = re_meta([d["y"] for d in tr], [d["s"] for d in tr])
        out["RE (trials only)"] = B.summarise_set(mt["lo"], mt["hi"], theta)
    else:
        out["RE (trials only)"] = dict(cover=False, width=None, unbounded=True, empty=False,
                                       benefit=False, harm=False, absent=True)
    out.update(comparators(theta, present, studies))
    t["established"] = time.perf_counter() - t0
    # sparse Bayesian decoder with magnitude priors (DECODE; the author's earlier design)
    t0 = time.perf_counter()
    Zm = np.array([[d["K"]] for d in studies])
    S = np.zeros((len(studies), len(B.PROCS)))
    for i, d in enumerate(studies):
        for p in d["procs"]:
            if p == "P" and d["design"] == "NU-L":
                continue
            S[i, B.PROCS.index(p)] = 1.0
    used = S.sum(axis=0) > 0
    L = Ledger(np.array([d["y"] for d in studies]), np.array([d["s"] for d in studies]), Zm, S[:, used],
               np.array(B.DEC_PI)[used], np.array(B.DEC_SIG)[used])
    lo, hi = decode(L, tau_grid=B.TAU).theta(0).interval()
    out["DECODE"] = B.summarise_set(lo, hi, theta)
    t["DECODE"] = time.perf_counter() - t0
    # SIGNET
    m, R = B.signet_model(present, studies)
    kw = dict(exact_rank=False, ngrid=24, tol=5e-4, outer_scan=0)
    for rn in ("R1", "R2", "R3"):
        t0 = time.perf_counter()
        r_ = signet_set(m, R[rn], [1.0], **kw)
        t[f"SIGNET {rn}"] = time.perf_counter() - t0
        out[f"SIGNET {rn}"] = B.summarise_set(r_["lo"], r_["hi"], theta, r_["empty"])
    t0 = time.perf_counter()
    tk = tolerant_set(m, R["R3"], [1.0], k=1, **kw)
    t["SIGNET R3 k=1"] = time.perf_counter() - t0
    out["SIGNET R3 k=1"] = B.summarise_set(tk["lo"], tk["hi"], theta, tk["empty"])
    mt3 = model_test(m, R["R3"], exact_rank=False)
    out["R3 rejected"] = bool(mt3["reject"])
    if not mt3["reject"]:
        work = R["R3"]; out["protocol rung"] = "R3"
    else:
        rp = repairs(m, R["R3"], kmax=1, droppable=[g for g in R["R3"] if g not in R["R1"]])
        if len(rp) == 1 and len(rp[0]) == 1:
            work = [g for g in R["R3"] if g not in rp[0]]; out["protocol rung"] = "R3 repaired"
        else:
            work = R["R2"]; out["protocol rung"] = "R2"
    r_ = signet_set(m, work, [1.0], **kw) if work != R["R3"] else None
    out["SIGNET protocol"] = out["SIGNET R3"] if r_ is None else B.summarise_set(r_["lo"], r_["hi"], theta, r_["empty"])
    # certificate sparsity at rung 3
    cb = cert_bounds(m, R["R3"], [1.0], a=0.0125)
    sp = []
    for side in ("hi", "lo"):
        c = cb["cert_" + side]
        if c is not None:
            sp.append(dict(side=side, cells=len(c["weights"]), assumptions=len(c["assumptions"]), n_cells=m.n))
    out["certificates R3"] = sp
    out["times"] = t
    return out


METHODS_FULL = ["RE (all designs)", "RE (trials only)", "FE (all designs)", "BA-U (unsigned bias priors)",
                "BA-S (elicited directions and magnitudes)", "NCO-cal (calibrated, RE)", "DECODE",
                "SIGNET R1", "SIGNET R2", "SIGNET R3", "SIGNET protocol", "SIGNET R3 k=1"]
METHODS_CMP = ["FE (all designs)", "BA-U (unsigned bias priors)", "BA-S (elicited directions and magnitudes)",
               "NCO-cal (calibrated, RE)"]


def summarise(reps, methods):
    s = {}
    for mth in methods:
        rr = [r[mth] for r in reps if not r[mth].get("absent")]
        w = [x["width"] for x in rr if x["width"] is not None]
        s[mth] = dict(n=len(rr),
                      coverage=float(np.mean([x["cover"] for x in rr])) if rr else None,
                      mc_se=float(np.sqrt(np.mean([x["cover"] for x in rr]) * (1 - np.mean([x["cover"] for x in rr])) / len(rr))) if rr else None,
                      median_width=float(np.median(w)) if w else None,
                      unbounded=float(np.mean([x["unbounded"] for x in rr])) if rr else None,
                      empty=float(np.mean([x["empty"] for x in rr])) if rr else None,
                      benefit=float(np.mean([x["benefit"] for x in rr])) if rr else None,
                      harm=float(np.mean([x["harm"] for x in rr])) if rr else None)
    if "R3 rejected" in reps[0]:
        s["R3 rejection rate"] = float(np.mean([r["R3 rejected"] for r in reps]))
    if "times" in reps[0]:
        s["median seconds"] = {k: float(np.median([r["times"][k] for r in reps])) for k in reps[0]["times"]}
        cs = [c for r in reps for c in r["certificates R3"]]
        if cs:
            s["certificate R3 sparsity"] = dict(mean_cells_used=float(np.mean([c["cells"] for c in cs])),
                                                mean_cells_available=float(np.mean([c["n_cells"] for c in cs])),
                                                mean_assumption_groups=float(np.mean([c["assumptions"] for c in cs])),
                                                n_certificates=len(cs))
    return s


# ----------------------------------------------------------------------
# magnitude-misspecification sweep (where does each method win?)
# ----------------------------------------------------------------------
SWEEP_MULTS = [0.5, 1.0, 2.0, 3.0, 4.0]
SWEEP_FAMILIES = ["mixed designs, true benefit", "one-signed literature, null effect"]


def literature_sweep(family, mult, rng):
    theta = B.TH if family.startswith("mixed") else 0.0
    b = dict(H=-abs(rng.normal(0.12, 0.04)) * mult, P=-abs(rng.normal(0.15, 0.05)) * mult,
             R=abs(rng.normal(0.10, 0.04)) * mult, C=rng.normal(0, 0.04),
             A=0.0 if rng.random() < 0.85 else rng.normal(0, 0.12 * mult),
             S=0.0 if rng.random() < 0.70 else rng.normal(0, 0.20 * mult),
             Pl=0.0 if rng.random() < 0.80 else rng.normal(0, 0.15 * mult))
    present = []
    while not any(K > 0 for (_n, K, *_r) in present):
        present = [d for d in B.DESIGNS if rng.random() < d[3]]
        if family.startswith("one-signed"):
            present = [d for d in present if d[0] in ("NU-U", "NU-L")]
    studies = []
    for (name, K, procs, _pr, (lo, hi)) in present:
        k = 1 + rng.poisson(2)
        mu = K * theta + sum(b[p] for p in procs if not (p == "P" and name == "NU-L"))
        for _ in range(k):
            se = rng.uniform(lo, hi)
            studies.append(dict(design=name, K=K, procs=procs, y=mu + rng.normal(0, se), s=se))
    return theta, b, present, studies


def one_sweep(args):
    family, mult, seed = args
    rng = np.random.default_rng(seed)
    theta, b, present, studies = literature_sweep(family, mult, rng)
    out = dict(theta=theta)
    main = [d for d in studies if d["K"] > 0]
    ma = re_meta([d["y"] for d in main], [d["s"] for d in main])
    out["RE (all designs)"] = B.summarise_set(ma["lo"], ma["hi"], theta)
    out.update(comparators(theta, present, studies))
    m, R = B.signet_model(present, studies)
    kw = dict(exact_rank=False, ngrid=24, tol=5e-4, outer_scan=0)
    r_ = signet_set(m, R["R3"], [1.0], **kw)
    out["SIGNET R3"] = B.summarise_set(r_["lo"], r_["hi"], theta, r_["empty"])
    tk = tolerant_set(m, R["R3"], [1.0], k=1, **kw)
    out["SIGNET R3 k=1"] = B.summarise_set(tk["lo"], tk["hi"], theta, tk["empty"])
    out.update(signet_bounded(m, R, theta))
    return out


METHODS_SWEEP = ["RE (all designs)", "FE (all designs)", "BA-U (unsigned bias priors)",
                 "BA-S (elicited directions and magnitudes)", "NCO-cal (calibrated, RE)", "SIGNET R3", "SIGNET R3 k=1",
                 "SIGNET-B (rung 3 + magnitude bounds)", "SIGNET-B k=1"]


# ----------------------------------------------------------------------
# redesign 2: an explicit, generous magnitude rung (bounds at three times the elicited typical size)
# ----------------------------------------------------------------------
BOUND_SIGNED, BOUND_UNSIGNED = 3 * SD_SIGNED, 3 * SD_UNSIGNED


def add_bounds(m, R):
    names = []
    for p in B.PROCS:
        if len(m.rows_with(p)):
            u = BOUND_SIGNED if SIGN[p] != 0.0 else BOUND_UNSIGNED
            nm = f"bound{p}"
            m.bound(nm, p, -u, u, rung=4)
            names.append(nm)
    return R["R3"] + names


def signet_bounded(m, R, theta):
    from signet.core import conf_set
    R4 = add_bounds(m, R)
    kw = dict(exact_rank=False, ngrid=24, tol=5e-4, outer_scan=0)
    cs = conf_set(m, R4, [1.0], 0.05, **kw)
    out = {"SIGNET-B (rung 3 + magnitude bounds)": B.summarise_set(cs.get("lo"), cs.get("hi"), theta, cs["empty"])}
    tk = tolerant_set(m, R4, [1.0], k=1, droppable=R["R3"], setfun=conf_set, **kw)
    out["SIGNET-B k=1"] = B.summarise_set(tk["lo"], tk["hi"], theta, tk["empty"])
    return out


def one_bounded(args):
    sc, seed = args
    rng = np.random.default_rng(seed)
    theta, b, present, studies = literature_v2(sc, rng)
    m, R = B.signet_model(present, studies)
    t0 = time.perf_counter()
    out = dict(theta=theta)
    out.update(signet_bounded(m, R, theta))
    out["t_bounded"] = time.perf_counter() - t0
    return out


METHODS_B = ["SIGNET-B (rung 3 + magnitude bounds)", "SIGNET-B k=1"]


# ----------------------------------------------------------------------
# redesign 3: a named allowance for between-population effect heterogeneity
# (cells that measure the effect get variance se^2 + tau_h^2; controls do not)
# ----------------------------------------------------------------------
TAU_H = 0.10


def one_hetallow(args):
    sc, seed = args
    rng = np.random.default_rng(seed)
    theta, b, present, studies = literature_v2(sc, rng)
    m, R = B.signet_model(present, studies)
    for i, r in enumerate(m.rows):
        if r["K"].sum() > 0:
            m.V[i, i] += TAU_H ** 2
    kw = dict(exact_rank=False, ngrid=24, tol=5e-4, outer_scan=0)
    out = dict(theta=theta)
    r_ = signet_set(m, R["R3"], [1.0], **kw)
    out["SIGNET R3 + heterogeneity allowance"] = B.summarise_set(r_["lo"], r_["hi"], theta, r_["empty"])
    tk = tolerant_set(m, R["R3"], [1.0], k=1, **kw)
    out["SIGNET R3 k=1 + heterogeneity allowance"] = B.summarise_set(tk["lo"], tk["hi"], theta, tk["empty"])
    return out


METHODS_H = ["SIGNET R3 + heterogeneity allowance", "SIGNET R3 k=1 + heterogeneity allowance"]


# ----------------------------------------------------------------------
# redesign 4: tolerance to one broken evidence cell as well as one failed group
# union of {sets omitting any one assumption group} and {sets omitting any one evidence cell}
# ----------------------------------------------------------------------
def cell_robust_set(present, studies, kw):
    m, R = B.signet_model(present, studies)
    t = tolerant_set(m, R["R3"], [1.0], k=1, **kw)
    lo, hi = (INF, -INF) if t["empty"] else (t["lo"], t["hi"])
    for d in present:
        rest = [e for e in present if e[0] != d[0]]
        if not any(e[1] > 0 for e in rest):
            continue
        st = [x for x in studies if x["design"] != d[0]]
        m2, R2 = B.signet_model(rest, st)
        r_ = signet_set(m2, R2["R3"], [1.0], **kw)
        if not r_["empty"]:
            lo, hi = min(lo, r_["lo"]), max(hi, r_["hi"])
    if lo == INF:
        return dict(lo=None, hi=None, empty=True)
    return dict(lo=lo, hi=hi, empty=False)


def one_cellrobust(args):
    sc, seed = args
    rng = np.random.default_rng(seed)
    theta, b, present, studies = literature_v2(sc, rng)
    kw = dict(exact_rank=False, ngrid=24, tol=5e-4, outer_scan=0)
    t0 = time.perf_counter()
    r_ = cell_robust_set(present, studies, kw)
    out = dict(theta=theta)
    out["SIGNET R3 k=1 (group or cell)"] = B.summarise_set(r_["lo"], r_["hi"], theta, r_["empty"])
    out["t_cellrobust"] = time.perf_counter() - t0
    return out


METHODS_CR = ["SIGNET R3 k=1 (group or cell)"]


# ----------------------------------------------------------------------
# iteration 3: cell-level tolerance combined with the heterogeneity allowance
# ----------------------------------------------------------------------
def _with_allowance(m):
    for i, r in enumerate(m.rows):
        if r["K"].sum() > 0:
            m.V[i, i] += TAU_H ** 2
    return m


def cell_robust_het_set(present, studies, kw):
    m, R = B.signet_model(present, studies)
    _with_allowance(m)
    t = tolerant_set(m, R["R3"], [1.0], k=1, **kw)
    lo, hi = (INF, -INF) if t["empty"] else (t["lo"], t["hi"])
    for d in present:
        rest = [e for e in present if e[0] != d[0]]
        if not any(e[1] > 0 for e in rest):
            continue
        st = [x for x in studies if x["design"] != d[0]]
        m2, R2 = B.signet_model(rest, st)
        _with_allowance(m2)
        r_ = signet_set(m2, R2["R3"], [1.0], **kw)
        if not r_["empty"]:
            lo, hi = min(lo, r_["lo"]), max(hi, r_["hi"])
    if lo == INF:
        return dict(lo=None, hi=None, empty=True)
    return dict(lo=lo, hi=hi, empty=False)


def one_cellhet(args):
    sc, seed = args
    rng = np.random.default_rng(seed)
    theta, b, present, studies = literature_v2(sc, rng)
    kw = dict(exact_rank=False, ngrid=24, tol=5e-4, outer_scan=0)
    t0 = time.perf_counter()
    r_ = cell_robust_het_set(present, studies, kw)
    out = dict(theta=theta)
    out["SIGNET R3 k=1 (group or cell) + heterogeneity allowance"] = B.summarise_set(r_["lo"], r_["hi"], theta, r_["empty"])
    out["t_cellhet"] = time.perf_counter() - t0
    return out


METHODS_CH = ["SIGNET R3 k=1 (group or cell) + heterogeneity allowance"]


if __name__ == "__main__":
    mode = sys.argv[1]
    nrep = int(sys.argv[2])
    os.chdir(HERE)
    if mode == "comparators":
        res = {}
        with Pool(2) as pool:
            for i, sc in enumerate(B.ALL_SCEN):
                reps = pool.map(one_comparators, [(sc, 7000 * (i + 1) + r) for r in range(nrep)], chunksize=16)
                res[sc] = summarise(reps, METHODS_CMP)
                print(sc, {k: (round(v["coverage"], 3), v["median_width"] and round(v["median_width"], 2),
                               round(v["benefit"], 2)) for k, v in res[sc].items() if isinstance(v, dict) and "coverage" in v},
                      flush=True)
        json.dump(dict(nrep=nrep, results=res), open("benchmark_v2_comparators.json", "w"), indent=1)
    elif mode == "bounded":
        res = {}
        jobs = [(sc, [7000 * (i + 1) + r for r in range(nrep)]) for i, sc in enumerate(B.ALL_SCEN)]
        jobs += [(sc, [91000 + 1000 * i + r for r in range(nrep)]) for i, sc in enumerate(STRESS)]
        with Pool(2) as pool:
            for sc, seeds in jobs:
                reps = pool.map(one_bounded, [(sc, sd) for sd in seeds], chunksize=4)
                res[sc] = summarise(reps, METHODS_B)
                res[sc]["median seconds"] = float(np.median([r["t_bounded"] for r in reps]))
                json.dump(dict(nrep=nrep, results=res), open("benchmark_v2_bounded.json", "w"), indent=1)
                print(sc, {k: (round(v["coverage"], 3), v["median_width"] and round(v["median_width"], 2), round(v["benefit"], 2))
                           for k, v in res[sc].items() if isinstance(v, dict) and "coverage" in v}, flush=True)
    elif mode == "hetallow":
        res = {}
        out_fn = "benchmark_v2_hetallow.json"
        if len(sys.argv) > 3:                      # sensitivity to the allowance: python3 benchmark_v2.py hetallow 400 0.05
            TAU_H = float(sys.argv[3])
            out_fn = f"benchmark_v2_hetallow_tau{sys.argv[3]}.json"
        jobs = [("reference", [7000 * 1 + r for r in range(nrep)]),
                ("heterogeneous populations", [91000 + 1000 * 2 + r for r in range(nrep)])]
        with Pool(2) as pool:
            for sc, seeds in jobs:
                reps = pool.map(one_hetallow, [(sc, sd) for sd in seeds], chunksize=4)
                res[sc] = summarise(reps, METHODS_H)
                json.dump(dict(nrep=nrep, tau_h=TAU_H, results=res), open(out_fn, "w"), indent=1)
                print(sc, {k: (round(v["coverage"], 3), v["median_width"] and round(v["median_width"], 2), round(v["benefit"], 2))
                           for k, v in res[sc].items() if isinstance(v, dict) and "coverage" in v}, flush=True)
    elif mode == "cellrobust":
        res = {}
        # stress scenarios first (contradictory evidence leads); seeds are fixed per scenario, so order does not change results
        jobs = [(sc, [91000 + 1000 * i + r for r in range(nrep)]) for i, sc in enumerate(STRESS)][::-1]
        jobs += [(sc, [7000 * (i + 1) + r for r in range(nrep)]) for i, sc in enumerate(B.ALL_SCEN)]
        with Pool(2) as pool:
            for sc, seeds in jobs:
                reps = pool.map(one_cellrobust, [(sc, sd) for sd in seeds], chunksize=4)
                res[sc] = summarise(reps, METHODS_CR)
                res[sc]["median seconds"] = float(np.median([r["t_cellrobust"] for r in reps]))
                json.dump(dict(nrep=nrep, results=res), open("benchmark_v2_cellrobust.json", "w"), indent=1)
                print(sc, {k: (round(v["coverage"], 3), v["median_width"] and round(v["median_width"], 2), round(v["benefit"], 2))
                           for k, v in res[sc].items() if isinstance(v, dict) and "coverage" in v}, flush=True)
    elif mode == "cellhet":
        res = {}
        jobs = [("contradictory evidence", [91000 + 1000 * 4 + r for r in range(nrep)]),
                ("heterogeneous populations", [91000 + 1000 * 2 + r for r in range(nrep)]),
                ("reference", [7000 * 1 + r for r in range(nrep)])]
        with Pool(2) as pool:
            for sc, seeds in jobs:
                reps = pool.map(one_cellhet, [(sc, sd) for sd in seeds], chunksize=4)
                res[sc] = summarise(reps, METHODS_CH)
                res[sc]["median seconds"] = float(np.median([r["t_cellhet"] for r in reps]))
                json.dump(dict(nrep=nrep, tau_h=TAU_H, results=res), open("benchmark_v2_cellhet.json", "w"), indent=1)
                print(sc, {k: (round(v["coverage"], 3), v["median_width"] and round(v["median_width"], 2), round(v["benefit"], 2))
                           for k, v in res[sc].items() if isinstance(v, dict) and "coverage" in v}, flush=True)
    elif mode == "sweep":
        res = {}
        with Pool(2) as pool:
            for fi, fam in enumerate(SWEEP_FAMILIES):
                for mi, mult in enumerate(SWEEP_MULTS):
                    reps = pool.map(one_sweep, [(fam, mult, 500000 + 100000 * fi + 10000 * mi + r) for r in range(nrep)], chunksize=4)
                    key = f"{fam} | x{mult}"
                    res[key] = summarise(reps, METHODS_SWEEP)
                    json.dump(dict(nrep=nrep, results=res), open("benchmark_v2_sweep.json", "w"), indent=1)
                    print(key, {k: (round(v["coverage"], 3), round(v["benefit"], 2)) for k, v in res[key].items()
                                if isinstance(v, dict) and "coverage" in v}, flush=True)
    elif mode == "stress":
        idx = [int(a) for a in sys.argv[3:]] or list(range(len(STRESS)))
        try:
            res = json.load(open("benchmark_v2_stress.json"))["results"]
        except FileNotFoundError:
            res = {}
        t_all = time.time()
        with Pool(2) as pool:
            for i in idx:
                sc = STRESS[i]
                fn = f"benchmark_v2_stress_reps_{i}.json"
                try:
                    reps = json.load(open(fn)); assert len(reps) == nrep
                except (FileNotFoundError, AssertionError):
                    reps = pool.map(one_full, [(sc, 91000 + 1000 * i + r) for r in range(nrep)], chunksize=4)
                    json.dump(reps, open(fn, "w"))
                res[sc] = summarise(reps, METHODS_FULL)
                json.dump(dict(nrep=nrep, results=res), open("benchmark_v2_stress.json", "w"), indent=1)
                print(sc, f"{time.time() - t_all:.0f}s", flush=True)
                for mth in METHODS_FULL:
                    v = res[sc][mth]
                    print(f"   {mth:<44} cov {B.fmt(v['coverage'])} width {B.fmt(v['median_width'])} "
                          f"unb {B.fmt(v['unbounded'], 2)} benefit {B.fmt(v['benefit'], 2)} harm {B.fmt(v['harm'], 2)}",
                          flush=True)
