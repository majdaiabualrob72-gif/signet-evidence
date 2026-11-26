"""
sim_audit.py  (signet 1.2.0; analyses added in response to the forensic audit of 22 November 2026)

Every job re-uses the archived generators, seeds and SIGNET settings of experiments/simulation, so each
record pairs, literature by literature, with the archived results.

Jobs (one JSON file per job and scenario in experiments/audit_v12/out):
  ablation   D10, D9, M3  rung 3 and the tolerant sets with the conditional test alone at the full level
                          (alpha = 0.05, no certificate intersection, no Bonferroni split), beside the
                          reported SIGNET sets; limits are stored, so bounded-set coverage can be computed.
                          400 literatures per scenario for rung 3 and group tolerance; the first 200 also
                          get the cell-tolerant set and hence the reporting default (package functions).
  groups     D3           group tolerance under two other readings of "one assumption group":
                          one group per constraint row, and one group per bias process.
  heldout    M3           the frozen default on literatures generated after it was fixed: new seeds for the
                          15 scenarios, and three new stress scenarios specified here before any was run.
  comparators M1, D2      random-effects intervals (Wald, Hartung-Knapp, modified Hartung-Knapp) and BA-S
                          with independent or shared bias terms, on the archived literatures.

Usage:
  python3 sim_audit.py worker        # claims and runs jobs until none is left (start one per core)
  python3 sim_audit.py comparators   # cheap, runs in the foreground
  python3 sim_audit.py status
"""
import json
import os
import sys
import time

import numpy as np
from scipy import stats as sst

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SIM = os.path.join(ROOT, "experiments", "simulation")
OUT = os.path.join(HERE, "out")
sys.path.insert(0, ROOT)
sys.path.insert(0, SIM)
import benchmark as B                                    # noqa: E402
import benchmark_v2 as B2                                # noqa: E402
from decode_core import re_meta                          # noqa: E402
from signet.core import (signet_set, conf_set, tolerant_set, cell_tolerant_set, _union,   # noqa: E402
                         pool_cell)

KW = dict(exact_rank=False, ngrid=24, tol=5e-4, outer_scan=0)      # the archived simulation settings
Z = float(sst.norm.ppf(0.975))
SCEN15 = list(B.ALL_SCEN) + list(B2.STRESS)


def seed_of(sc, r):
    """Archived seeds: base scenario i uses 7000 (i + 1) + r, stress scenario j uses 91000 + 1000 j + r."""
    if sc in B.ALL_SCEN:
        return 7000 * (B.ALL_SCEN.index(sc) + 1) + r
    return 91000 + 1000 * B2.STRESS.index(sc) + r


# ----------------------------------------------------------------------
# held-out literatures (M3): specified on 22 November 2026, before any of them was generated
# ----------------------------------------------------------------------
HELD_NEW = {
    "held-out: registry against its direction": dict(must=["REG", "RCT"], shift={"REG": -0.30}),
    "held-out: new-user cohort against its direction, larger shift": dict(must=["NU-U", "RCT"], shift={"NU-U": +0.50}),
    "held-out: two designs against their directions": dict(must=["NU-L", "REG", "RCT"], shift={"NU-L": +0.30, "REG": -0.30}),
}
HELD_ALL = SCEN15 + list(HELD_NEW)


def held_seed(sc, r):
    """New seeds, disjoint from every archived seed (archived seeds are below 1,000,000)."""
    return 5_000_000 + 10_000 * HELD_ALL.index(sc) + r


def literature_any(sc, rng):
    if sc not in HELD_NEW:
        return B2.literature_v2(sc, rng)
    spec = HELD_NEW[sc]
    theta = B.TH
    b = B.draw_biases("reference", rng)
    present = []
    while not any(K > 0 for (_n, K, *_r) in present):
        present = [d for d in B.DESIGNS if rng.random() < d[3]]
        names = [d[0] for d in present]
        for must in spec["must"]:
            if must not in names:
                present.append([d for d in B.DESIGNS if d[0] == must][0])
    studies = []
    for (name, K, procs, _pr, (lo, hi)) in present:
        k = 1 + rng.poisson(2)
        mu = K * theta + sum(b[p] for p in procs if not (p == "P" and name == "NU-L")) + spec["shift"].get(name, 0.0)
        for _ in range(k):
            se = rng.uniform(lo, hi)
            studies.append(dict(design=name, K=K, procs=procs, y=mu + rng.normal(0, se), s=se))
    return theta, b, present, studies


# ----------------------------------------------------------------------
# records
# ----------------------------------------------------------------------
def rec(s, theta):
    if s.get("empty") or s.get("lo") is None:
        return dict(lo=None, hi=None, empty=True, cover=False, benefit=False, harm=False, bounded=False, upper=False, lower=False)
    lo, hi = float(s["lo"]), float(s["hi"])
    return dict(lo=lo if np.isfinite(lo) else None, hi=hi if np.isfinite(hi) else None, empty=False,
                cover=bool(lo <= theta <= hi), benefit=bool(hi < 0), harm=bool(lo > 0),
                bounded=bool(np.isfinite(lo) and np.isfinite(hi)), upper=bool(np.isfinite(hi)), lower=bool(np.isfinite(lo)))


def rec_interval(lo, hi, theta):
    return rec(dict(lo=lo, hi=hi, empty=False), theta)


# ----------------------------------------------------------------------
# comparators (M1, D2)
# ----------------------------------------------------------------------
def re_variants(y, s):
    """DerSimonian-Laird estimate with three intervals: Wald (normal quantile), Hartung-Knapp (variance
    factor q, t quantile) and modified Hartung-Knapp (factor max(1, q), t quantile; the interval used by
    the archived comparator, decode_core.re_meta). One study: all three are the Wald interval."""
    y = np.asarray(y, float); s = np.asarray(s, float); k = len(y)
    w = 1 / s ** 2
    mu_fe = np.sum(w * y) / np.sum(w)
    Q = np.sum(w * (y - mu_fe) ** 2)
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - (k - 1)) / C) if k > 1 else 0.0
    ws = 1 / (s ** 2 + tau2)
    mu = float(np.sum(ws * y) / np.sum(ws))
    se = float(np.sqrt(1 / np.sum(ws)))
    out = dict(est=mu, tau2=float(tau2), k=int(k), se_dl=se)
    out["wald"] = (mu - Z * se, mu + Z * se, float(2 * sst.norm.sf(abs(mu / se))))
    if k > 1:
        q = float(np.sum(ws * (y - mu) ** 2) / (k - 1))
        tq = float(sst.t.ppf(0.975, k - 1))
        for nm, f in (("hk", q), ("mhk", max(1.0, q))):
            se_ = se * np.sqrt(f)
            p = 1.0 if se_ == 0 else float(2 * sst.t.sf(abs(mu / se_), k - 1))
            out[nm] = (mu - tq * se_, mu + tq * se_, p)
        out["q"] = q
    else:
        out["hk"] = out["mhk"] = out["wald"]
        out["q"] = None
    return out


def ba_s(cl, shared):
    """BA-S on the cells that measure the effect. Prior for each bias process: signed, N(direction x 0.15,
    0.15^2); unsigned, N(0, 0.20^2); trial cells, N(0, 0.05^2).
    shared=False (archived specification): each cell has its own, independent bias.
    shared=True: a process carried by several cells is one term common to them, so their prior biases are
    correlated (generalised least squares with the implied prior covariance)."""
    cells = [c for c in cl if c["K"] != 0]
    n = len(cells)
    m = np.zeros(n); Om = np.zeros((n, n))
    for i, c in enumerate(cells):
        mi, si = B2.bias_prior(c, True)
        m[i] = mi; Om[i, i] = si ** 2
    if shared:
        for i, ci in enumerate(cells):
            for j, cj in enumerate(cells):
                if i < j and ci["name"] != "RCT" and cj["name"] != "RCT":
                    common = set(ci["procs"]) & set(cj["procs"])
                    cov = sum((B2.SD_SIGNED if B2.SIGN[p] != 0.0 else B2.SD_UNSIGNED) ** 2 for p in common)
                    Om[i, j] = Om[j, i] = cov
    S = np.diag([c["s"] ** 2 for c in cells]) + Om
    Si = np.linalg.inv(S)
    one = np.ones(n)
    var = 1.0 / float(one @ Si @ one)
    mu = var * float(one @ Si @ (np.array([c["y"] for c in cells]) - m))
    se = float(np.sqrt(var))
    return mu - Z * se, mu + Z * se


def comparators_one(sc, seed, gen=None):
    rng = np.random.default_rng(seed)
    theta, b, present, studies = (gen or literature_any)(sc, rng) if not isinstance(sc, tuple) else gen(sc, rng)
    main = [d for d in studies if d["K"] > 0]
    rv = re_variants([d["y"] for d in main], [d["s"] for d in main])
    cl = B2.cells(present, studies)
    out = dict(seed=int(seed), k_studies=rv["k"], tau2=rv["tau2"])
    for nm in ("wald", "hk", "mhk"):
        out["RE " + nm] = rec_interval(rv[nm][0], rv[nm][1], theta)
    out["BA-S independent"] = rec_interval(*ba_s(cl, False), theta)
    out["BA-S shared"] = rec_interval(*ba_s(cl, True), theta)
    out["n_effect_cells"] = int(sum(1 for c in cl if c["K"] != 0))
    out["designs"] = [c["name"] for c in cl]
    return out


# ----------------------------------------------------------------------
# group readings (D3)
# ----------------------------------------------------------------------
PROC_OF = {"H-": "H", "shH": "H", "P-": "P", "lagP0": "P", "R+": "R", "shC": "C", "vA": "A", "vS": "S", "vPl": "Pl"}


def regroup(m, names, reading):
    """Replace the assumption groups of model m (the archived, statement-level groups listed in names)
    by another reading of 'one group', keeping every constraint row.
      'row'      each constraint row is its own group
      'process'  every statement about one bias process forms one group (signs, links and validity)
    Returns the new list of group names."""
    old = [g for g in m.groups if g["name"] in names]
    assert len(old) == len(names)
    new = []
    if reading == "row":
        for g in old:
            for k, con in enumerate(g["cons"]):
                new.append(dict(name=f"{g['name']}#{k}", rung=g["rung"], kind=g["kind"], cons=[con]))
    elif reading == "process":
        by = {}
        for g in old:
            by.setdefault(PROC_OF[g["name"]], []).extend(g["cons"])
        for p, cons in by.items():
            new.append(dict(name=f"process {p}", rung=3, kind="merged", cons=cons))
    else:
        raise ValueError(reading)
    m.groups = new
    return [g["name"] for g in new]


# ----------------------------------------------------------------------
# jobs
# ----------------------------------------------------------------------
def job_ablation(sc, r):
    seed = seed_of(sc, r)
    rng = np.random.default_rng(seed)
    theta, b, present, studies = B2.literature_v2(sc, rng)
    m, R = B.signet_model(present, studies)
    out = dict(seed=int(seed), n_cells=int(m.n), n_effect_cells=int(sum(1 for x in m.rows if x["K"].sum() > 0)),
               n_groups=len(R["R3"]))
    out["R3"] = rec(signet_set(m, R["R3"], [1.0], **KW), theta)
    out["R3 test only"] = rec(conf_set(m, R["R3"], [1.0], 0.05, **KW), theta)
    k1t = tolerant_set(m, R["R3"], [1.0], k=1, setfun=conf_set, **KW)
    out["K1 test only"] = rec(k1t, theta)
    if r < 200:
        k1 = tolerant_set(m, R["R3"], [1.0], k=1, **KW)
        out["K1"] = rec(k1, theta)
        ct = cell_tolerant_set(m, R["R3"], [1.0], alpha=0.05, **KW)
        ctt = cell_tolerant_set(m, R["R3"], [1.0], alpha=0.05, setfun=conf_set, **KW)
        out["cell tolerant"] = rec(ct, theta)
        out["cell tolerant test only"] = rec(ctt, theta)
        out["default"] = rec(_union([k1, ct]), theta)
        out["default test only"] = rec(_union([k1t, ctt]), theta)
    return out


def job_groups(sc, r):
    seed = seed_of(sc, r)
    out = dict(seed=int(seed))
    for reading in ("row", "process"):
        rng = np.random.default_rng(seed)
        theta, b, present, studies = B2.literature_v2(sc, rng)
        m, R = B.signet_model(present, studies)
        names = regroup(m, R["R3"], reading)
        out[f"K1 {reading}"] = rec(tolerant_set(m, names, [1.0], k=1, **KW), theta)
        out[f"n_groups {reading}"] = len(names)
    return out


def job_heldout(sc, r):
    seed = held_seed(sc, r)
    rng = np.random.default_rng(seed)
    theta, b, present, studies = literature_any(sc, rng)
    m, R = B.signet_model(present, studies)
    out = comparators_one(sc, seed)
    out["n_cells"] = int(m.n)
    out["R3"] = rec(signet_set(m, R["R3"], [1.0], **KW), theta)
    k1 = tolerant_set(m, R["R3"], [1.0], k=1, **KW)
    ct = cell_tolerant_set(m, R["R3"], [1.0], alpha=0.05, **KW)
    out["K1"] = rec(k1, theta)
    out["cell tolerant"] = rec(ct, theta)
    out["default"] = rec(_union([k1, ct]), theta)
    return out


JOBS = ([("ablation", sc, 400) for sc in SCEN15] +
        [("groups", sc, 200) for sc in SCEN15] +
        [("heldout", sc, 200) for sc in HELD_ALL] +
        [("groups400", sc, 400) for sc in ("contradictory evidence", "reference")])
FUN = dict(ablation=job_ablation, groups=job_groups, heldout=job_heldout, groups400=job_groups)


def fname(kind, sc):
    return os.path.join(OUT, f"{kind}__{sc.replace(' ', '_').replace(',', '').replace(':', '').replace('%', 'pct')}.json")


def _clear_stale_locks():
    for f in os.listdir(OUT):
        if f.endswith(".lock"):
            try:
                pid = int(open(os.path.join(OUT, f)).read().strip() or 0)
                os.kill(pid, 0)                      # raises if that process is gone
            except (ValueError, ProcessLookupError, PermissionError, OSError):
                os.remove(os.path.join(OUT, f))


def worker():
    os.makedirs(OUT, exist_ok=True)
    _clear_stale_locks()
    while True:
        todo = None
        for kind, sc, n in JOBS:
            fn = fname(kind, sc)
            if os.path.exists(fn):
                continue
            try:
                fd = os.open(fn + ".lock", os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                os.write(fd, str(os.getpid()).encode()); os.close(fd)
            except FileExistsError:
                continue
            todo = (kind, sc, n, fn)
            break
        if todo is None:
            return
        kind, sc, n, fn = todo
        t0 = time.time()
        start = 200 if kind == "groups400" else 0
        reps = [FUN[kind](sc, r) for r in range(start, n)]
        json.dump(dict(kind=kind, scenario=sc, n=len(reps), seconds=time.time() - t0, records=reps), open(fn + ".tmp", "w"))
        os.replace(fn + ".tmp", fn)
        os.remove(fn + ".lock")
        print(f"{time.strftime('%H:%M:%S')} done {kind} | {sc} | {len(reps)} literatures | {time.time() - t0:.0f}s", flush=True)


def run_comparators():
    os.makedirs(OUT, exist_ok=True)
    res = {}
    for sc in SCEN15:
        res[sc] = [comparators_one(sc, seed_of(sc, r)) for r in range(400)]
    json.dump(dict(kind="comparators", n=400, records=res), open(os.path.join(OUT, "comparators__15_scenarios.json"), "w"))
    sw = {}
    for fi, fam in enumerate(B2.SWEEP_FAMILIES):
        for mi, mult in enumerate(B2.SWEEP_MULTS):
            key = f"{fam} | x{mult}"
            recs = []
            for r in range(200):
                seed = 500000 + 100000 * fi + 10000 * mi + r
                rng = np.random.default_rng(seed)
                theta, b, present, studies = B2.literature_sweep(fam, mult, rng)
                main = [d for d in studies if d["K"] > 0]
                rv = re_variants([d["y"] for d in main], [d["s"] for d in main])
                cl = B2.cells(present, studies)
                o = dict(seed=int(seed))
                for nm in ("wald", "hk", "mhk"):
                    o["RE " + nm] = rec_interval(rv[nm][0], rv[nm][1], theta)
                o["BA-S independent"] = rec_interval(*ba_s(cl, False), theta)
                o["BA-S shared"] = rec_interval(*ba_s(cl, True), theta)
                recs.append(o)
            sw[key] = recs
    json.dump(dict(kind="comparators sweep", n=200, records=sw), open(os.path.join(OUT, "comparators__sweep.json"), "w"))
    print("comparators written")


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "worker":
        worker()
    elif mode == "comparators":
        run_comparators()
    elif mode == "status":
        for kind, sc, n in JOBS:
            fn = fname(kind, sc)
            st = "done" if os.path.exists(fn) else ("running" if os.path.exists(fn + ".lock") else "waiting")
            print(f"{st:8s} {kind:10s} {sc}")
