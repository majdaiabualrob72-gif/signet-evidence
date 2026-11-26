"""
benchmark.py
Synthetic literatures with known truth. Each literature mixes up to nine
designs whose biases have known directions (healthy-user and protopathic
bias favour the drug; channelling in registries opposes it), unknown-sign
confounding in active-comparator cohorts, and possible failures of trial,
surrogate and Mendelian-randomisation validity. Methods: random-effects
pooling of all designs, pooling of trials only, DECODE (sparse Bayesian
decoding with magnitude priors) and SIGNET at three assumption rungs, with
and without tolerance for one failed assumption.
Usage: python3 benchmark.py [nrep]
"""
import json, sys, time
import numpy as np
from multiprocessing import Pool
from decode_core import Ledger, decode, re_meta
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
sys.path.insert(0, _os.path.abspath(_os.path.join(_HERE, '..', '..')))
sys.path.insert(0, _HERE)
from signet.core import *

TH = float(np.log(0.75))
DESIGNS = [  # name, K, processes, presence probability, study SE range
    ("NU-U", 1.0, ["H", "P"], 0.9, (0.06, 0.16)),
    ("NU-L", 1.0, ["H", "P"], 0.7, (0.06, 0.16)),
    ("NU-NCO", 0.0, ["H"], 0.6, (0.06, 0.16)),
    ("AC", 1.0, ["C"], 0.7, (0.06, 0.16)),
    ("AC-NCO", 0.0, ["C"], 0.5, (0.06, 0.16)),
    ("REG", 1.0, ["R"], 0.5, (0.08, 0.18)),
    ("RCT", 1.0, ["A"], 0.6, (0.10, 0.25)),
    ("SUR", 1.0, ["S"], 0.5, (0.05, 0.10)),
    ("MR", 1.0, ["Pl"], 0.5, (0.10, 0.20)),
]
PROCS = ["H", "P", "C", "R", "A", "S", "Pl"]
DEC_PI = [0.9, 0.9, 0.5, 0.9, 0.15, 0.3, 0.2]
DEC_SIG = [0.15, 0.18, 0.06, 0.12, 0.12, 0.20, 0.15]
TAU = np.array([0.0, 0.03, 0.06, 0.10, 0.16])
SCEN = ["reference", "null effect", "magnitudes doubled", "dense small biases",
        "15% of signs wrong", "omitted bias", "no design opposing benefit", "no bias"]


def draw_biases(sc, rng):
    mult = 2.0 if sc in ("magnitudes doubled", "one-signed literature, null effect, magnitudes doubled") else 1.0
    b = dict(H=-abs(rng.normal(0.12, 0.04)) * mult, P=-abs(rng.normal(0.15, 0.05)) * mult,
             R=abs(rng.normal(0.10, 0.04)) * mult, C=rng.normal(0, 0.04),
             A=0.0 if rng.random() < 0.85 else rng.normal(0, 0.12 * mult),
             S=0.0 if rng.random() < 0.70 else rng.normal(0, 0.20 * mult),
             Pl=0.0 if rng.random() < 0.80 else rng.normal(0, 0.15 * mult))
    if sc == "15% of signs wrong":
        for p in ("H", "P", "R"):
            if rng.random() < 0.15:
                b[p] = -b[p]
    if sc == "no bias":
        b = {p: 0.0 for p in PROCS}
    return b


def literature(sc, rng):
    theta = 0.0 if sc in ("null effect", "one-signed literature, null effect",
                          "one-signed literature, null effect, magnitudes doubled") else TH
    b = draw_biases(sc, rng)
    present = []
    while not any(K > 0 for (_n, K, *_r) in present):
        present = [d for d in DESIGNS if rng.random() < d[3]]
        if sc == "no design opposing benefit":
            present = [d for d in present if d[0] in ("NU-U", "NU-L", "NU-NCO", "AC", "AC-NCO")]
        if sc.startswith("one-signed literature"):
            present = [d for d in present if d[0] in ("NU-U", "NU-L")]
    studies = []
    for (name, K, procs, _pr, (lo, hi)) in present:
        k = 1 + rng.poisson(2)
        mu = K * theta + sum(b[p] for p in procs if not (p == "P" and name == "NU-L"))
        if sc == "omitted bias" and name in ("NU-U", "NU-L", "AC", "REG"):
            mu += 0.10
        for _ in range(k):
            se = rng.uniform(lo, hi)
            extra = -abs(rng.normal(0, 0.05)) if (sc == "dense small biases" and name in
                                                   ("NU-U", "NU-L", "NU-NCO", "AC", "AC-NCO", "REG")) else 0.0
            studies.append(dict(design=name, K=K, procs=procs, y=mu + extra + rng.normal(0, se), s=se))
    return theta, b, present, studies


def signet_model(present, studies):
    m = Model(1)
    ys, ss = [], []
    for (name, K, procs, *_r) in present:
        st = [d for d in studies if d["design"] == name]
        y_, s_, _ = pool_cell([d["y"] for d in st], [d["s"] for d in st])
        m.add_row(name, [K], procs); ys.append(y_); ss.append(s_)
    m.set_data(ys, ss)
    names = [r["name"] for r in m.rows]
    groups = {1: [], 2: [], 3: []}

    def has(p):
        return len(m.rows_with(p)) > 0
    if has("H"):
        m.sign("H-", "H", -1, rung=1); groups[1].append("H-")
    if has("P"):
        m.sign("P-", "P", -1, rung=1); groups[1].append("P-")
        if "NU-L" in names:
            m.zero("lagP0", "P", ["NU-L"], rung=1); groups[1].append("lagP0")
    if has("R"):
        m.sign("R+", "R", +1, rung=1); groups[1].append("R+")
    if len(m.rows_with("H")) > 1:
        m.share("shH", "H", rung=2); groups[2].append("shH")
    if len(m.rows_with("C")) > 1:
        m.share("shC", "C", rung=2); groups[2].append("shC")
    for p, g in (("A", "vA"), ("S", "vS"), ("Pl", "vPl")):
        if has(p):
            m.zero(g, p, rung=3); groups[3].append(g)
    r1 = groups[1]; r2 = r1 + groups[2]; r3 = r2 + groups[3]
    return m, dict(R1=r1, R2=r2, R3=r3)


def summarise_set(lo, hi, theta, empty=False):
    if empty or lo is None:
        return dict(cover=False, width=None, unbounded=False, empty=True, benefit=False, harm=False)
    w = hi - lo
    return dict(cover=bool(lo <= theta <= hi), width=float(w) if np.isfinite(w) else None,
                unbounded=bool(not np.isfinite(w)), empty=False, benefit=bool(hi < 0), harm=bool(lo > 0))


def one(args):
    sc, seed = args
    rng = np.random.default_rng(seed)
    theta, b, present, studies = literature(sc, rng)
    out = dict(theta=theta)
    main = [d for d in studies if d["K"] > 0]
    ma = re_meta([d["y"] for d in main], [d["s"] for d in main])
    out["RE (all designs)"] = summarise_set(ma["lo"], ma["hi"], theta)
    tr = [d for d in studies if d["design"] == "RCT"]
    if tr:
        mt = re_meta([d["y"] for d in tr], [d["s"] for d in tr])
        out["RE (trials only)"] = summarise_set(mt["lo"], mt["hi"], theta)
    else:
        out["RE (trials only)"] = dict(cover=False, width=None, unbounded=True, empty=False, benefit=False, harm=False, absent=True)
    # DECODE at study level
    Z = np.array([[d["K"]] for d in studies])
    S = np.zeros((len(studies), len(PROCS)))
    for i, d in enumerate(studies):
        for p in d["procs"]:
            if p == "P" and d["design"] == "NU-L":
                continue
            S[i, PROCS.index(p)] = 1.0
    used = S.sum(axis=0) > 0
    L = Ledger(np.array([d["y"] for d in studies]), np.array([d["s"] for d in studies]), Z, S[:, used],
               np.array(DEC_PI)[used], np.array(DEC_SIG)[used])
    post = decode(L, tau_grid=TAU)
    lo, hi = post.theta(0).interval()
    out["DECODE"] = summarise_set(lo, hi, theta)
    # SIGNET
    m, R = signet_model(present, studies)
    kw = dict(exact_rank=False, ngrid=24, tol=5e-4, outer_scan=0)
    for rn in ("R1", "R2", "R3"):
        r_ = signet_set(m, R[rn], [1.0], **kw)
        out[f"SIGNET {rn}"] = summarise_set(r_["lo"], r_["hi"], theta, r_["empty"])
        bx = box_set(m, R[rn], [1.0])
        out[f"box {rn}"] = summarise_set(bx.get("lo"), bx.get("hi"), theta, bx["empty"])
        cb = cert_bounds(m, R[rn], [1.0], a=0.025)
        out[f"cert {rn}"] = summarise_set(cb["lo"], cb["hi"], theta)
    cc = conf_set(m, R["R3"], [1.0], 0.05, **kw)
    out["cc R3"] = summarise_set(cc.get("lo"), cc.get("hi"), theta, cc["empty"])
    tk = tolerant_set(m, R["R3"], [1.0], k=1, **kw)
    out["SIGNET R3 k=1"] = summarise_set(tk["lo"], tk["hi"], theta, tk["empty"])
    mt3 = model_test(m, R["R3"], exact_rank=False)
    out["R3 rejected"] = mt3["reject"]
    # the protocol of record: test rung 3; adopt a unique single repair, else fall back to rung 2
    if not mt3["reject"]:
        work = R["R3"]; out["protocol rung"] = "R3"
    else:
        rp = repairs(m, R["R3"], kmax=1, droppable=[g for g in R["R3"] if g not in R["R1"]])
        if len(rp) == 1 and len(rp[0]) == 1:
            work = [g for g in R["R3"] if g not in rp[0]]; out["protocol rung"] = "R3 repaired"
        else:
            work = R["R2"]; out["protocol rung"] = "R2"
    r_ = signet_set(m, work, [1.0], **kw) if work != R["R3"] else None
    out["SIGNET protocol"] = out["SIGNET R3"] if r_ is None else summarise_set(r_["lo"], r_["hi"], theta, r_["empty"])
    out["R3 has failure"] = bool(any(abs(b[p]) > 0 for p in ("A", "S", "Pl") if len(m.rows_with(p))))
    return out


METHODS = ["RE (all designs)", "RE (trials only)", "DECODE", "SIGNET R1", "SIGNET R2", "SIGNET R3",
           "SIGNET protocol", "SIGNET R3 k=1", "box R1", "box R2", "box R3", "cert R1", "cert R2", "cert R3", "cc R3"]


def summarise(reps):
    s = {}
    for mth in METHODS:
        rr = [r[mth] for r in reps if not r[mth].get("absent")]
        w = [x["width"] for x in rr if x["width"] is not None]
        th = reps[0]["theta"]
        s[mth] = dict(n=len(rr), coverage=float(np.mean([x["cover"] for x in rr])) if rr else None,
                      median_width=float(np.median(w)) if w else None,
                      unbounded=float(np.mean([x["unbounded"] for x in rr])) if rr else None,
                      empty=float(np.mean([x["empty"] for x in rr])) if rr else None,
                      benefit=float(np.mean([x["benefit"] for x in rr])) if rr else None,
                      harm=float(np.mean([x["harm"] for x in rr])) if rr else None)
    s["R3 rejection rate"] = float(np.mean([r["R3 rejected"] for r in reps]))
    s["protocol rungs"] = {k: float(np.mean([r["protocol rung"] == k for r in reps])) for k in ("R3", "R3 repaired", "R2")}
    f = [r for r in reps if r["R3 has failure"]]
    s["R3 rejection when a validity zero failed"] = float(np.mean([r["R3 rejected"] for r in f])) if f else None
    nf = [r for r in reps if not r["R3 has failure"]]
    s["R3 rejection when none failed"] = float(np.mean([r["R3 rejected"] for r in nf])) if nf else None
    return s


def fmt(v, d=3):
    return "NA" if v is None else f"{v:.{d}f}"


ALL_SCEN = SCEN + ["one-signed literature, null effect", "one-signed literature, null effect, magnitudes doubled"]

if __name__ == "__main__":
    # usage: benchmark.py nrep [index ...]; results are merged into benchmark_results.json
    nrep = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    idx = [int(a) for a in sys.argv[2:]] or list(range(len(ALL_SCEN)))
    t0 = time.time()
    try:
        res = json.load(open("benchmark_results.json"))["results"]
    except FileNotFoundError:
        res = {}
    with Pool(2) as pool:
        for i in idx:
            sc = ALL_SCEN[i]
            try:
                reps = json.load(open(f"benchmark_reps_{i}.json"))
                assert len(reps) == nrep
            except (FileNotFoundError, AssertionError):
                reps = pool.map(one, [(sc, 7000 * (i + 1) + r) for r in range(nrep)], chunksize=4)
                json.dump(reps, open(f"benchmark_reps_{i}.json", "w"))
            res[sc] = summarise(reps)
            json.dump(dict(nrep=nrep, results=res), open("benchmark_results.json", "w"), indent=1)
            print(sc, f"{time.time() - t0:.0f}s", flush=True)
            for mth in METHODS:
                v = res[sc][mth]
                print(f"   {mth:<18} cov {fmt(v['coverage'])}  width {fmt(v['median_width'])}"
                      f"  unb {fmt(v['unbounded'], 2)}  empty {fmt(v['empty'], 2)}  benefit {fmt(v['benefit'], 2)}  harm {fmt(v['harm'], 2)}", flush=True)
