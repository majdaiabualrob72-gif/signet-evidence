"""
validation_v11.py: simulation checks requested by the November 2026 audit (signet v1.1.0).

On the archived simulated literatures of the main benchmark (benchmark.py; same generator and
seeds 7000 (i + 1) + r), SIGNET rung 1 and rung 3 sets are recomputed under:
  legacy    the archived settings (estimated V, approximate rank, 24-point grid, tol 5e-4,
            no dense check); must reproduce benchmark_reps_{i}.json
  hardened  exact rank, 40-point grid, tol 1e-4, 16-point dense check, bracketing to 1e4
  outer     the n-degree outer set (exact end points, conservative by construction)
  frozen    certificate weights chosen with the within-study (fixed-effect) variances, which do
            not depend on the effect estimates; everything else as legacy
  knownV    cell estimates pooled with the true heterogeneity variance and analysed with
            their exact sampling covariance (the setting of the theory); otherwise as legacy
Per replicate it also records whether each certificate is unique (then its choice cannot depend
on V or y), whether the dense check extended a hull, and whether a v1.0.0 search would have
stopped before bracketing an outer limit.

Usage: python3 experiments/validation/validation_v11.py NREP [scenario index ...]
Output: experiments/validation/validation_reps_{i}.json (one record per replicate)
"""
import json
import os
import sys
import time

import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SIM = os.path.join(ROOT, "experiments", "simulation")
sys.path.insert(0, ROOT)
sys.path.insert(0, SIM)
import benchmark as B                                    # noqa: E402
from signet.core import Model, signet_set, certificate_range   # noqa: E402

LEGACY = dict(exact_rank=False, ngrid=24, tol=5e-4, outer_scan=0)
HARD = dict(exact_rank=True, ngrid=40, tol=1e-4, outer_scan=16)
DENSE = ("NU-U", "NU-L", "NU-NCO", "AC", "AC-NCO", "REG")
TAU2_DENSE = 0.05 ** 2 * (1 - 2 / np.pi)       # variance of -|N(0, 0.05^2)|


def fe_var(studies, name):
    s = np.array([d["s"] for d in studies if d["design"] == name])
    return float(1 / np.sum(1 / s ** 2))


def known_model(sc, present, studies):
    """Cells pooled with weights 1 / (s^2 + tau2_true); V is their exact sampling covariance."""
    m0, R = B.signet_model(present, studies)       # same assumption groups and rungs
    ys, vs = [], []
    for (name, K, procs, *_r) in present:
        st = [d for d in studies if d["design"] == name]
        t2 = TAU2_DENSE if (sc == "dense small biases" and name in DENSE) else 0.0
        w = 1 / (np.array([d["s"] for d in st]) ** 2 + t2)
        ys.append(float(np.sum(w * np.array([d["y"] for d in st])) / np.sum(w)))
        vs.append(float(1 / np.sum(w)))
    m0.set_data(ys, np.sqrt(vs))
    return m0, R


def rec(s, theta):
    if s["empty"]:
        return dict(lo=None, hi=None, cover=False, empty=True)
    lo, hi = s["lo"], s["hi"]
    return dict(lo=lo if np.isfinite(lo) else str(lo), hi=hi if np.isfinite(hi) else str(hi),
                cover=bool(lo <= theta <= hi), empty=False, benefit=bool(hi < 0), harm=bool(lo > 0))


def one(args):
    sc, seed = args
    rng = np.random.default_rng(seed)
    theta, b, present, studies = B.literature(sc, rng)
    m, R = B.signet_model(present, studies)
    out = dict(seed=seed, theta=theta, n_cells=m.n)
    t0 = time.perf_counter()
    Vfe = np.diag([fe_var(studies, r["name"]) for r in m.rows])
    mk, _ = known_model(sc, present, studies)
    for rn in ("R1", "R3"):
        inc = R[rn]
        out[f"legacy {rn}"] = rec(signet_set(m, inc, [1.0], **LEGACY), theta)
        out[f"frozen {rn}"] = rec(signet_set(m, inc, [1.0], cert_V=Vfe, **LEGACY), theta)
        out[f"knownV {rn}"] = rec(signet_set(mk, inc, [1.0], **LEGACY), theta)
        cu = {}
        for side in ("hi", "lo"):
            cr = certificate_range(m, inc, [1.0], side)
            cu[side] = None if cr is None else cr["unique"]
        out[f"cert unique {rn}"] = cu
    h = signet_set(m, R["R3"], [1.0], **HARD)
    out["hardened R3"] = rec(h, theta)
    cc = h.get("cc", {})
    out["hardened flags"] = dict(hull_extended=cc.get("hull_extended"), bracketed=cc.get("bracketed"),
                                 legacy_truncation=cc.get("legacy_truncation"), warnings=cc.get("warnings"))
    out["outer R3"] = rec(signet_set(m, R["R3"], [1.0], inversion="outer", tol=5e-4), theta)
    out["seconds"] = time.perf_counter() - t0
    return out


if __name__ == "__main__":
    nrep = int(sys.argv[1])
    idx = [int(a) for a in sys.argv[2:]] or list(range(len(B.ALL_SCEN)))
    with Pool(2) as pool:
        for i in idx:
            sc = B.ALL_SCEN[i]
            fn = os.path.join(HERE, f"validation_reps_{i}.json")
            t = time.time()
            reps = pool.map(one, [(sc, 7000 * (i + 1) + r) for r in range(nrep)], chunksize=4)
            json.dump(dict(scenario=sc, nrep=nrep, reps=reps), open(fn, "w"))
            print(sc, f"{time.time() - t:.0f}s", flush=True)
