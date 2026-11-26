"""
validation_default_v11.py: the reporting default (tolerant to one failed assumption group or one
broken evidence cell) recomputed with the archived settings and with the hardened settings of
signet v1.1.0 (exact rank, 40-point grid, tol 1e-4, 16-point dense check), on the first NREP
archived literatures of each of the 15 scenarios (base seeds 7000 (i + 1) + r; stress seeds
91000 + 1000 i + r).

Usage: python3 experiments/validation/validation_default_v11.py NREP
Output: experiments/validation/validation_default_reps.json
"""
import json
import os
import sys
import time

import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "simulation"))
import benchmark as B          # noqa: E402
import benchmark_v2 as B2      # noqa: E402

LEGACY = dict(exact_rank=False, ngrid=24, tol=5e-4, outer_scan=0)
HARD = dict(exact_rank=True, ngrid=40, tol=1e-4, outer_scan=16)


def rec(r, theta):
    lo, hi = r["lo"], r["hi"]
    if r["empty"] or lo is None:
        return dict(cover=False, empty=True, lo=None, hi=None, benefit=False)
    f = lambda v: float(v) if np.isfinite(v) else str(v)
    return dict(cover=bool(lo <= theta <= hi), empty=False, lo=f(lo), hi=f(hi), benefit=bool(hi < 0), harm=bool(lo > 0))


def one(args):
    sc, seed = args
    rng = np.random.default_rng(seed)
    theta, b, present, studies = B2.literature_v2(sc, rng)
    t0 = time.perf_counter()
    a = B2.cell_robust_set(present, studies, LEGACY)
    t1 = time.perf_counter()
    h = B2.cell_robust_set(present, studies, HARD)
    t2 = time.perf_counter()
    return dict(seed=seed, theta=theta, legacy=rec(a, theta), hardened=rec(h, theta), t_legacy=t1 - t0, t_hard=t2 - t1)


if __name__ == "__main__":
    nrep = int(sys.argv[1])
    jobs = [(sc, [7000 * (i + 1) + r for r in range(nrep)]) for i, sc in enumerate(B.ALL_SCEN)]
    jobs += [(sc, [91000 + 1000 * i + r for r in range(nrep)]) for i, sc in enumerate(B2.STRESS)]
    res = {}
    t = time.time()
    with Pool(2) as pool:
        for sc, seeds in jobs:
            res[sc] = pool.map(one, [(sc, s) for s in seeds], chunksize=2)
            json.dump(res, open(os.path.join(HERE, "validation_default_reps.json"), "w"))
            print(sc, f"{time.time() - t:.0f}s", flush=True)
