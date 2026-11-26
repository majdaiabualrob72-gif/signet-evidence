"""
paired_reruns.py: re-runs of archived simulation parts that kept only summaries, now saving one
record per simulated literature so that methods can be compared in pairs (McNemar).

Each part uses the archived generator, settings and seeds, so its summaries must reproduce the
archived files (checked by paired_contrasts.py):
  cellrobust   benchmark_v2.one_cellrobust, 15 scenarios x 200 (benchmark_v2_cellrobust.json)
  comparators  benchmark_v2.one_comparators, 10 scenarios x 400 (benchmark_v2_comparators.json)
  sweep        benchmark_v2.one_sweep, the two cells used in the contrasts, 200 each
  hetallow     benchmark_v2.one_hetallow, heterogeneous populations, 400 (tau_h = 0.10)
  bounded      benchmark_v2.one_bounded, reference, 400

Usage: python3 experiments/validation/paired_reruns.py PART [PART ...]
Output: experiments/validation/reps_{PART}.json
"""
import json
import os
import sys
import time

from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "simulation"))
import benchmark as B          # noqa: E402
import benchmark_v2 as B2      # noqa: E402


def jobs(part):
    if part == "cellrobust":
        J = [(sc, [91000 + 1000 * i + r for r in range(200)]) for i, sc in enumerate(B2.STRESS)]
        J += [(sc, [7000 * (i + 1) + r for r in range(200)]) for i, sc in enumerate(B.ALL_SCEN)]
        return B2.one_cellrobust, [(sc, [(sc, s) for s in seeds]) for sc, seeds in J]
    if part == "comparators":
        return B2.one_comparators, [(sc, [(sc, 7000 * (i + 1) + r) for r in range(400)]) for i, sc in enumerate(B.ALL_SCEN)]
    if part == "sweep":
        out = []
        for fi, fam in enumerate(B2.SWEEP_FAMILIES):
            mult = 4.0 if fi == 0 else 3.0
            mi = B2.SWEEP_MULTS.index(mult)
            out.append((f"{fam} | x{mult}", [(fam, mult, 500000 + 100000 * fi + 10000 * mi + r) for r in range(200)]))
        return B2.one_sweep, out
    if part == "hetallow":
        return B2.one_hetallow, [("heterogeneous populations", [("heterogeneous populations", 91000 + 1000 * 2 + r) for r in range(400)])]
    if part == "bounded":
        return B2.one_bounded, [("reference", [("reference", 7000 + r) for r in range(400)])]
    raise SystemExit("unknown part " + part)


if __name__ == "__main__":
    for part in sys.argv[1:]:
        fun, J = jobs(part)
        res = {}
        t = time.time()
        with Pool(2) as pool:
            for key, args in J:
                res[key] = pool.map(fun, args, chunksize=4)
                json.dump(res, open(os.path.join(HERE, f"reps_{part}.json"), "w"))
                print(part, key, f"{time.time() - t:.0f}s", flush=True)
