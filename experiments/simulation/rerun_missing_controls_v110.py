"""
rerun_missing_controls_v110.py: regenerates the 'missing controls and variances' scenario after the
v1.1.0 fix to its generator (benchmark_v2.literature_v2 iterated over a Python set of design names
when imputing variances, so the random draws depended on the per-process string-hash order and the
scenario was not reproducible literature by literature). Same seeds and settings as the archived runs.
The stress part is re-run with `python benchmark_v2.py stress 400 3`; this script re-runs the
magnitude rung (benchmark_v2_bounded.json, 400) and the cell-level redesign
(benchmark_v2_cellrobust.json, 200), merges them, refreshes all_scenarios_summary.json from
benchmark_v2_stress.json, and replaces the scenario in ../validation/reps_cellrobust.json.
"""
import json
import os
import sys

import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
sys.path.insert(0, HERE)
os.chdir(HERE)
import benchmark_v2 as B2   # noqa: E402

SC = "missing controls and variances"
I = B2.STRESS.index(SC)

if __name__ == "__main__":
    with Pool(2) as pool:
        reps_b = pool.map(B2.one_bounded, [(SC, 91000 + 1000 * I + r) for r in range(400)], chunksize=4)
        reps_c = pool.map(B2.one_cellrobust, [(SC, 91000 + 1000 * I + r) for r in range(200)], chunksize=4)
    d = json.load(open("benchmark_v2_bounded.json"))
    d["results"][SC] = B2.summarise(reps_b, B2.METHODS_B)
    d["results"][SC]["median seconds"] = float(np.median([r["t_bounded"] for r in reps_b]))
    json.dump(d, open("benchmark_v2_bounded.json", "w"), indent=1)
    d = json.load(open("benchmark_v2_cellrobust.json"))
    d["results"][SC] = B2.summarise(reps_c, B2.METHODS_CR)
    d["results"][SC]["median seconds"] = float(np.median([r["t_cellrobust"] for r in reps_c]))
    json.dump(d, open("benchmark_v2_cellrobust.json", "w"), indent=1)
    st = json.load(open("benchmark_v2_stress.json"))["results"][SC]
    a = json.load(open("all_scenarios_summary.json"))
    a[SC] = st
    json.dump(a, open("all_scenarios_summary.json", "w"), indent=1)
    v = os.path.join(HERE, "..", "validation", "reps_cellrobust.json")
    if os.path.exists(v):
        r = json.load(open(v)); r[SC] = reps_c; json.dump(r, open(v, "w"))
    for name, reps, key in (("bounded", reps_b, "SIGNET-B k=1"), ("cellrobust", reps_c, "SIGNET R3 k=1 (group or cell)")):
        print(name, "coverage", np.mean([x[key]["cover"] for x in reps]), "benefit", np.mean([x[key]["benefit"] for x in reps]), flush=True)
