"""GENESIS VI: does the package's default_set reproduce the simulation script's group-or-cell set?
Compares signet.core.default_set (row omission via _OmitRows) with benchmark_v2.cell_robust_set (models rebuilt
without each design) on archived contradictory-evidence and reference literatures."""
import os, sys
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SIM = os.path.join(ROOT, "experiments", "simulation"); sys.path.insert(0, ROOT); sys.path.insert(0, SIM); os.chdir(SIM)
import benchmark_v2 as B2
import benchmark as B
from signet.core import default_set
kw = dict(exact_rank=False, ngrid=24, tol=5e-4)
worst = 0.0; n = same_claim = same_cov = 0
for sc, base in (("contradictory evidence", 91000 + 1000 * B2.STRESS.index("contradictory evidence")), ("reference", 7000 * (B.ALL_SCEN.index("reference") + 1))):
    for r in range(12):
        rng = np.random.default_rng(base + r)
        theta, b, present, studies = B2.literature_v2(sc, rng)
        a = B2.cell_robust_set(present, studies, kw)
        m, R = B.signet_model(present, studies)
        d = default_set(m, R["R3"], [1.0], **kw)
        n += 1
        if a["empty"] or d["empty"]:
            same_claim += a["empty"] == d["empty"]; same_cov += a["empty"] == d["empty"]; continue
        for u, v in ((a["lo"], d["lo"]), (a["hi"], d["hi"])):
            if np.isfinite(u) and np.isfinite(v): worst = max(worst, abs(u - v))
            elif u != v: worst = np.inf
        same_claim += (a["hi"] < 0) == (d["hi"] < 0)
        same_cov += (a["lo"] <= theta <= a["hi"]) == (d["lo"] <= theta <= d["hi"])
msg = f"default_set vs simulation script on {n} archived literatures: largest difference in limits {worst:.2e}; same benefit claim {same_claim}/{n}; same coverage {same_cov}/{n}"
print(msg); open(os.path.join(ROOT, "audit", "default_equivalence_output.txt"), "w").write(msg + "\n")
