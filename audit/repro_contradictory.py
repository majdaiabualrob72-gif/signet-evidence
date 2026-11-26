"""GENESIS VI: re-run the default (group-or-cell) set on the 200 contradictory-evidence literatures with the archived
seeds, check that the archived coverage (0.925) is reproduced, and record whether each non-covering literature also
contained a failed validity assumption (trial, surrogate or Mendelian randomisation), which the default does not cover."""
import os, sys, json
from multiprocessing import Pool
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SIM = os.path.join(ROOT, "experiments", "simulation"); sys.path.insert(0, ROOT); sys.path.insert(0, SIM)
os.chdir(SIM)
import benchmark_v2 as B2
import benchmark as B
i = B2.STRESS.index("contradictory evidence")
seeds = [91000 + 1000 * i + r for r in range(200)]
def one(seed):
    rng = np.random.default_rng(seed)
    theta, b, present, studies = B2.literature_v2("contradictory evidence", rng)
    r_ = B2.cell_robust_set(present, studies, dict(exact_rank=False, ngrid=24, tol=5e-4))
    cov = (not r_["empty"]) and r_["lo"] - 1e-9 <= theta <= r_["hi"] + 1e-9
    names = [d[0] for d in present]
    procs = set(p for d in present for p in d[2]) if len(present[0]) > 2 else set()
    vfail = any(abs(b.get(p, 0.0)) > 0 for p in ("A", "S", "Pl") if any(p in d[2] for d in present))
    return dict(seed=seed, cov=bool(cov), empty=bool(r_["empty"]), vfail=bool(vfail))
if __name__ == "__main__":
    with Pool(2) as pool:
        res = pool.map(one, seeds, chunksize=5)
    cov = np.mean([r["cov"] for r in res]); nc = [r for r in res if not r["cov"]]
    vf_nc = sum(r["vfail"] for r in nc); vf_all = sum(r["vfail"] for r in res)
    cov_novf = np.mean([r["cov"] for r in res if not r["vfail"]]) if any(not r["vfail"] for r in res) else float("nan")
    msg = (f"contradictory evidence, default set, 200 archived seeds: coverage {cov:.3f} (archived 0.925); "
           f"non-covering literatures {len(nc)}, of which {vf_nc} also had a failed validity assumption; "
           f"literatures with a validity failure {vf_all}/200; coverage when no validity assumption failed {cov_novf:.3f}")
    print(msg); open(os.path.join(ROOT, "audit", "repro_contradictory_output.txt"), "w").write(msg + "\n")
