"""GENESIS VI: the archived simulation of the default (benchmark_v2.cell_robust_set) skipped the analysis that omits
a design when no other design measured the effect. Proposition S4 requires that analysis, which is unbounded, so the
correct default is unbounded in such literatures. This script regenerates the 200 archived literatures per scenario,
finds the affected ones, and reports the corrected coverage and benefit rates (exact for every literature: unaffected
literatures keep their archived result; affected ones are recomputed)."""
import os, sys, json
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SIM = os.path.join(ROOT, "experiments", "simulation"); sys.path.insert(0, ROOT); sys.path.insert(0, SIM); os.chdir(SIM)
import benchmark_v2 as B2
import benchmark as B
from signet.core import tolerant_set
arch = json.load(open("benchmark_v2_cellrobust.json"))["results"]
KEY = "SIGNET R3 k=1 (group or cell)"
jobs = [(sc, [91000 + 1000 * i + r for r in range(200)]) for i, sc in enumerate(B2.STRESS)]
jobs += [(sc, [7000 * (i + 1) + r for r in range(200)]) for i, sc in enumerate(B.ALL_SCEN)]
kw = dict(exact_rank=False, ngrid=24, tol=5e-4)
lines = ["| Scenario | Affected literatures | Archived coverage | Corrected coverage | Archived benefit | Corrected benefit | Group-only benefit lost |", "|---|---|---|---|---|---|---|"]
summary = {}
for sc, seeds in jobs:
    n_aff = cov_b = ben_b = gben = 0
    for sd in seeds:
        rng = np.random.default_rng(sd)
        theta, b, present, studies = B2.literature_v2(sc, rng)
        eff = [d for d in present if d[1] > 0]
        if len(eff) != 1:
            continue
        n_aff += 1
        r_ = B2.cell_robust_set(present, studies, kw)                # archived behaviour for this literature
        cov_b += (not r_["empty"]) and r_["lo"] - 1e-9 <= theta <= r_["hi"] + 1e-9
        ben_b += (not r_["empty"]) and r_["hi"] < 0
        m, R = B.signet_model(present, studies)
        t = tolerant_set(m, R["R3"], [1.0], k=1, **kw)
        gben += (not t["empty"]) and t["hi"] < 0
    a = arch[sc][KEY]; n = a["n"]
    cov_c = (a["coverage"] * n - cov_b + n_aff) / n
    ben_c = (a["benefit"] * n - ben_b) / n
    summary[sc] = dict(affected=n_aff, coverage_archived=a["coverage"], coverage_corrected=cov_c,
                       benefit_archived=a["benefit"], benefit_corrected=ben_c, group_only_benefit_lost=gben)
    lines.append(f"| {sc} | {n_aff}/200 | {a['coverage']:.3f} | {cov_c:.3f} | {a['benefit']:.3f} | {ben_c:.3f} | {gben} |")
    print(lines[-1], flush=True)
open(os.path.join(ROOT, "audit", "cellrobust_edgecase_output.md"), "w").write("\n".join(lines) + "\n")
json.dump(summary, open(os.path.join(ROOT, "audit", "cellrobust_edgecase_output.json"), "w"), indent=1)
