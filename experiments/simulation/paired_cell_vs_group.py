"""Paired comparison on the same simulated literatures: tolerance to one failed group (archived per-literature
outputs, first 200 literatures of each scenario) against tolerance to one failed group or one broken cell
(benchmark_v2_cellrobust.json, the same 200 seeds). Writes paired_cell_vs_group.md."""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import benchmark as B                     # noqa: E402

STRESS = ["extreme outlier", "publication bias", "heterogeneous populations",
          "missing controls and variances", "contradictory evidence"]
cr = json.load(open("benchmark_v2_cellrobust.json"))["results"]
rows = ["| Scenario | Group: coverage | Group or cell: coverage | Group: benefit | Group or cell: benefit |",
        "|---|---|---|---|---|"]
for i, sc in enumerate(B.ALL_SCEN):
    fn = f"benchmark_reps_{i}.json"
    reps = json.load(open(fn))[:200]
    g = [r["SIGNET R3 k=1"] for r in reps]
    c = cr.get(sc, {}).get("SIGNET R3 k=1 (group or cell)")
    if c is None:
        continue
    rows.append(f"| {sc} | {sum(x['cover'] for x in g) / 200:.3f} | {c['coverage']:.3f} | "
                f"{sum(x['benefit'] for x in g) / 200:.3f} | {c['benefit']:.3f} |")
for i, sc in enumerate(STRESS):
    reps = json.load(open(f"benchmark_v2_stress_reps_{i}.json"))[:200]
    g = [r["SIGNET R3 k=1"] for r in reps]
    c = cr.get(sc, {}).get("SIGNET R3 k=1 (group or cell)")
    if c is None:
        continue
    rows.append(f"| {sc} | {sum(x['cover'] for x in g) / 200:.3f} | {c['coverage']:.3f} | "
                f"{sum(x['benefit'] for x in g) / 200:.3f} | {c['benefit']:.3f} |")
txt = "\n".join(rows) + "\n"
open("paired_cell_vs_group.md", "w").write("# Same 200 literatures per scenario: group tolerance vs group-or-cell tolerance\n\n" + txt)
print(txt)
