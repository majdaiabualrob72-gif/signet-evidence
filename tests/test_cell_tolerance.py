"""GENESIS VI checks of the reporting default (cell tolerance, Proposition S4). Run: python3 tests/test_cell_tolerance.py"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from signet.core import Model, signet_set, cell_tolerant_set, default_set
ok = True
def check(name, cond, val=""):
    global ok
    ok &= bool(cond); print(("PASS " if cond else "FAIL ") + name, val)
# 1. one design measuring the effect: omitting it leaves nothing, so the default must be unbounded
m = Model(1); m.add_row("cohort", [1.0], ["H"]); m.add_row("control", [0.0], ["H"])
m.sign("H down", "H", -1); m.share("H shared", "H", rung=2); m.set_data([-0.4, -0.1], [0.1, 0.1])
d = default_set(m, ["H down", "H shared"], [1.0])
check("single effect design gives an unbounded default", d["lo"] == -np.inf and d["hi"] == np.inf, (d["lo"], d["hi"]))
# 2. coverage with one broken cell: two cohorts with opposite signed biases and a trial; the trial is broken
def build():
    m = Model(1)
    m.add_row("cohort A", [1.0], ["H"]); m.add_row("cohort B", [1.0], ["R"]); m.add_row("trial", [1.0], ["T"])
    m.sign("H down", "H", -1); m.sign("R up", "R", +1); m.zero("trial valid", "T", rung=3)
    return m
rng = np.random.default_rng(5); theta = -0.2; s = np.array([0.08, 0.08, 0.06]); R = 400
cov_full = cov_cell = 0
for r in range(R):
    mu = np.array([theta - 0.10, theta + 0.10, theta + 0.35])      # the trial is shifted: broken cell
    y = mu + s * rng.standard_normal(3)
    m = build(); m.set_data(y, s); inc = ["H down", "R up", "trial valid"]
    f = signet_set(m, inc, [1.0]); c = cell_tolerant_set(m, inc, [1.0])
    cov_full += (not f["empty"]) and f["lo"] <= theta <= f["hi"]
    cov_cell += (not c["empty"]) and c["lo"] <= theta <= c["hi"]
check("cell-tolerant set covers with one broken cell (>= 0.93 in 400)", cov_cell / R >= 0.93, f"{cov_cell/R:.3f} (full analysis {cov_full/R:.3f})")
print("ALL PASS" if ok else "SOME FAILED")
