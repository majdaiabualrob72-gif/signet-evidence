"""GENESIS VI: frequency of accepted values outside the reported conditional chi-squared hull (40-point grid),
and whether outer_scan=40 captures them. Worked example, true assumption set, Gaussian cells with known variances."""
import os, sys
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "experiments", "simulation"))
from signet.core import *            # noqa
import worked_example as W           # noqa
m = W.build(); y0, s0, _, info = W.simulate(np.random.default_rng(W.SEED)); V0 = np.diag(s0 ** 2)
TRUE = [g for g in W.RUNG3 if g != "Surrogate valid"]
mu = np.array([d["mu_true"] for d in info]); cf = W.FUNCS["Early-stage population"]; c = m.cvec(cf)
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 777)
N = int(sys.argv[2]) if len(sys.argv) > 2 else 120
n = found = captured = 0
for r in range(N):
    y = mu + s0 * rng.standard_normal(len(s0))
    a = conf_set(m, TRUE, cf, 0.025, y, V0)
    if a.get("lo") is None or not (np.isfinite(a["lo"]) and np.isfinite(a["hi"])): continue
    n += 1; F = Fitter(m, TRUE, y, V0); w = a["hi"] - a["lo"]
    g = np.linspace(a["lo"] - 3 * w, a["hi"] + 3 * w, 300)
    out = [t for t in g if (t < a["lo"] - 1e-6 or t > a["hi"] + 1e-6) and F.accept(c, t, 0.025)[0]]
    if out:
        found += 1
        b = conf_set(m, TRUE, cf, 0.025, y, V0, outer_scan=40)
        captured += all(b["lo"] - 1e-6 <= t <= b["hi"] + 1e-6 for t in out)
        print(f"lit {r}: hull {a['lo']:.3f} to {a['hi']:.3f}; accepted outside {min(out):.3f} to {max(out):.3f}; outer_scan hull {b['lo']:.3f} to {b['hi']:.3f}", flush=True)
msg = f"island check (seed {sys.argv[1] if len(sys.argv)>1 else 777}): {n} literatures, accepted values outside the 40-point hull in {found}, captured by outer_scan=40 in {captured}"
print(msg); open(os.path.join(ROOT, "audit", "island_check_output.txt"), "a").write(msg + "\n")
