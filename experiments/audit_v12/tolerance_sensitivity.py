"""
tolerance_sensitivity.py  (N10; signet 1.2.0)

The three numerical thresholds of the inference code (signet.core.TOL) are absolute:
  T_zero    T(theta0) at or below this is treated as zero                         default 1e-9
  slack     slack below this makes a constraint a candidate active constraint     default 1e-4
  implicit  largest attainable slack below this makes it an implicit equality     default 1e-7
This script varies each by a factor of 100 in both directions and recomputes
  (a) every SIGNET set and directional p-value of the real-data series, and
  (b) rung-3 sets for the first 20 archived literatures of each of the 15 simulated scenarios, with the
      archived settings (single-programme rank) and with exact ranks;
and (c) rescales the worked example (estimates and standard errors multiplied by c), to show the range of
scales over which the absolute settings (thresholds, bracketing step and bisection tolerance) are adequate.
Output: tolerance_sensitivity.json and tolerance_sensitivity.md
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "experiments", "simulation"))
import signet.core as C                                   # noqa: E402
from signet.core import Model, signet_set, claim_pvalue   # noqa: E402
import benchmark as B                                     # noqa: E402
import benchmark_v2 as B2                                 # noqa: E402
import sim_audit as S                                     # noqa: E402
import realdata_audit_models as RM                        # noqa: E402

BASE = dict(C.TOL)
SETTINGS = [("default", {})]
for k in ("T_zero", "slack", "implicit"):
    for f in (1e-2, 1e2):
        SETTINGS.append((f"{k} x {f:g}", {k: BASE[k] * f}))


def with_tol(over):
    C.TOL.update(BASE); C.TOL.update(over)


def fin(v):
    return None if (v is None or not np.isfinite(v)) else float(v)


def real_all():
    out = {}
    for lab, (M, R), cf in RM.ANALYSES():
        for rn, inc in R.items():
            s = signet_set(M, inc, cf, 0.05)
            pb = claim_pvalue(M, inc, cf, "benefit", empty=bool(s["empty"]))["p"]
            ph = claim_pvalue(M, inc, cf, "harm", empty=bool(s["empty"]))["p"]
            out[f"{lab} {rn}"] = dict(empty=bool(s["empty"]), lo=fin(s["lo"]), hi=fin(s["hi"]), p_benefit=pb, p_harm=ph)
    return out


def sim_all(exact):
    kw = dict(S.KW); kw["exact_rank"] = exact
    out = {}
    for sc in S.SCEN15:
        for r in range(20):
            seed = S.seed_of(sc, r)
            rng = np.random.default_rng(seed)
            theta, b, present, studies = B2.literature_v2(sc, rng)
            m, R = B.signet_model(present, studies)
            s = signet_set(m, R["R3"], [1.0], **kw)
            out[f"{sc}|{r}"] = S.rec(s, theta)
    return out


def maxdiff(a, b, keys):
    worst, changed = 0.0, 0
    for k in a:
        for q in keys:
            x, y = a[k].get(q), b[k].get(q)
            if isinstance(x, bool) or isinstance(y, bool):
                changed += int(x != y); continue
            if (x is None) != (y is None):
                changed += 1; continue
            if x is not None:
                worst = max(worst, abs(x - y))
    return worst, changed


res = dict(base=BASE, real={}, sim_fast={}, sim_exact={})
L = ["# Sensitivity of the results to the numerical thresholds (N10; signet 1.2.0)", "",
     f"Defaults: T_zero = {BASE['T_zero']:g}, slack = {BASE['slack']:g}, implicit = {BASE['implicit']:g} (signet.core.TOL). "
     "Each threshold was multiplied by 0.01 and by 100, one at a time.", "",
     "| Setting | Real data: largest change in a set limit (log scale) | Real data: largest change in a p-value | Real data: sets changing between empty, bounded and unbounded | Simulation, archived rank: largest change in a limit | coverage or claim outcomes changed (of 300) | Simulation, exact rank: largest change in a limit | outcomes changed (of 300) |",
     "|---|---|---|---|---|---|---|---|"]
ref_real = ref_fast = ref_exact = None
for name, over in SETTINGS:
    with_tol(over)
    r, f, e = real_all(), sim_all(False), sim_all(True)
    if name == "default":
        ref_real, ref_fast, ref_exact = r, f, e
    dl, dc1 = maxdiff(ref_real, r, ("lo", "hi")); dp, dc2 = maxdiff(ref_real, r, ("p_benefit", "p_harm"))
    _, dc3 = maxdiff(ref_real, r, ("empty",))
    fl, fc = maxdiff(ref_fast, f, ("lo", "hi")); _, fo = maxdiff(ref_fast, f, ("cover", "benefit", "harm", "empty"))
    el, ec = maxdiff(ref_exact, e, ("lo", "hi")); _, eo = maxdiff(ref_exact, e, ("cover", "benefit", "harm", "empty"))
    res["real"][name] = dict(max_limit_change=dl, max_p_change=dp, state_changes=dc1 + dc2 + dc3)
    res["sim_fast"][name] = dict(max_limit_change=fl, state_changes=fc, outcome_changes=fo)
    res["sim_exact"][name] = dict(max_limit_change=el, state_changes=ec, outcome_changes=eo)
    L.append(f"| {name} | {dl:.1e} | {dp:.1e} | {dc1 + dc2 + dc3} | {fl:.1e} | {fo} | {el:.1e} | {eo} |")
with_tol({})
res["n_real_analyses"] = len(ref_real)
# exact versus archived rank at the default thresholds
dl_, _ = maxdiff(ref_fast, ref_exact, ("lo", "hi")); _, do_ = maxdiff(ref_fast, ref_exact, ("cover", "benefit", "harm", "empty"))
res["fast_vs_exact_default"] = dict(max_limit_change=dl_, outcome_changes=do_)
L += ["", f"Real-data analyses: {len(ref_real)} sets (every rung of the four questions, frozen and post hoc encodings of the method-agreement case). "
      f"Simulation: 300 literatures (20 per scenario). At the default thresholds, exact and archived ranks differed by at most {dl_:.1e} in a limit and in {do_} outcomes.", ""]

# (c) scale
L += ["## Rescaling the worked example", "",
      "Estimates and standard errors of the sigmoidoscopy example (rung 2) were multiplied by c. A method without absolute settings would return limits",
      "multiplied by c and unchanged p-values. The bracketing step (0.05), the bisection tolerance (1e-4) and the slack thresholds are absolute.", "",
      "| c | Typical standard error | Lower limit / c | Upper limit / c | p for benefit | Largest relative change in a limit |", "|---|---|---|---|---|---|"]
D = RM.D
scale = {}
ref = None
for c in (1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3):
    d, p = D["B_distal"], D["B_proximal"]
    M = Model(G=1, field_names=["distal"])
    M.add_row("B_distal", [1], ["H"]); M.add_row("B_proximal", [0], ["H"])
    M.set_data([d["y"] * c, p["y"] * c], [d["se"] * c, p["se"] * c])
    M.sign("H_toward_lower_risk", "H", -1, rung=1); M.share("H_shared", "H", rows=["B_distal", "B_proximal"], rung=2)
    inc = ["H_toward_lower_risk", "H_shared"]
    try:
        s = signet_set(M, inc, [1], 0.05)
        pb = claim_pvalue(M, inc, [1], "benefit", empty=bool(s["empty"]), span=3.0 * c)["p"]
        rec = dict(empty=bool(s["empty"]), lo=None if s["empty"] else s["lo"] / c, hi=None if s["empty"] else s["hi"] / c, p=pb)
    except Exception as ex:                                  # noqa: BLE001
        rec = dict(error=str(ex)[:80])
    scale[str(c)] = rec
    if c == 1.0:
        ref = rec
for c, rec in scale.items():
    if "error" in rec or rec.get("empty"):
        L.append(f"| {float(c):g} | {D['B_distal']['se'] * float(c):.2g} | {'error' if 'error' in rec else 'empty'} | | | |"); continue
    rel = max(abs(rec["lo"] - ref["lo"]) / abs(ref["lo"]), abs(rec["hi"] - ref["hi"]) / abs(ref["hi"])) if np.isfinite(rec["hi"]) and np.isfinite(rec["lo"]) else float("nan")
    rec["rel_change"] = rel
    L.append(f"| {float(c):g} | {D['B_distal']['se'] * float(c):.2g} | {rec['lo']:.5f} | {rec['hi']:.5f} | {rec['p']:.5f} | {rel:.1e} |")
res["scale"] = scale
json.dump(res, open(os.path.join(HERE, "tolerance_sensitivity.json"), "w"), indent=1, default=lambda o: None if (isinstance(o, float) and not np.isfinite(o)) else float(o))
open(os.path.join(HERE, "tolerance_sensitivity.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
