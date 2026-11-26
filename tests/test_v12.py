"""Tests for the v1.2.0 changes (response to the forensic audit of 22 November 2026).

Run with:  python tests/test_v12.py     (or: python -m pytest -q tests/test_v12.py)
Findings addressed: D5 (empty sets), M3 (capacity of the default), N1 (kernel condition of the corollary),
N10 (named thresholds), D4 (non-nested rungs), M1 and D2 (comparator specifications), D3 (group readings).
"""
import json
import os
import sys

import numpy as np
from scipy import stats as sst

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "simulation"))
sys.path.insert(0, os.path.join(ROOT, "experiments", "audit_v12"))
import signet.core as C                                                                   # noqa: E402
from signet.core import (Model, capacity, signet_set, claim_pvalue, default_capacity,     # noqa: E402
                         default_set, tolerant_set, cert_bounds)
from signet.spec import run, report_md, scale_warnings                                    # noqa: E402

Z = float(sst.norm.ppf(0.975))


def _se(lo, hi):
    return (np.log(hi) - np.log(lo)) / (2 * Z)


def _sigmoidoscopy():
    M = Model(G=1, field_names=["distal"])
    M.add_row("distal", [1], ["H"]); M.add_row("proximal", [0], ["H"])
    M.set_data([np.log(0.41), np.log(0.96)], [_se(0.25, 0.69), _se(0.61, 1.50)])
    M.sign("H-", "H", -1, rung=1); M.share("shH", "H", rows=["distal", "proximal"], rung=2)
    return M, ["H-"], ["H-", "shH"]


def _conflict():
    """Two valid cells that measure the same effect and disagree by 40 standard errors."""
    M = Model(G=1)
    M.add_row("a", [1], ["A"]); M.add_row("b", [1], ["B"])
    M.set_data([-1.0, 1.0], [0.05, 0.05])
    M.zero("vA", "A", rung=3); M.zero("vB", "B", rung=3)
    return M, ["vA", "vB"]


def test_empty_set_reports_no_directional_p_value():
    M, inc = _conflict()
    s = signet_set(M, inc, [1.0])
    assert s["empty"]
    for claim in ("benefit", "harm"):
        c = claim_pvalue(M, inc, [1.0], claim)
        assert c["set_empty"] is True and c["p"] is None and c["p_upper"] is None
        # the uncorrected formula would have been small for both directions at once
        assert min(1.0, 4 * c["p_cert"], 2 * c["p_cc"]) < 1e-6
        assert "empty" in c["note"]
        # the caller may pass the emptiness it has already computed
        assert claim_pvalue(M, inc, [1.0], claim, empty=True)["p"] is None


def test_non_empty_set_p_value_unchanged():
    M, R1, R2 = _sigmoidoscopy()
    c = claim_pvalue(M, R2, [1.0], "benefit")
    assert c["set_empty"] is False
    assert abs(c["p"] - min(1.0, 4 * c["p_cert"], 2 * c["p_cc"])) < 1e-15
    assert abs(c["p"] - 0.0279112) < 1e-6
    assert claim_pvalue(M, R2, [1.0], "benefit", empty=False)["p"] == c["p"]


def test_spec_run_marks_empty_set_claims_as_not_reported():
    spec = {"name": "conflict", "field": ["theta"],
            "cells": [{"id": "a", "estimate": float(np.exp(-1.0)), "se": 0.05, "processes": ["A"]},
                      {"id": "b", "estimate": float(np.exp(1.0)), "se": 0.05, "processes": ["B"]}],
            "assumptions": [{"name": "vA", "type": "zero", "process": "A", "rung": 3, "source": "test"},
                            {"name": "vB", "type": "zero", "process": "B", "rung": 3, "source": "test"}],
            "analysis": {"rungs": [3], "tolerance": [0, 1]}}
    out = run(spec)
    rec = out["targets"]["theta"]["rung3"]
    assert rec["set"] is None
    assert rec["claims"]["benefit"]["p"] is None and rec["claims"]["harm"]["p"] is None
    assert "not reported (empty set)" in report_md(out)


def test_default_capacity_single_effect_cell_is_unbounded_whatever_the_data():
    M, R1, R2 = _sigmoidoscopy()
    assert capacity(M, R2, [1.0]) == dict(above=True, below=True)          # the structural rung can bound the effect
    dc = default_capacity(M, R2, [1.0])
    assert dc["above"] is False and dc["below"] is False
    assert "without cell distal" in dc["limiting"]["above"]
    d = default_set(M, R2, [1.0])
    assert not np.isfinite(d["lo"]) and not np.isfinite(d["hi"])
    # the data cannot change this: any estimates give an unbounded default
    M.set_data([np.log(0.05), np.log(1.0)], [0.01, 0.01])
    d = default_set(M, R2, [1.0])
    assert not np.isfinite(d["lo"]) and not np.isfinite(d["hi"])


def test_default_capacity_with_three_routes_can_be_bounded():
    M = Model(G=1)
    for nm, p in (("t1", "A1"), ("t2", "A2"), ("t3", "A3")):
        M.add_row(nm, [1], [p])
    M.set_data([-0.3, -0.25, -0.35], [0.1, 0.1, 0.1])
    inc = []
    for p in ("A1", "A2", "A3"):
        M.zero("v" + p, p, rung=3); inc.append("v" + p)
    dc = default_capacity(M, inc, [1.0])
    assert dc["above"] and dc["below"] and dc["limiting"] == dict(above=[], below=[])
    d = default_set(M, inc, [1.0])
    assert np.isfinite(d["lo"]) and np.isfinite(d["hi"])


def test_corollary_needs_non_negative_kernels():
    """N1. One-signed cells cannot certify benefit when every kernel is non-negative; a cell with a negative
    kernel and its own non-positive bias bounds the effect from above, so the condition is needed."""
    M = Model(G=1)
    M.add_row("c1", [1.0], ["H1"]); M.add_row("c2", [1.0], ["H2"])
    M.set_data([-0.4, -0.3], [0.1, 0.1])
    M.sign("H1-", "H1", -1); M.sign("H2-", "H2", -1)
    cap = capacity(M, ["H1-", "H2-"], [1.0])
    assert cap["above"] is False and cap["below"] is True
    N = Model(G=1)
    N.add_row("c1", [1.0], ["H1"]); N.add_row("c2", [-1.0], ["H2"])
    N.set_data([-0.4, 0.3], [0.1, 0.1])
    N.sign("H1-", "H1", -1); N.sign("H2-", "H2", -1)
    assert capacity(N, ["H1-", "H2-"], [1.0])["above"] is True


def test_named_thresholds_default_values_and_insensitivity():
    assert C.TOL == dict(T_zero=1e-9, slack=1e-4, implicit=1e-7)
    M, R1, R2 = _sigmoidoscopy()
    base = signet_set(M, R2, [1.0])
    try:
        for k, f in (("slack", 100.0), ("slack", 0.01), ("T_zero", 100.0), ("implicit", 100.0)):
            C.TOL[k] = dict(T_zero=1e-9, slack=1e-4, implicit=1e-7)[k] * f
            s = signet_set(M, R2, [1.0])
            assert abs(s["lo"] - base["lo"]) < 1e-9 and abs(s["hi"] - base["hi"]) < 1e-9
            C.TOL.update(dict(T_zero=1e-9, slack=1e-4, implicit=1e-7))
    finally:
        C.TOL.update(dict(T_zero=1e-9, slack=1e-4, implicit=1e-7))


def test_rungs_need_not_be_nested_and_the_reason():
    """D4. The rule prefers the exact link at rung 2, whose certificate has the larger variance."""
    M, R1, R2 = _sigmoidoscopy()
    s1, s2 = signet_set(M, R1, [1.0]), signet_set(M, R2, [1.0])
    assert abs(np.exp(s1["lo"]) - 0.229) < 1e-3 and abs(np.exp(s2["lo"]) - 0.203) < 1e-3
    assert s2["lo"] < s1["lo"]
    c1, c2 = cert_bounds(M, R1, [1.0], 0.0125), cert_bounds(M, R2, [1.0], 0.0125)
    assert c2["cert_lo"]["sd"] > c1["cert_lo"]["sd"]
    w2 = c2["cert_lo"]["weights"]
    assert abs(w2["distal"] - 1) < 1e-6 and abs(w2["proximal"] + 1) < 1e-6


def test_scale_warning():
    spec = {"cells": [{"id": "x", "estimate": 1.0, "lower": 0.98, "upper": 1.01}, {"id": "y", "estimate": 0.8, "lower": 0.6, "upper": 1.07}]}
    w = scale_warnings(spec)
    assert len(w) == 1 and "cell x" in w[0]


def test_report_prints_structural_rung_beside_default():
    spec = json.load(open(os.path.join(ROOT, "examples", "sigmoidoscopy_distal.json")))
    txt = report_md(run(spec))
    assert "Beside it, rung2 without tolerance" in txt and "0.203 to 0.928" in txt
    assert "Capacity of the default: upper limit unbounded whatever the data" in txt
    assert "One cell measures this target" in txt


def test_random_effects_interval_variants():
    """M1. Vitamin E, two cohorts: Wald, Hartung-Knapp and modified Hartung-Knapp intervals."""
    import sim_audit as S
    from decode_core import re_meta
    y = [np.log(0.59), np.log(0.63)]; s = [_se(0.38, 0.91), _se(0.47, 0.84)]
    v = S.re_variants(y, s)
    r = lambda t: (round(float(np.exp(t[0])), 2), round(float(np.exp(t[1])), 2))
    assert r(v["wald"]) == (0.48, 0.79) and r(v["hk"]) == (0.42, 0.91) and r(v["mhk"]) == (0.13, 2.96)
    assert abs(v["hk"][2] - 0.040) < 1e-3 and abs(v["mhk"][2] - 0.16) < 5e-3 and v["wald"][2] < 1e-4
    a = re_meta(y, s)                                   # the archived simulation comparator is the modified interval
    assert abs(a["lo"] - v["mhk"][0]) < 1e-12 and abs(a["hi"] - v["mhk"][1]) < 1e-12
    one = S.re_variants([0.1], [0.2])
    assert one["wald"] == one["hk"] == one["mhk"]


def test_ba_s_shared_closed_form():
    """D2. Two cells sharing one bias term: mean = fixed effect + 0.15, variance = fixed-effect variance + 0.15^2."""
    import benchmark_v2 as B2
    import sim_audit as S
    y = np.array([np.log(0.59), np.log(0.63)]); s = np.array([_se(0.38, 0.91), _se(0.47, 0.84)])
    cl = [dict(name="NU-U", K=1.0, procs=["H"], y=y[0], s=s[0]), dict(name="NU-L", K=1.0, procs=["H"], y=y[1], s=s[1])]
    lo, hi = S.ba_s(cl, True)
    w = 1 / s ** 2; fe = np.sum(w * y) / np.sum(w); vfe = 1 / np.sum(w)
    m = fe + 0.15; se = np.sqrt(vfe + 0.15 ** 2)
    assert abs(lo - (m - Z * se)) < 1e-12 and abs(hi - (m + Z * se)) < 1e-12
    assert (round(float(np.exp(m)), 2), round(float(np.exp(lo)), 2), round(float(np.exp(hi)), 2)) == (0.72, 0.49, 1.05)
    lo_i, hi_i = S.ba_s(cl, False)
    assert (round(float(np.exp(lo_i)), 2), round(float(np.exp(hi_i)), 2)) == (0.52, 0.99)
    a_lo, a_hi = B2.bias_adjusted(cl, True)             # the archived function is the independent specification
    assert abs(a_lo - lo_i) < 1e-12 and abs(a_hi - hi_i) < 1e-12


def test_group_readings_are_nested():
    """D3. Finer groups give a smaller union, coarser groups a larger one."""
    import benchmark as B
    import benchmark_v2 as B2
    import sim_audit as S
    for r in (0, 1, 2):
        sets = {}
        for reading in ("row", "statement", "process"):
            rng = np.random.default_rng(S.seed_of("contradictory evidence", r))
            theta, b, present, studies = B2.literature_v2("contradictory evidence", rng)
            m, R = B.signet_model(present, studies)
            n_rows = sum(len(g["cons"]) for g in m.groups)
            names = R["R3"] if reading == "statement" else S.regroup(m, R["R3"], reading)
            assert sum(len(g["cons"]) for g in m.groups) == n_rows            # no constraint is lost
            sets[reading] = tolerant_set(m, names, [1.0], k=1, **S.KW)
        lo = {k: (np.inf if v["empty"] else v["lo"]) for k, v in sets.items()}
        hi = {k: (-np.inf if v["empty"] else v["hi"]) for k, v in sets.items()}
        assert lo["process"] <= lo["statement"] + 2e-3 and lo["statement"] <= lo["row"] + 2e-3
        assert hi["process"] >= hi["statement"] - 2e-3 and hi["statement"] >= hi["row"] - 2e-3


if __name__ == "__main__":
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print("PASS", name)
            except Exception as e:      # noqa: BLE001
                fails += 1
                print("FAIL", name, repr(e))
    print("ALL PASS" if not fails else f"{fails} FAILED")
    sys.exit(1 if fails else 0)
