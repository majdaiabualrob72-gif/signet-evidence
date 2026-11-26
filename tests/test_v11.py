"""Tests for the v1.1.0 changes (response to the GENESIS audit of November 2026).

Run with:  python -m pytest -q tests/test_v11.py
"""
import os
import sys

import numpy as np
from scipy import stats as sst

try:
    import pytest
except ImportError:                      # minimal stand-ins, so that `python tests/test_v11.py` runs without pytest
    class _Approx:
        __array_ufunc__ = None           # numpy defers the comparison to this object

        def __init__(self, v, rel=1e-6, abs=1e-12):
            self.v, self.rel, self.abs = float(v), rel, abs

        def __eq__(self, o):
            return abs(float(o) - self.v) <= max(self.rel * abs(self.v), self.abs)

    class _Raises:
        def __init__(self, exc):
            self.exc = exc

        def __enter__(self):
            return self

        def __exit__(self, et, ev, tb):
            if et is None:
                raise AssertionError(f"{self.exc.__name__} not raised")
            return issubclass(et, self.exc)

    class pytest:                        # noqa: N801
        approx = _Approx
        raises = _Raises

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from signet.core import (Model, Fitter, conf_set, signet_set, certificate, cert_bounds,  # noqa: E402
                         check_certificate, certificate_range, claim_pvalue, INF)
from signet.spec import build, validate, SpecError  # noqa: E402


def three_cell(scale_order=1.0, as_linear=True):
    """Trial, confounded cohort, negative-control cohort; an ordering assumption on the
    confounding bias, optionally written with a scaled coefficient row."""
    M = Model(1)
    M.add_row("RCT", [1.0], ["A"])
    M.add_row("COH", [1.0], ["H"])
    M.add_row("NCO", [0.0], ["H2"])
    M.set_data([-0.20, -0.35, -0.10], [0.10, 0.06, 0.05])
    M.zero("A0", "A", rung=1)
    M.sign("Hneg", "H", -1, rung=1)
    M.sign("H2neg", "H2", -1, rung=1)
    # the control's bias is at least as negative as the cohort's: b_NCO - b_COH <= 0
    M.linear("ordH", {"H2@NCO": scale_order, "H@COH": -scale_order}, "le", 0.0, rung=2)
    return M


INC = ["A0", "Hneg", "H2neg", "ordH"]


def two_route(scale=1.0):
    """Two routes to an upper bound for f: the sign assumption (f <= mu_COH) or the ordering
    assumption against a negative control (f <= mu_COH - mu_NCO). Under the v1.0.0 cost the two
    tie when the ordering row has unit coefficients, and scaling that row breaks the tie."""
    M = Model(1)
    M.add_row("COH", [1.0], ["H"])
    M.add_row("NCO", [0.0], ["H2"])
    M.set_data([-0.30, -0.10], [0.10, 0.05])
    M.sign("Hpos", "H", +1, rung=1)
    M.linear("ord", {"H2@NCO": scale, "H@COH": -scale}, "le", 0.0, rung=1)
    return M


def test_row_scaling_invariance_of_certificates_and_sets():
    inc = ["Hpos", "ord"]
    M1, M7 = two_route(1.0), two_route(7.0)
    for side in ("hi", "lo"):
        c1, c7 = certificate(M1, inc, [1.0], side), certificate(M7, inc, [1.0], side)
        assert (c1 is None) == (c7 is None)
        if c1 is not None:
            assert np.allclose(c1["lam"], c7["lam"], atol=1e-6)
    s1 = signet_set(M1, inc, [1.0], 0.05)
    s7 = signet_set(M7, inc, [1.0], 0.05)
    for k in ("lo", "hi"):
        assert (not np.isfinite(s1[k]) and not np.isfinite(s7[k])) or abs(s1[k] - s7[k]) < 1e-3


def test_without_normalisation_scaling_changes_the_choice():
    # documents why v1.1.0 normalises: the v1.0.0 first-stage cost was not scale invariant
    inc = ["Hpos", "ord"]
    M1, M7 = two_route(1.0), two_route(7.0)
    M1.normalise_rows = M7.normalise_rows = False
    c1 = certificate(M1, inc, [1.0], "hi")
    c7 = certificate(M7, inc, [1.0], "hi")
    assert not np.allclose(c1["lam"], c7["lam"], atol=1e-3)


def test_outer_set_contains_hull_and_has_exact_end_points():
    M = three_cell()
    hull = conf_set(M, INC, [1.0], 0.05)
    out = conf_set(M, INC, [1.0], 0.05, inversion="outer")
    assert out["lo"] <= hull["lo"] + 1e-6 and hull["hi"] <= out["hi"] + 1e-6
    F = Fitter(M, INC)
    cmax = sst.chi2.ppf(0.95, M.n)
    c = M.cvec([1.0])
    for t in (out["lo"], out["hi"]):
        if np.isfinite(t):
            assert abs(F.fit(c, t)["T"] - cmax) < 1e-2
    assert hull["bracketed"]["lo"] in (True, None) and hull["bracketed"]["hi"] in (True, None)


def test_unbracketed_limit_is_reported_as_unbounded_with_a_warning():
    M = three_cell()
    s = conf_set(M, INC, [1.0], 0.05, max_span=0.2)
    assert s["warnings"]
    assert not np.isfinite(s["lo"]) or not np.isfinite(s["hi"])


def test_dense_check_does_not_shrink_the_hull():
    M = three_cell()
    a = conf_set(M, INC, [1.0], 0.05, outer_scan=0)
    b = conf_set(M, INC, [1.0], 0.05, outer_scan=32)
    assert b["lo"] <= a["lo"] + 1e-9 and a["hi"] <= b["hi"] + 1e-9


def test_validated_p_value_bound_dominates_grid_value():
    M = three_cell()
    for claim in ("benefit", "harm"):
        r = claim_pvalue(M, INC, [1.0], claim)
        assert r["p_cc_upper"] >= r["p_cc"] - 1e-12
        assert r["p_upper"] >= r["p"] - 1e-12


def test_frozen_certificate_is_valid_and_uses_data_variance():
    M = three_cell()
    cb = cert_bounds(M, INC, [1.0], 0.0125, cert_V=np.eye(M.n))
    for side in ("hi", "lo"):
        ct = cb["cert_" + side]
        if ct is None:
            continue
        lam = ct["lam"] if side == "hi" else -ct["lam"]
        assert check_certificate(M, INC, [1.0], lam, side)["valid"]
        assert ct["sd"] == pytest.approx(float(np.sqrt(ct["lam"] @ M.V @ ct["lam"])))


def test_check_certificate_rejects_a_non_certificate():
    M = three_cell()
    assert not check_certificate(M, INC, [1.0], np.array([0.0, 0.0, 1.0]), "hi")["valid"]


def test_non_homogeneous_certificate_offset():
    # theta = f, one cell y = f + b with -0.1 <= b <= 0.3: f = mu - b <= mu + 0.1, f >= mu - 0.3
    M = Model(1)
    M.add_row("OBS", [1.0], ["B"])
    M.set_data([0.0], [0.1])
    M.bound("Bb", "B", -0.1, 0.3, rung=4)
    hi = check_certificate(M, ["Bb"], [1.0], np.array([1.0]), "hi")
    lo = check_certificate(M, ["Bb"], [1.0], np.array([1.0]), "lo")
    assert hi["valid"] and hi["offset"] == pytest.approx(0.1)
    assert lo["valid"] and lo["offset"] == pytest.approx(0.3)


def test_certificate_uniqueness():
    M = Model(G=2, field_names=["distal", "proximal"])
    M.add_row("d", [1, 0], ["H"]); M.add_row("p", [0, 1], ["H"])
    M.set_data([-0.5, -0.05], [0.08, 0.07])
    M.sign("Hneg", "H", -1, rung=1); M.share("Hsh", "H", rung=2)
    r = certificate_range(M, ["Hneg", "Hsh"], [1, -1], "hi")
    assert r["unique"]
    assert certificate_range(M, ["Hneg", "Hsh"], [1, 0], "hi") is None


def test_approximate_rank_never_below_exact_rank():
    rng = np.random.default_rng(11)
    M = three_cell()
    c = M.cvec([1.0])
    for _ in range(25):
        y = np.array([-0.2, -0.35, -0.1]) + rng.normal(0, 0.15, 3)
        Fe = Fitter(M, INC, y=y, exact_rank=True)
        Fa = Fitter(M, INC, y=y, exact_rank=False)
        for t in np.linspace(-1.0, 0.6, 9):
            fr = Fe.fit(c, t)
            if fr["ok"] and fr["T"] > 1e-9:
                assert Fa.rank(Fa.fit(c, t)) >= Fe.rank(fr)


SPEC = {
    "name": "correlated cells",
    "field": ["theta"],
    "cells": [
        {"id": "a", "estimate": 0.8, "lower": 0.6, "upper": 1.07, "processes": ["H"], "source": "x"},
        {"id": "b", "estimate": 0.9, "lower": 0.7, "upper": 1.16, "processes": ["H"], "source": "x"},
        {"id": "c", "estimate": 1.0, "se": 0.1, "scale": "log", "kernel": [0.0], "processes": ["N"], "source": "x"},
    ],
    "assumptions": [{"name": "Hneg", "type": "sign", "process": "H", "direction": "negative", "rung": 1, "source": "x"}],
    "correlations": [{"cells": ["a", "b"], "rho": 0.5}],
}


def test_spec_correlations_enter_the_covariance():
    M, _ = build(SPEC)
    s = np.sqrt(np.diag(M.V))
    assert M.V[0, 1] == pytest.approx(0.5 * s[0] * s[1])
    assert M.V[0, 2] == 0.0


def test_spec_rejects_bad_correlations():
    bad = dict(SPEC, correlations=[{"cells": ["a", "b"], "rho": 1.0}])
    with pytest.raises(SpecError):
        validate(bad)
    bad = dict(SPEC, correlations=[{"cells": ["a", "z"], "rho": 0.2}])
    with pytest.raises(SpecError):
        validate(bad)
    bad = dict(SPEC, correlations=[{"cells": ["a", "b"], "rho": 0.9}, {"cells": ["a", "c"], "rho": 0.9},
                                   {"cells": ["b", "c"], "rho": -0.9}])
    with pytest.raises(SpecError):
        validate(bad)


def test_spec_certificate_se_freezes_the_choice():
    from signet.spec import run
    spec = dict(SPEC, certificate_se={"a": 0.2, "b": 0.2, "c": 0.1})
    M, rungs = build(spec)
    assert np.allclose(np.diag(M.V) > 0, True)
    assert M.cert_V[0, 0] == pytest.approx(0.04) and M.cert_V[0, 1] == pytest.approx(0.5 * 0.04)
    out = run(spec)
    assert out["certificate_selection"].startswith("frozen")
    bad = dict(SPEC, certificate_se={"a": 0.2, "b": 0.2})
    with pytest.raises(SpecError):
        validate(bad)


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
