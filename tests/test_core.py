"""Closed-form and Monte Carlo checks of signet_core."""
import numpy as np
from scipy import stats as sst
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from signet.core import *

z = sst.norm.ppf(0.975)
ok = True


def check(name, cond, detail=""):
    global ok
    ok &= bool(cond)
    print(("PASS " if cond else "FAIL ") + name + ("  " + detail if detail else ""))


# 1. sharing with a negative control point-identifies the effect
m = Model(1)
m.add_row("lagged", [1.0], ["H"]); m.add_row("nco", [0.0], ["H"])
m.share("shareH", "H", rung=2)
y = np.array([-0.40, -0.12]); s = np.array([0.05, 0.04]); m.set_data(y, s)
cs = conf_set(m, ["shareH"], [1.0])
lo, hi = y[0] - y[1] - z * np.hypot(*s), y[0] - y[1] + z * np.hypot(*s)
check("point identification by a shared negative control", abs(cs["lo"] - lo) < 2e-4 and abs(cs["hi"] - hi) < 2e-4,
      f"{cs['lo']:.4f},{cs['hi']:.4f} vs {lo:.4f},{hi:.4f}")

# 2. bracketing by two designs biased in opposite directions
m = Model(1)
m.add_row("harmward", [1.0], ["R"]); m.add_row("benefitward", [1.0], ["H"])
m.sign("R+", "R", +1); m.sign("H-", "H", -1)
y = np.array([-0.10, -0.45]); s = np.array([0.06, 0.05]); m.set_data(y, s)
cs = conf_set(m, ["R+", "H-"], [1.0])
check("bracket from opposite signs", abs(cs["hi"] - (y[0] + z * s[0])) < 2e-4 and abs(cs["lo"] - (y[1] - z * s[1])) < 2e-4,
      f"{cs['lo']:.4f},{cs['hi']:.4f}")
idf = identified(m, ["R+", "H-"], [1.0], y)
check("certificate of the upper bound", idf["cert_hi"]["weights"] == {"harmward": 1.0} and idf["cert_hi"]["assumptions"] == ["R+"],
      str(idf["cert_hi"]))

# 3. capacity: one design biased towards harm bounds the effect from above only
cap = capacity(m, ["R+"], [1.0])
check("cone criterion, one-sided design", cap == {"above": True, "below": False}, str(cap))
cap = capacity(m, [], [1.0])
check("cone criterion, no sign knowledge", cap == {"above": False, "below": False}, str(cap))

# 4. tolerance: a bound that rests on one assumption is lost at k = 1;
#    a bound with two independent routes survives
ts = tolerant_set(m, ["R+", "H-"], [1.0], k=1)
check("tolerance dial removes one-route bounds", ts["hi"] == INF and ts["lo"] == -INF, str(ts))
m2 = Model(1)
m2.add_row("reg1", [1.0], ["R1"]); m2.add_row("reg2", [1.0], ["R2"]); m2.add_row("nu", [1.0], ["H"])
m2.sign("R1+", "R1", +1); m2.sign("R2+", "R2", +1); m2.sign("H-", "H", -1)
m2.set_data([-0.10, -0.12, -0.45], [0.06, 0.06, 0.05])
ts = tolerant_set(m2, ["R1+", "R2+", "H-"], [1.0], k=1)
check("tolerance dial keeps a two-route bound", np.isfinite(ts["hi"]) and ts["lo"] == -INF, str(ts))

# 5. intersection bounds: two upper-bounding designs, both unbiased, truth on the boundary
rng = np.random.default_rng(1)
m = Model(1)
m.add_row("a", [1.0], ["Ra"]); m.add_row("b", [1.0], ["Rb"])
m.sign("Ra+", "Ra", +1); m.sign("Rb+", "Rb", +1)
F = None; cover = []; prof = []
for r in range(4000):
    y = rng.normal(0, 1, 2)
    F = Fitter(m, ["Ra+", "Rb+"], y, np.eye(2))
    a, T, _ = F.accept(np.r_[1.0, 0, 0], 0.0)
    cover.append(a)
    Tmin = F.fit()["T"]
    prof.append(T - Tmin <= sst.chi2.ppf(0.95, 1))
check("conditional chi-squared covers at an intersection point", np.mean(cover) >= 0.945,
      f"CC {np.mean(cover):.4f} (theory 0.9625); profile chi2(1) {np.mean(prof):.4f} (theory 0.938)")

# 6. over-identified equality model: exact coverage
m = Model(1)
m.add_row("lag", [1.0], ["H"]); m.add_row("nco", [0.0], ["H"]); m.add_row("rct", [1.0], [])
m.share("shH", "H")
cover = []
for r in range(4000):
    y = np.array([-0.3 - 0.1, -0.1, -0.3]) + rng.normal(0, 1, 3) * np.array([0.05, 0.04, 0.10])
    F = Fitter(m, ["shH"], y, np.diag([0.05, 0.04, 0.10]) ** 2)
    cover.append(F.accept(np.r_[1.0, 0, 0], -0.3)[0])
check("exact size in an over-identified model", abs(np.mean(cover) - 0.95) < 0.012, f"{np.mean(cover):.4f}")

# 7. repairs: a single invalid route is isolated
m = Model(1)
m.add_row("r1", [1.0], ["B1"]); m.add_row("r2", [1.0], ["B2"]); m.add_row("r3", [1.0], ["B3"])
m.zero("v1", "B1", rung=3); m.zero("v2", "B2", rung=3); m.zero("v3", "B3", rung=3)
m.set_data([-0.30, -0.29, 0.20], [0.05, 0.05, 0.05])
mt = model_test(m, ["v1", "v2", "v3"])
rp = repairs(m, ["v1", "v2", "v3"], kmax=1)
check("falsification and minimal repair", mt["reject"] and rp == [("v3",)], f"T={mt['T']:.1f}, r={mt['r']}, repairs={rp}")

# 8. fast and exact rank agree on random problems
agree = 0; tot = 0
m = Model(1)
for nm, pr in (("u", ["H", "P"]), ("l", ["H"]), ("n", ["H"]), ("g", ["R"]), ("t", ["A"])):
    m.add_row(nm, [0.0 if nm == "n" else 1.0], pr)
m.sign("H-", "H", -1); m.sign("P-", "P", -1); m.sign("R+", "R", +1); m.share("shH", "H"); m.zero("vA", "A")
inc = ["H-", "P-", "R+", "shH", "vA"]
for r in range(300):
    y = rng.normal(-0.3, 0.2, 5)
    Fe = Fitter(m, inc, y, np.diag(np.full(5, 0.05 ** 2)), True)
    Ff = Fitter(m, inc, y, np.diag(np.full(5, 0.05 ** 2)), False)
    th = rng.normal(-0.3, 0.3)
    fe = Fe.fit(m.cvec([1.0]), th); ff = Ff.fit(m.cvec([1.0]), th)
    if fe["ok"] and fe["T"] > 1e-6:
        tot += 1
        re, rf = Fe.rank(fe), Ff.rank(ff)
        agree += (re == rf)
        assert rf >= re, (re, rf)
check("fast rank equals exact rank (fast never smaller)", agree == tot, f"{agree}/{tot}")
print("ALL PASS" if ok else "SOME FAILED")

# 9. certificate bounds: closed form in the shared negative-control design
m = Model(1)
m.add_row("lagged", [1.0], ["H"]); m.add_row("nco", [0.0], ["H"]); m.add_row("reg", [1.0], ["R"])
m.share("shareH", "H"); m.sign("H-", "H", -1); m.sign("R+", "R", +1)
y = np.array([-0.40, -0.12, -0.20]); s = np.array([0.05, 0.04, 0.08]); m.set_data(y, s)
cb = cert_bounds(m, ["shareH", "H-", "R+"], [1.0], a=0.025)
# an exact (equality-only) route exists, so both sides use lagged - nco and no sign
v1 = s[0] ** 2 + s[1] ** 2
check("certificate prefers the exact route (upper)", abs(cb["hi"] - (y[0] - y[1] + z * np.sqrt(v1))) < 1e-6
      and cb["cert_hi"]["slack_weight"] < 1e-9, f"{cb['hi']:.5f}; {cb['cert_hi']['weights']}")
check("certificate prefers the exact route (lower)", abs(cb["lo"] - (y[0] - y[1] - z * np.sqrt(v1))) < 1e-6,
      f"{cb['lo']:.5f}")
# without sharing, the upper side must rest on the registry sign and the lower on the healthy-user sign
cb2 = cert_bounds(m, ["H-", "R+"], [1.0], a=0.025)
check("sign-based certificates without sharing", abs(cb2["hi"] - (y[2] + z * s[2])) < 1e-6 and abs(cb2["lo"] - (y[0] - z * s[0])) < 1e-6
      and cb2["cert_hi"]["assumptions"] == ["R+"], f"{cb2['lo']:.4f}, {cb2['hi']:.4f}")
# Monte Carlo: certificate bound covers at its nominal one-sided level even when the registry bias is large
rng = np.random.default_rng(7)
mu = np.array([-0.3 - 0.1, -0.1, -0.3 + 0.25]); hits = 0
for r in range(4000):
    yy = mu + rng.normal(0, 1, 3) * s
    hits += cert_bounds(m, ["shareH", "H-", "R+"], [1.0], a=0.025, y=yy)["hi"] >= -0.3
check("certificate bound is valid whatever the bias size", hits / 4000 >= 0.975 - 0.006, f"{hits / 4000:.4f}")
print("ALL PASS" if ok else "SOME FAILED")
