"""Specification layer: validation errors are caught and the CLI reproduces the real-data analysis."""
import json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
from signet.spec import validate, run, SpecError

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.join(HERE, "..", "examples")
ok = True

def check(name, cond, info=""):
    global ok
    print(("PASS " if cond else "FAIL ") + name, info)
    ok &= bool(cond)

bad = {"cells": [{"id": "x", "estimate": -1}], "assumptions": [{"name": "a", "type": "sign", "process": "Q"}]}
try:
    validate(bad); check("invalid spec rejected", False)
except SpecError as e:
    check("invalid spec rejected", "process Q" in str(e) and "positive" in str(e), str(e)[:80])

r = run(json.load(open(os.path.join(EX, "hormone_therapy_cvd.json"))))
t = r["targets"]["cvd"]
check("rung 1 cannot certify benefit", t["rung1"]["capacity"]["above"] is False)
check("rung 3 identified at estimates equals log(0.57/0.83)",
      abs(t["rung3"]["identified_at_estimates"]["hi"] - np.log(0.57 / 0.83)) < 1e-6,
      t["rung3"]["identified_at_estimates"]["hi"])
check("k = 1 removes the rung-3 benefit claim", t["rung3"]["tolerance_k1"]["hi"] == "inf")
r2 = run(json.load(open(os.path.join(EX, "sigmoidoscopy_distal.json"))))
s = r2["targets"]["distal"]["rung2"]
check("anatomical control certifies benefit at rung 2", s["set"]["hi"] < 0, s["set"]["ratio_scale"])
check("benefit p-value below 0.05", s["claims"]["benefit"]["p"] < 0.05, s["claims"]["benefit"]["p"])
print("ALL PASS" if ok else "SOME FAILED")
sys.exit(0 if ok else 1)
