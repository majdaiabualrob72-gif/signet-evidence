"""Model builders of the real-data series (protocol encodings and the post hoc encoding of the method-agreement case).
Imported by realdata_audit.py and tolerance_sensitivity.py."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "realdata"))
import run_realdata as RD                 # noqa: E402
from signet.core import Model             # noqa: E402

D = RD.load()


def model_A(f_cvd=1.0, f_cancer=1.0):
    a, c = D["A_cvd"], D["A_cancer"]
    M = Model(G=2, field_names=["cvd", "cancer"])
    M.add_row("A_cvd", [1, 0], ["H"]); M.add_row("A_cancer", [0, 1], ["H", "D"])
    M.set_data([a["y"], c["y"]], [a["se"] * f_cvd, c["se"] * f_cancer])
    M.sign("H_toward_lower_risk", "H", -1, rung=1); M.sign("D_toward_higher_risk", "D", +1, rung=1)
    M.fieldsign("control_signed", 1, +1, rung=2); M.share("H_shared", "H", rows=["A_cvd", "A_cancer"], rung=2)
    M.linear("control_exact", {"cancer": 1.0}, "eq", 0.0, rung=3); M.zero("D_absent", "D", rung=3)
    R1 = ["H_toward_lower_risk", "D_toward_higher_risk"]; R2 = R1 + ["control_signed", "H_shared"]
    return M, dict(R1=R1, R2=R2, R3=R2 + ["control_exact", "D_absent"])


def model_B():
    d, p = D["B_distal"], D["B_proximal"]
    M = Model(G=1, field_names=["distal"])
    M.add_row("B_distal", [1], ["H"]); M.add_row("B_proximal", [0], ["H"])
    M.set_data([d["y"], p["y"]], [d["se"], p["se"]])
    M.sign("H_toward_lower_risk", "H", -1, rung=1); M.share("H_shared", "H", rows=["B_distal", "B_proximal"], rung=2)
    return M, dict(R1=["H_toward_lower_risk"], R2=["H_toward_lower_risk", "H_shared"])


def model_C():
    n, h = D["C_nhs"], D["C_hpfs"]
    M = Model(G=1, field_names=["chd"])
    M.add_row("C_nhs", [1], ["H"]); M.add_row("C_hpfs", [1], ["H"])
    M.set_data([n["y"], h["y"]], [n["se"], h["se"]])
    M.sign("H_toward_lower_risk", "H", -1, rung=1)
    return M, dict(R1=["H_toward_lower_risk"])


def model_D(trial_named):
    r, m, v = D["D_rct"], D["D_mr"], D["D_mvr"]
    M = Model(G=1, field_names=["chd"])
    M.add_row("D_rct", [1], ["A"] if trial_named else []); M.add_row("D_mr", [1], ["P"]); M.add_row("D_mvr", [1], ["H", "L"])
    M.set_data([r["y"], m["y"], v["y"]], [r["se"], m["se"], v["se"]])
    M.sign("H_toward_lower_risk", "H", -1, rung=1)
    M.zero("P_absent", "P", rung=2); M.zero("L_absent", "L", rung=2)
    R1 = ["H_toward_lower_risk"]; R2 = R1 + ["P_absent", "L_absent"]
    R = dict(R1=R1, R2=R2)
    if trial_named:
        M.zero("trial_valid", "A", rung=3); R["R3"] = R2 + ["trial_valid"]
    return M, R




def ANALYSES():
    """Every real-data analysis: (label, (model, rungs), target coefficients)."""
    return [("hormone therapy, cardiovascular disease", model_A(), [1, 0]),
            ("hormone therapy, total cancer", model_A(), [0, 1]),
            ("sigmoidoscopy, distal cancer death", model_B(), [1]),
            ("vitamin E, coronary disease", model_C(), [1]),
            ("beta-carotene, frozen encoding", model_D(False), [1]),
            ("beta-carotene, trial validity named", model_D(True), [1])]
