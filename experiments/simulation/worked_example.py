"""
worked_example.py
Hypothetical drug X for disease N. The effect of X (log hazard ratio) is a
field over disease duration at initiation, in five bins. Fourteen evidence
cells, pooled from 41 studies, measure population mixtures of that field
plus signed bias processes. One assumption (surrogate validity) is false
by construction. Usage: python3 worked_example.py [main|design]
"""
import json, sys, time
import numpy as np
from scipy import stats as sst
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
sys.path.insert(0, _os.path.abspath(_os.path.join(_HERE, '..', '..')))
sys.path.insert(0, _HERE)
from signet.core import *

SEED = 20261117
BINS = ["0-3 y", "3-6 y", "6-9 y", "9-12 y", "12-15 y"]
KE = np.array([0.45, 0.35, 0.15, 0.05, 0.00])
KM = np.array([0.10, 0.25, 0.35, 0.20, 0.10])
KL = np.array([0.00, 0.05, 0.15, 0.35, 0.45])
F_TRUE = np.array([-0.30, -0.30, -0.15, 0.00, 0.00])
B_TRUE = dict(H=-0.15, P=-0.22, C=0.0, R=0.10, A=0.0, S=-0.25, Pl=0.0)
Z0 = np.zeros(5)

# name, mixture, processes, studies, study SE range, true bias overrides
CELLS = [
    ("Non-user cohort, unlagged, early", KE, ["H", "P"], 6, (0.08, 0.15), {}),
    ("Non-user cohort, lagged, early",   KE, ["H", "P"], 4, (0.09, 0.16), {"P": 0.0}),
    ("Non-user cohort, unlagged, mid",   KM, ["H", "P"], 5, (0.08, 0.15), {}),
    ("Non-user cohort, lagged, mid",     KM, ["H", "P"], 3, (0.09, 0.16), {"P": 0.0}),
    ("Non-user negative control",        Z0, ["H"],      4, (0.08, 0.14), {}),
    ("Active comparator, early",         KE, ["C"],      3, (0.10, 0.18), {}),
    ("Active comparator, mid",           KM, ["C"],      3, (0.10, 0.18), {}),
    ("Active-comparator negative control", Z0, ["C"],    3, (0.08, 0.15), {}),
    ("Registry (channelled), early",     KE, ["R"],      2, (0.10, 0.16), {}),
    ("Registry (channelled), mid",       KM, ["R"],      2, (0.10, 0.16), {}),
    ("Trial, mid duration",              KM, ["A"],      2, (0.12, 0.20), {}),
    ("Trial, late duration",             KL, ["A"],      3, (0.10, 0.16), {}),
    ("Biomarker chain, mid",             KM, ["S"],      1, (0.05, 0.05), {}),
    ("Mendelian randomisation",          KE, ["Pl"],     2, (0.10, 0.15), {}),
]

RUNG1 = ["H favours X", "P favours X", "R opposes X", "Lag removes P"]
RUNG2 = RUNG1 + ["H shared", "P shared", "C shared", "R shared", "Monotone field"]
RUNG3 = RUNG2 + ["Trials valid", "Surrogate valid", "MR valid"]


def build(extra=()):
    m = Model(5, BINS)
    for c in list(CELLS) + list(extra):
        m.add_row(c[0], c[1], c[2])
    lag = [r["name"] for r in m.rows if "lagged" in r["name"] and "unlagged" not in r["name"]]
    m.sign("H favours X", "H", -1, rung=1)
    m.sign("P favours X", "P", -1, rung=1)
    m.sign("R opposes X", "R", +1, rung=1)
    m.zero("Lag removes P", "P", lag, rung=1)
    m.share("H shared", "H", rung=2)
    m.share("P shared", "P", [r for r in m.rows_with("P") if m.rows[r]["name"] not in lag], rung=2)
    m.share("C shared", "C", rung=2)
    m.share("R shared", "R", rung=2)
    m.monotone("Monotone field", +1, rung=2)
    m.zero("Trials valid", "A", rung=3)
    m.zero("Surrogate valid", "S", rung=3)
    m.zero("MR valid", "Pl", rung=3)
    return m


def true_mean(K, procs, over):
    return float(K @ F_TRUE + sum(over.get(p, B_TRUE[p]) for p in procs))


def simulate(rng):
    studies, cell_y, cell_s, info = [], [], [], []
    for (name, K, procs, k, (a, b), over) in CELLS:
        mu = true_mean(K, procs, over)
        s = rng.uniform(a, b, k)
        ys = mu + rng.normal(0, s)
        m_, se_, inf_ = pool_cell(ys, s)
        cell_y.append(m_); cell_s.append(se_)
        info.append(dict(cell=name, k=k, mu_true=mu, y=m_, se=se_, Q=inf_["Q"], tau2=inf_["tau2"]))
        for yy, ss in zip(ys, s):
            studies.append(dict(cell=name, y=float(yy), s=float(ss), main=bool(K.sum() > 0)))
    return np.array(cell_y), np.array(cell_s), studies, info


FUNCS = {
    "Early-stage population": KE, "Mid-stage population": KM, "Late-stage population": KL,
    "0-3 y": np.eye(5)[0], "3-6 y": np.eye(5)[1], "6-9 y": np.eye(5)[2], "9-12 y": np.eye(5)[3],
    "12-15 y": np.eye(5)[4], "Late minus early population": KL - KE,
}
CLAIMS = [
    ("Benefit in early-stage populations", "Early-stage population", "hi", "<", 0.0),
    ("Benefit in the 0-3 y bin", "0-3 y", "hi", "<", 0.0),
    ("Benefit smaller in late than early populations", "Late minus early population", "lo", ">", 0.0),
    ("Any 12-15 y benefit is weaker than HR 0.90", "12-15 y", "lo", ">", float(np.log(0.90))),
    ("Harm above HR 1.10 excluded at 12-15 y", "12-15 y", "hi", "<", float(np.log(1.10))),
]


def fmt(v):
    if v is None:
        return "none"
    if v == INF:
        return "+inf"
    if v == -INF:
        return "-inf"
    return f"{v:+.3f}"


def hr(v):
    if v is None:
        return "none"
    if not np.isfinite(v):
        return "unbounded"
    return f"{np.exp(v):.2f}"


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    y, s, studies, info = simulate(rng)
    m = build()
    m.set_data(y, s)
    out = dict(seed=SEED, cells=info, studies=len(studies))
    # conventional comparators
    from_main = [d for d in studies if d["main"]]
    ya = np.array([d["y"] for d in from_main]); sa = np.array([d["s"] for d in from_main])
    from decode_core import re_meta, Ledger, decode
    re_all = re_meta(ya, sa)
    tr = [d for d in studies if d["cell"].startswith("Trial")]
    re_tr = re_meta(np.array([d["y"] for d in tr]), np.array([d["s"] for d in tr]))
    out["re_all"] = re_all; out["re_trials"] = re_tr
    # model tests by rung
    tests = {}
    for nm, inc in (("rung1", RUNG1), ("rung2", RUNG2), ("rung3", RUNG3)):
        tests[nm] = model_test(m, inc)
    out["tests"] = tests
    rep = repairs(m, RUNG3, kmax=2, droppable=RUNG3)
    out["repairs_rung3"] = [list(r) for r in rep]
    if not tests["rung3"]["reject"]:
        WORK = list(RUNG3)
    else:
        WORK = [g for g in RUNG3 if g not in rep[0]] if len(rep) == 1 else list(RUNG2)
    out["work_test"] = model_test(m, WORK)
    sets = {"rung1": RUNG1, "rung2": RUNG2, "rung3r": WORK}
    out["working_set"] = WORK
    # capacity (data-free cone criterion) and SIGNET confidence sets
    cap, cs, parts = {}, {}, {}
    for fn, cf in FUNCS.items():
        cap[fn], cs[fn], parts[fn] = {}, {}, {}
        for nm, inc in sets.items():
            cap[fn][nm] = capacity(m, inc, cf)
            r_ = signet_set(m, inc, cf)
            cs[fn][nm] = dict(lo=r_["lo"], hi=r_["hi"], empty=r_["empty"])
            cb = r_["cert"]; cc = r_.get("cc", {})
            parts[fn][nm] = dict(cert_lo=cb["lo"], cert_hi=cb["hi"], cc_lo=cc.get("lo"), cc_hi=cc.get("hi"),
                                 box=box_set(m, inc, cf))
    out["capacity"] = cap; out["cs"] = cs; out["parts"] = parts
    # tolerance at every rung for the claim functionals
    tol = {}
    for fn in sorted(set(c[1] for c in CLAIMS) | {"Mid-stage population", "Late-stage population"}):
        tol[fn] = {}
        for nm, inc in sets.items():
            for k in (1, 2):
                r_ = tolerant_set(m, inc, FUNCS[fn], k=k)
                tol[fn][f"{nm}/k{k}"] = dict(lo=r_["lo"], hi=r_["hi"], empty=r_["empty"])
    out["tolerance"] = tol
    # frontier: which (rung, k) certify each claim
    frontier = {}
    for (cname, fn, side, op, thr) in CLAIMS:
        frontier[cname] = {}
        for nm in sets:
            for k in (0, 1, 2):
                r_ = cs[fn][nm] if k == 0 else tol[fn][f"{nm}/k{k}"]
                v = r_[side]
                # a margin above the bisection tolerance, so that a bound implied by
                # the assumptions alone (for example by monotonicity) is not read as data
                holds = (v is not None) and ((v < thr - 1e-3) if op == "<" else (v > thr + 1e-3))
                frontier[cname][f"{nm}/k{k}"] = bool(holds)
    out["frontier"] = frontier
    # data-independent certificates for the claims at the working rung
    certs = {}
    for (cname, fn, side, op, thr) in CLAIMS:
        ct = certificate(m, WORK, FUNCS[fn], side)
        certs[cname] = None if ct is None else dict(weights=ct["weights"], assumptions=ct["assumptions"],
                                                    slack_weight=ct["slack_weight"], sd=ct["sd"],
                                                    estimate=float(ct["lam"] @ y) * (1 if side == "hi" else -1))
    out["certificates"] = certs
    # identified sets at the TRUE cell means (what infinite precision would give)
    mu_true = np.array([d["mu_true"] for d in info])
    ident = {fn: {nm: identified(m, inc, cf, mu_true) for nm, inc in sets.items()} for fn, cf in FUNCS.items()}
    out["identified_true"] = {fn: {nm: dict(lo=v["lo"], hi=v["hi"], feasible=v["feasible"]) for nm, v in d.items()}
                              for fn, d in ident.items()}
    # DECODE on the same cells, for comparison (three population estimands)
    procs = ["H", "P", "C", "R", "A", "S", "Pl"]
    Z = np.zeros((14, 3)); S = np.zeros((14, 7))
    for i, (name, K, pr, *_r) in enumerate(CELLS):
        if K.sum() > 0:
            j = [np.allclose(K, KE), np.allclose(K, KM), np.allclose(K, KL)].index(True)
            Z[i, j] = 1.0
        for p in pr:
            if p == "P" and "lagged" in name and "unlagged" not in name:
                continue
            S[i, procs.index(p)] = 1.0
    L = Ledger(y, s, Z, S, [0.8, 0.8, 0.3, 0.8, 0.15, 0.3, 0.2], [0.25] * 7,
               estimand_names=["E", "M", "L"], assumption_names=procs)
    post = decode(L)
    out["decode"] = dict(pip=dict(zip(procs, post.pip().tolist())),
                         theta={nm: [post.theta(j).mean(), *post.theta(j).interval()] for j, nm in enumerate(["E", "M", "L"])})
    out["runtime_s"] = time.time() - t0
    json.dump(out, open("worked_example.json", "w"), indent=1, default=float)
    report(out)


def report(o):
    lines = []
    P = lines.append
    P(f"Seed {o['seed']}; {o['studies']} studies in 14 cells")
    for c in o["cells"]:
        P(f"  {c['cell']:<38} k={c['k']}  y={c['y']:+.3f} (SE {c['se']:.3f})  true mean {c['mu_true']:+.3f}  Q={c['Q']:.2f}")
    ra, rt = o["re_all"], o["re_trials"]
    P(f"RE all designs: HR {np.exp(ra['est']):.2f} ({np.exp(ra['lo']):.2f} to {np.exp(ra['hi']):.2f}), p={ra['p']:.2g}, I2={ra['I2']:.2f}, k={ra['k']}")
    P(f"RE trials only: HR {np.exp(rt['est']):.2f} ({np.exp(rt['lo']):.2f} to {np.exp(rt['hi']):.2f}), p={rt['p']:.2g}, k={rt['k']}")
    for nm, t in o["tests"].items():
        P(f"model test {nm}: T={t['T']:.2f}, rank={t['r']}, crit={t['crit']:.2f}, p={t['p']:.3g}, reject={t['reject']}")
    P(f"minimal repairs at rung 3: {o['repairs_rung3']}")
    t = o["work_test"]; P(f"working set test: T={t['T']:.2f}, rank={t['r']}, p={t['p']:.3g}")
    P("95% confidence sets (log HR; HR in brackets)")
    for fn, d in o["cs"].items():
        row = []
        for nm, v in d.items():
            row.append(f"{nm}: [{fmt(v['lo'])}, {fmt(v['hi'])}] ({hr(v['lo'])}-{hr(v['hi'])})")
        P(f"  {fn:<28} " + " | ".join(row))
    P("capacity (bounded below/above) by rung")
    for fn, d in o["capacity"].items():
        P(f"  {fn:<28} " + " | ".join(f"{nm}: {int(v['below'])}{int(v['above'])}" for nm, v in d.items()))
    P("tolerance")
    for fn, d in o["tolerance"].items():
        P(f"  {fn:<28} " + " | ".join(f"{k}: [{fmt(v['lo'])}, {fmt(v['hi'])}]" for k, v in d.items()))
    P("certificate frontier")
    for cn, d in o["frontier"].items():
        P(f"  {cn:<48} " + " ".join(f"{k}:{'Y' if v else '.'}" for k, v in d.items()))
    P("components at each rung: certificate bounds | conditional chi-squared | Bonferroni box")
    for fn, d in o["parts"].items():
        P(f"  {fn:<28} " + " | ".join(f"{nm}: cert [{fmt(v['cert_lo'])}, {fmt(v['cert_hi'])}] cc [{fmt(v['cc_lo'])}, {fmt(v['cc_hi'])}] box [{fmt(v['box'].get('lo'))}, {fmt(v['box'].get('hi'))}]" for nm, v in d.items()))
    P("certificates for the claims (working rung)")
    for cn, d in o["certificates"].items():
        P(f"  {cn}: {d}")
    P("identified sets at the true means")
    for fn, d in o["identified_true"].items():
        P(f"  {fn:<28} " + " | ".join(f"{nm}: [{fmt(v['lo'])}, {fmt(v['hi'])}]" for nm, v in d.items()))
    P(f"DECODE pips: { {k: round(v, 3) for k, v in o['decode']['pip'].items()} }")
    P(f"DECODE theta: { {k: [round(x, 3) for x in v] for k, v in o['decode']['theta'].items()} }")
    P(f"runtime {o['runtime_s']:.1f} s")
    open("worked_example_report.txt", "w").write("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__" and (len(sys.argv) == 1 or sys.argv[1] == "main"):
    main()


# ----------------------------------------------------------------------
# Repeated sampling: does the falsification step find the planted error,
# and does the tolerance dial protect the claims when it does not?
# ----------------------------------------------------------------------
def one_rep(seed):
    rng = np.random.default_rng(seed)
    y, s, studies, info = simulate(rng)
    m = build(); m.set_data(y, s)
    mt = model_test(m, RUNG3, exact_rank=False)
    rep = repairs(m, RUNG3, kmax=1) if mt["reject"] else [()]
    if not mt["reject"]:
        work = list(RUNG3)
    elif len(rep) == 1 and len(rep[0]) == 1:
        work = [g for g in RUNG3 if g not in rep[0]]
    else:
        work = list(RUNG2)
    res = dict(reject=mt["reject"], repairs=[list(r) for r in rep], work_drops=[g for g in RUNG3 if g not in work])
    for fn in ["Early-stage population", "Mid-stage population", "Late-stage population"]:
        truth = float(FUNCS[fn] @ F_TRUE)
        for k in ((0, 1) if fn != "Late-stage population" else (0,)):
            r_ = signet_set(m, work, FUNCS[fn], exact_rank=False) if k == 0 else \
                tolerant_set(m, work, FUNCS[fn], k=1, exact_rank=False)
            lo, hi = (r_["lo"], r_["hi"]) if not r_["empty"] else (None, None)
            res[f"{fn}/k{k}"] = dict(lo=lo, hi=hi, cover=bool(lo is not None and lo <= truth <= hi))
    return res


def replicates(nrep=200):
    from multiprocessing import Pool
    seeds = [SEED + 1000 + i for i in range(nrep)]
    with Pool(2) as pool:
        reps = pool.map(one_rep, seeds)
    summ = dict(n=nrep, reject_rate=float(np.mean([r["reject"] for r in reps])))
    rej = [r for r in reps if r["reject"]]
    summ["surrogate_unique_repair_given_reject"] = float(np.mean([r["repairs"] == [["Surrogate valid"]] for r in rej])) if rej else None
    summ["surrogate_in_some_repair_given_reject"] = float(np.mean([any("Surrogate valid" in x for x in r["repairs"]) for r in rej])) if rej else None
    for fn in ["Early-stage population", "Mid-stage population", "Late-stage population"]:
        for k in ((0, 1) if fn != "Late-stage population" else (0,)):
            key = f"{fn}/k{k}"
            summ[key] = dict(coverage=float(np.mean([r[key]["cover"] for r in reps])),
                             median_width=float(np.median([(r[key]["hi"] - r[key]["lo"]) if r[key]["lo"] is not None else np.inf for r in reps])),
                             benefit_certified=float(np.mean([(r[key]["hi"] is not None and r[key]["hi"] < 0) for r in reps])))
    json.dump(dict(summary=summ, reps=reps), open("worked_example_replicates.json", "w"), indent=1, default=float)
    print(json.dumps(summ, indent=1))


# ----------------------------------------------------------------------
# Design: which new evidence changes what can be concluded?
# ----------------------------------------------------------------------
CANDIDATES = [
    ("Trial restricted to 12-15 y", np.eye(5)[4], ["A"], 0.10, {"A": 0.0}),
    ("Registry (channelled), late", KL, ["R"], 0.09, {}),
    ("Active comparator, late", KL, ["C"], 0.08, {}),
    ("Non-user cohort, lagged, late", KL, ["H", "P"], 0.06, {"P": 0.0}),
    ("Non-user cohort, lagged, early (replicate)", KE, ["H", "P"], 0.06, {"P": 0.0}),
    ("Registry negative control", Z0, ["R"], 0.06, {}),
]


def design(nsim=40):
    rng = np.random.default_rng(SEED)
    y, s, studies, info = simulate(rng)
    m0 = build(); m0.set_data(y, s)
    mt = model_test(m0, RUNG3)
    work = list(RUNG3) if not mt["reject"] else [g for g in RUNG3 if g != "Surrogate valid"]
    fr = Fitter(m0, work).fit()
    x_hat = fr["x"]
    names = m0.param_names()
    plug = {nm: v for nm, v in zip(names, x_hat)}
    f_hat = x_hat[:5]
    rows = []
    rng2 = np.random.default_rng(SEED + 7)
    for (nm, K, procs, se, over) in CANDIDATES:
        m = build(extra=[(nm, K, procs, 1, (se, se), over)])
        # plug-in mean of the new cell: field part plus the fitted shared bias of each process
        mean_new = float(K @ f_hat)
        for p in procs:
            if p in over:
                continue
            vals = [plug[f"{p}@{r['name']}"] for r in m0.rows if f"{p}@{r['name']}" in plug and
                    not ("lagged" in r["name"] and "unlagged" not in r["name"] and p == "P")]
            mean_new += float(np.mean(vals)) if vals else 0.0
        caps = {}
        for fn in ["Late-stage population", "12-15 y"]:
            for rn, inc in (("rung1", RUNG1), ("rung2", RUNG2), ("rung3", work)):
                caps[f"{fn}/{rn}"] = capacity(m, inc, FUNCS[fn])
        widths = {"Late-stage population": [], "12-15 y": []}
        for _ in range(nsim):
            y_new = mean_new + rng2.normal(0, se)
            m.set_data(np.r_[y, y_new], np.r_[s, se])
            for fn in widths:
                r_ = signet_set(m, work, FUNCS[fn], exact_rank=False)
                widths[fn].append((r_["hi"] - r_["lo"]) if not r_["empty"] else np.nan)
        rows.append(dict(candidate=nm, se=se, plug_mean=mean_new, capacity=caps,
                         width={fn: float(np.nanmedian(v)) for fn, v in widths.items()}))
    base = {fn: float(signet_set(m0, work, FUNCS[fn])["hi"] - signet_set(m0, work, FUNCS[fn])["lo"])
            for fn in ["Late-stage population", "12-15 y"]}
    json.dump(dict(base_width=base, rows=rows), open("worked_example_design.json", "w"), indent=1, default=float)
    print("base widths", base)
    for r in rows:
        print(f"{r['candidate']:<44} plug {r['plug_mean']:+.3f}  widths " +
              ", ".join(f"{k}: {v:.3f}" for k, v in r["width"].items()))


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "rep":
    replicates(int(sys.argv[2]) if len(sys.argv) > 2 else 200)
if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "design":
    design()
