"""N11 of the forensic audit: how often does the certificate rule give zero weight to the most precise cell
that measures the effect, and what does the rule cost in standard deviation?

Rule of the paper (Theorem 2): among all certificates, first minimise the total weight on inequality assumptions,
then minimise the variance. The comparison certificate here minimises the variance alone.
Run from the repository root:  python3 experiments/audit_v12/certificate_weights.py
Writes certificate_weights.json and certificate_weights.md in this folder."""
import json
import os
import sys

import numpy as np
import cvxpy as cp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "experiments", "simulation"))
import benchmark as B                      # noqa: E402
import benchmark_v2 as B2                  # noqa: E402
import sim_audit as S                      # noqa: E402
import realdata_audit_models as RM         # noqa: E402
from signet.core import certificate       # noqa: E402
from scipy import stats as sst             # noqa: E402

Z = float(sst.norm.ppf(1 - 0.0125))        # the certificate limits of the reported set use alpha / 4 on each side
EQ = 1e-4                                  # limits closer than this on the log scale are counted as equal


def min_variance_certificate(model, include, cf, side):
    """Certificate of least variance, without the first stage of the rule. Returns (weights, standard deviation) or None."""
    c_ = model.compile(include)
    c = model.cvec(cf) * (1.0 if side == "hi" else -1.0)
    A, D, E, V = c_["A"], c_["D"], c_["E"], model.V
    n, mD, mE = A.shape[0], D.shape[0], E.shape[0]
    lam = cp.Variable(n)
    expr = A.T @ lam
    cons = []
    if mD:
        nu = cp.Variable(mD, nonneg=True); expr = expr + D.T @ nu
    if mE:
        kap = cp.Variable(mE); expr = expr + E.T @ kap
    cons.append(expr == c)
    Lc = np.linalg.cholesky(V + 1e-14 * np.eye(n))
    prob = cp.Problem(cp.Minimize(cp.sum_squares(Lc.T @ lam)), cons)
    try:
        prob.solve(solver=cp.CLARABEL)
    except Exception:
        return None
    if prob.status not in ("optimal", "optimal_inaccurate"):
        return None
    l = np.array(lam.value).ravel()
    return l, float(np.sqrt(max(l @ V @ l, 0.0)))


def one_model(m, include, cf):
    """For each side with a certificate: is the most precise effect-measuring cell given zero weight, and the
    ratio of the rule's standard deviation to the least possible."""
    K = np.array([float(np.dot(r["K"], cf)) for r in m.rows])
    sd_cell = np.sqrt(np.diag(m.V))
    eff = [i for i in range(m.n) if abs(K[i]) > 0]
    if not eff:
        return {}
    best = min(eff, key=lambda i: sd_cell[i] / abs(K[i]))
    out = {}
    for side in ("hi", "lo"):
        ce = certificate(m, include, cf, side)
        if ce is None:
            continue
        mv = min_variance_certificate(m, include, cf, side)
        # one-sided limit at level alpha / 4 = 0.0125 under each rule, on the side's own scale (smaller is tighter)
        lim_rule = float(ce["lam"] @ m.y + Z * ce["sd"])
        lim_mv = float(mv[0] @ m.y + Z * mv[1]) if mv is not None else None
        out[side] = dict(zero_on_most_precise=bool(abs(ce["lam"][best]) <= 1e-6), sd_rule=ce["sd"],
                         sd_minvar=(mv[1] if mv is not None else None), slack_weight=ce["slack_weight"],
                         limit_rule=lim_rule, limit_minvar=lim_mv,
                         n_effect_cells=len(eff), most_precise=m.rows[best]["name"],
                         weights=ce["weights"])
    return out


def simulation(n_rep=400):
    res = {}
    for sc in S.SCEN15:
        rows = {rn: dict(hi=[], lo=[]) for rn in ("R1", "R2", "R3")}
        for r in range(n_rep):
            rng = np.random.default_rng(S.seed_of(sc, r))
            theta, b, present, studies = B2.literature_v2(sc, rng)
            m, R = B.signet_model(present, studies)
            for rn in rows:
                o = one_model(m, R[rn], [1.0])
                for side in o:
                    if o[side]["n_effect_cells"] >= 2:
                        rows[rn][side].append((o[side]["zero_on_most_precise"],
                                               o[side]["sd_rule"] / o[side]["sd_minvar"] if o[side]["sd_minvar"] else None,
                                               (o[side]["limit_rule"] - o[side]["limit_minvar"]) if o[side]["limit_minvar"] is not None else None))
        res[sc] = {}
        for rn in rows:
            for side in ("hi", "lo"):
                x = rows[rn][side]
                ratios = [q for (_z, q, _d) in x if q is not None]
                diffs = [d for (_z, _q, d) in x if d is not None]
                res[sc][f"{rn} {side}"] = dict(
                    n=len(x), zero=int(sum(z for (z, _q, _d) in x)),
                    prop_zero=(float(np.mean([z for (z, _q, _d) in x])) if x else None),
                    # limit under the rule minus limit under least variance; positive means the rule gave the looser limit
                    mean_limit_diff=(float(np.mean(diffs)) if diffs else None),
                    median_limit_diff=(float(np.median(diffs)) if diffs else None),
                    # a difference below 1e-4 on the log scale (the bisection tolerance of the package) counts as equal
                    prop_rule_tighter=(float(np.mean([d < -EQ for d in diffs])) if diffs else None),
                    prop_rule_looser=(float(np.mean([d > EQ for d in diffs])) if diffs else None),
                    prop_same_limit=(float(np.mean([abs(d) <= EQ for d in diffs])) if diffs else None),
                    records=[[int(z), (None if q is None else round(q, 6)), (None if d is None else round(d, 6))] for (z, q, d) in x],
                    median_ratio=(float(np.median(ratios)) if ratios else None),
                    p90_ratio=(float(np.quantile(ratios, 0.9)) if ratios else None),
                    prop_ratio_above_1_05=(float(np.mean([q > 1.05 for q in ratios])) if ratios else None))
        print(sc, {k: (v["n"], v["zero"]) for k, v in res[sc].items()}, flush=True)
    return res


def realdata():
    out = {}
    for (label, (M, R), cf) in RM.ANALYSES():
        if "total cancer" in label:        # control outcome: bounded by assumption at rungs 2 and 3, no data weight
            continue
        out[label] = {}
        for rn in R:
            o = one_model(M, R[rn], cf)
            out[label][rn] = {s: dict(zero_on_most_precise=v["zero_on_most_precise"], most_precise=v["most_precise"],
                                      weights=v["weights"], sd_rule=v["sd_rule"], sd_minvar=v["sd_minvar"],
                                      limit_rule=v["limit_rule"], limit_minvar=v["limit_minvar"], slack_weight=v["slack_weight"],
                                      n_effect_cells=v["n_effect_cells"]) for s, v in o.items()}
    return out


def main():
    n_rep = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    out = dict(n_rep=n_rep, simulation=simulation(n_rep), realdata=realdata())
    L = ["# Certificate weights (finding N11)", "",
         "Rule of the paper: least total weight on inequality assumptions first, least variance second.",
         "Comparison: the certificate of least variance alone. Literatures with at least two cells measuring the effect.",
         f"Simulated literatures per scenario: {n_rep} (archived seeds).", "",
         "## Simulation, rung 3, upper limit (the side that supports a claim of benefit)", "",
         "The last three columns compare the one-sided limits (level 0.0125) of the two rules on the log scale: mean of",
         "(limit under the rule minus limit under least variance); a negative value means that the rule gave the tighter limit on average.",
         "Limits closer than 1e-4 on the log scale are counted as the same.", "",
         "| Scenario | Literatures with a certificate | Most precise cell given zero weight | Median SD ratio (rule / least variance) | 90th centile | Ratio above 1.05 | Mean difference in limit | Rule tighter | Same limit | Rule looser |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for sc, v in out["simulation"].items():
        x = v["R3 hi"]
        if x["n"]:
            L.append(f"| {sc} | {x['n']} | {x['zero']} ({x['prop_zero']:.3f}) | {x['median_ratio']:.3f} | {x['p90_ratio']:.3f} | {x['prop_ratio_above_1_05']:.3f} | {x['mean_limit_diff']:+.3f} | {x['prop_rule_tighter']:.3f} | {x['prop_same_limit']:.3f} | {x['prop_rule_looser']:.3f} |")
        else:
            L.append(f"| {sc} | 0 | not applicable | | | | | | | |")
    L += ["", "## Pooled over the 15 scenarios", "",
          "| Rung and side | Literatures with a certificate | Most precise cell given zero weight | SD ratio above 1.05 | Mean difference in limit | Rule tighter | Same limit | Rule looser |", "|---|---|---|---|---|---|---|---|"]
    out["pooled"] = {}
    for k in ("R1 hi", "R1 lo", "R2 hi", "R2 lo", "R3 hi", "R3 lo"):
        recs = [r for v in out["simulation"].values() for r in v[k]["records"]]
        dd = [r[2] for r in recs if r[2] is not None]; qq = [r[1] for r in recs if r[1] is not None]
        pz = dict(n=len(recs), zero=int(sum(r[0] for r in recs)), prop_zero=float(np.mean([r[0] for r in recs])),
                  prop_ratio_above_1_05=float(np.mean([q > 1.05 for q in qq])), median_ratio=float(np.median(qq)),
                  mean_limit_diff=float(np.mean(dd)), prop_rule_tighter=float(np.mean([d < -EQ for d in dd])),
                  prop_same_limit=float(np.mean([abs(d) <= EQ for d in dd])), prop_rule_looser=float(np.mean([d > EQ for d in dd])))
        out["pooled"][k] = pz
        L.append(f"| {k.replace('hi', 'upper').replace('lo', 'lower').replace('R', 'Rung ')} | {pz['n']} | {pz['zero']} ({pz['prop_zero']:.3f}) | {pz['prop_ratio_above_1_05']:.3f} | {pz['mean_limit_diff']:+.3f} | {pz['prop_rule_tighter']:.3f} | {pz['prop_same_limit']:.3f} | {pz['prop_rule_looser']:.3f} |")
    L += ["", "## Simulation, all rungs and sides", "",
          "| Scenario | Rung and side | Literatures with a certificate | Zero weight on the most precise cell | Median SD ratio | Ratio above 1.05 | Mean difference in limit |", "|---|---|---|---|---|---|---|"]
    for sc, v in out["simulation"].items():
        for k, x in v.items():
            if x["n"]:
                L.append(f"| {sc} | {k} | {x['n']} | {x['zero']} ({x['prop_zero']:.3f}) | {x['median_ratio']:.3f} | {x['prop_ratio_above_1_05']:.3f} | {x['mean_limit_diff']:+.3f} |")
    L += ["", "## Real-data series", "",
          "The control outcome of the hormone-therapy question is omitted: at rungs 2 and 3 it is bounded by assumption.", "",
          "| Question | Rung | Side | Cells measuring the target | Most precise | Zero weight on it | Weights | SD under the rule | Least SD |", "|---|---|---|---|---|---|---|---|---|"]
    for label, v in out["realdata"].items():
        for rn, o in v.items():
            for side, x in o.items():
                w = ", ".join(f"{k}: {val:+.3f}" for k, val in x["weights"].items())
                L.append(f"| {label} | {rn} | {'upper' if side == 'hi' else 'lower'} | {x['n_effect_cells']} | {x['most_precise']} | {'yes' if x['zero_on_most_precise'] else 'no'} | {w} | {x['sd_rule']:.4f} | {x['sd_minvar']:.4f} |")
    json.dump(out, open(os.path.join(HERE, "certificate_weights.json"), "w"))
    open(os.path.join(HERE, "certificate_weights.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L[:30]))


if __name__ == "__main__":
    main()
