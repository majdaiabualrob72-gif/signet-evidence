"""
paired_contrasts.py: paired re-analysis of the method contrasts (manuscript Table 4) and summary
of the v1.1.0 validation experiment.

Every contrast compares two methods on the same simulated literatures (same generator and seeds),
now using one record per literature: a = both succeed, b = only A succeeds, c = only B succeeds.
Difference in proportions (b - c) / n with the Newcombe paired-data interval (method 10, which
combines Wilson intervals for the two marginal proportions with a correlation correction) and
the exact two-sided McNemar test (binomial test on the b + c discordant literatures).

Inputs : experiments/simulation/benchmark_reps_{i}.json, benchmark_v2_stress_reps_{i}.json (archived)
         experiments/validation/reps_{part}.json (paired_reruns.py) and validation_reps_{i}.json
Outputs: experiments/validation/paired_contrasts.md, paired_contrasts.json, validation_summary.md
"""
import json
import math
import os
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SIM = os.path.join(ROOT, "experiments", "simulation")
sys.path.insert(0, SIM)
import benchmark as B          # noqa: E402

Z = 1.959963984540054
MAIN = {sc: i for i, sc in enumerate(B.ALL_SCEN)}
STRESS = ["extreme outlier", "publication bias", "heterogeneous populations",
          "missing controls and variances", "contradictory evidence"]
K1, RE, PR, BAS = "SIGNET R3 k=1", "RE (all designs)", "SIGNET protocol", "BA-S (elicited directions and magnitudes)"


def wilson(x, n):
    p = x / n
    d = 1 + Z ** 2 / n
    c = p + Z ** 2 / (2 * n)
    h = Z * math.sqrt(p * (1 - p) / n + Z ** 2 / (4 * n ** 2))
    return (c - h) / d, (c + h) / d


def newcombe_paired(a, b, c, d):
    """Newcombe (1998), paired data, method 10."""
    n = a + b + c + d
    p1, p2 = (a + b) / n, (a + c) / n
    l1, u1 = wilson(a + b, n)
    l2, u2 = wilson(a + c, n)
    den = (a + b) * (c + d) * (a + c) * (b + d)
    if den == 0:
        phi = 0.0
    else:
        num = a * d - b * c
        phi = (max(num - n / 2, 0.0) if num > 0 else num) / math.sqrt(den)
    th = p1 - p2
    lo = th - math.sqrt(max((p1 - l1) ** 2 - 2 * phi * (p1 - l1) * (u2 - p2) + (u2 - p2) ** 2, 0.0))
    hi = th + math.sqrt(max((p2 - l2) ** 2 - 2 * phi * (p2 - l2) * (u1 - p1) + (u1 - p1) ** 2, 0.0))
    return th, lo, hi


def mcnemar_exact(b, c):
    if b + c == 0:
        return 1.0
    return float(min(1.0, 2 * stats.binom.cdf(min(b, c), b + c, 0.5)))


def load(fn):
    return json.load(open(fn))


def main_reps(sc):
    return load(os.path.join(SIM, f"benchmark_reps_{MAIN[sc]}.json"))


def stress_reps(sc):
    return load(os.path.join(SIM, f"benchmark_v2_stress_reps_{STRESS.index(sc)}.json"))


def val(part):
    fn = os.path.join(HERE, f"reps_{part}.json")
    return load(fn) if os.path.exists(fn) else None


def outcome(rec, what):
    if rec is None or rec.get("absent"):
        return None
    return bool(rec.get(what, False))


def contrast(label, A, Bv, what):
    """A, Bv: lists of per-literature records (aligned by seed) for methods A and B."""
    n = min(len(A), len(Bv))
    pairs = [(outcome(A[i], what), outcome(Bv[i], what)) for i in range(n)]
    pairs = [p for p in pairs if None not in p]
    a = sum(1 for x, y in pairs if x and y)
    b = sum(1 for x, y in pairs if x and not y)
    c = sum(1 for x, y in pairs if y and not x)
    d = sum(1 for x, y in pairs if not x and not y)
    m = a + b + c + d
    th, lo, hi = newcombe_paired(a, b, c, d)
    p = mcnemar_exact(b, c)
    return dict(label=label, n=m, xA=a + b, xB=a + c, a=a, b=b, c=c, d=d, diff=th, lo=lo, hi=hi, p=p)


def fmt_row(r):
    ps = "< 0.0001" if r["p"] < 1e-4 else f"{r['p']:.4f}"
    return (f"| {r['label']} | {r['xA']}/{r['n']} ({r['xA'] / r['n']:.3f}) | {r['xB']}/{r['n']} ({r['xB'] / r['n']:.3f}) | "
            f"{r['b']} / {r['c']} | {r['diff']:+.3f} ({r['lo']:+.3f} to {r['hi']:+.3f}) | {ps} |")


def col(reps, key):
    return [r[key] for r in reps]


def contrasts():
    rows = []
    r0 = main_reps("reference")
    rows.append(contrast("Coverage, reference: SIGNET k=1 vs RE", col(r0, K1), col(r0, RE), "cover"))
    rows.append(contrast("Coverage, reference: SIGNET k=1 vs adaptive protocol", col(r0, K1), col(r0, PR), "cover"))
    r4 = main_reps("15% of signs wrong")
    rows.append(contrast("Coverage, 15% of signs wrong: SIGNET k=1 vs adaptive protocol", col(r4, K1), col(r4, PR), "cover"))
    sh = stress_reps("heterogeneous populations")
    rows.append(contrast("Coverage, heterogeneous populations: SIGNET k=1 vs adaptive protocol", col(sh, K1), col(sh, PR), "cover"))
    cmp_ = val("comparators")
    if cmp_:
        rows.append(contrast("Benefit certified, reference (true benefit): SIGNET k=1 vs BA-S",
                             col(r0, K1), col(cmp_["reference"], BAS), "benefit"))
    r1 = main_reps("null effect")
    rows.append(contrast("False benefit, null effect: SIGNET k=1 vs RE", col(r1, K1), col(r1, RE), "benefit"))
    sc2 = "one-signed literature, null effect, magnitudes doubled"
    r9 = main_reps(sc2)
    if cmp_:
        rows.append(contrast("False benefit, one-signed literature, doubled biases: SIGNET k=1 vs BA-S",
                             col(r9, K1), col(cmp_[sc2], BAS), "benefit"))
    rows.append(contrast("False benefit, one-signed literature, doubled biases: SIGNET k=1 vs RE", col(r9, K1), col(r9, RE), "benefit"))
    sw = val("sweep")
    if sw:
        k4, k3 = "mixed designs, true benefit | x4.0", "one-signed literature, null effect | x3.0"
        rows.append(contrast("Coverage, biases x4 (mixed designs): SIGNET k=1 vs BA-S", col(sw[k4], K1), col(sw[k4], BAS), "cover"))
        rows.append(contrast("False benefit, biases x3 (one-signed): SIGNET k=1 vs BA-S", col(sw[k3], K1), col(sw[k3], BAS), "benefit"))
    sc_ = stress_reps("contradictory evidence")
    rows.append(contrast("Coverage, contradictory evidence: SIGNET k=1 vs RE", col(sc_, K1), col(sc_, RE), "cover"))
    ha = val("hetallow")
    if ha:
        rows.append(contrast("Coverage, heterogeneous populations: SIGNET R3 with vs without allowance",
                             col(ha["heterogeneous populations"], "SIGNET R3 + heterogeneity allowance"), col(sh, "SIGNET R3"), "cover"))
    bd = val("bounded")
    if bd:
        rows.append(contrast("Benefit certified, reference: SIGNET-B k=1 vs SIGNET k=1",
                             col(bd["reference"], "SIGNET-B k=1"), col(r0, K1), "benefit"))
    cr = val("cellrobust")
    if cr:
        m = "SIGNET R3 k=1 (group or cell)"
        rows.append(contrast("Coverage, contradictory evidence: group-or-cell tolerance vs group tolerance",
                             col(cr["contradictory evidence"], m), col(sc_, K1), "cover"))
        rows.append(contrast("Coverage, contradictory evidence: group-or-cell tolerance vs RE",
                             col(cr["contradictory evidence"], m), col(sc_, RE), "cover"))
        rows.append(contrast("Benefit certified, reference: group-or-cell tolerance vs group tolerance",
                             col(cr["reference"], m), col(r0, K1), "benefit"))
        rows.append(contrast("False benefit, null effect: group-or-cell tolerance vs RE",
                             col(cr["null effect"], m), col(r1, RE), "benefit"))
    return rows


def reproduction_checks():
    """Do the re-runs reproduce the archived summaries?"""
    out = []

    def agg(reps, key, what):
        v = [r[key][what] for r in reps if not r[key].get("absent")]
        return sum(v), len(v)
    cmp_ = val("comparators")
    if cmp_:
        arch = load(os.path.join(SIM, "benchmark_v2_comparators.json"))["results"]
        for sc, reps in cmp_.items():
            for mth in ("FE (all designs)", "BA-U (unsigned bias priors)", BAS, "NCO-cal (calibrated, RE)"):
                x, n = agg(reps, mth, "cover")
                out.append(("comparators", sc, mth, x / n, arch[sc][mth]["coverage"]))
    cr = val("cellrobust")
    V100 = os.path.join(SIM, "v1.0.0_missing_controls")      # the version 1.0.0 archive, before the generator fix
    if cr:
        arch = load(os.path.join(V100, "benchmark_v2_cellrobust.json"))["results"]
        for sc, reps in cr.items():
            for what, key in (("cover", "coverage"), ("benefit", "benefit")):
                x, n = agg(reps, "SIGNET R3 k=1 (group or cell)", what)
                out.append(("cellrobust", sc, f"group or cell, {key}", x / n, arch[sc]["SIGNET R3 k=1 (group or cell)"][key]))
    sw = val("sweep")
    if sw:
        arch = load(os.path.join(SIM, "benchmark_v2_sweep.json"))["results"]
        for key, reps in sw.items():
            for mth in (K1, BAS):
                x, n = agg(reps, mth, "cover")
                out.append(("sweep", key, mth, x / n, arch[key][mth]["coverage"]))
    ha = val("hetallow")
    if ha:
        arch = load(os.path.join(SIM, "benchmark_v2_hetallow.json"))["results"]
        x, n = agg(ha["heterogeneous populations"], "SIGNET R3 + heterogeneity allowance", "cover")
        out.append(("hetallow", "heterogeneous populations", "R3 + allowance", x / n,
                    arch["heterogeneous populations"]["SIGNET R3 + heterogeneity allowance"]["coverage"]))
    bd = val("bounded")
    if bd:
        arch = load(os.path.join(V100, "benchmark_v2_bounded.json"))["results"]
        x, n = agg(bd["reference"], "SIGNET-B k=1", "benefit")
        out.append(("bounded", "reference", "SIGNET-B k=1 benefit", x / n, arch["reference"]["SIGNET-B k=1"]["benefit"]))
    return out


def validation_summary():
    L = ["# Validation of signet v1.1.0 on the archived simulated literatures (generated by paired_contrasts.py)", "",
         "Settings: legacy = archived (estimated V, approximate rank, 24-point grid, tol 5e-4, no dense check); "
         "hardened = exact rank, 40-point grid, tol 1e-4, 16-point dense check, bracketing to 1e4; outer = n-degree outer set; "
         "frozen = certificate weights chosen with within-study (fixed-effect) variances; knownV = cells pooled with the true "
         "heterogeneity and analysed with their exact covariance. Differences: paired, Newcombe method 10 interval, exact McNemar p.", ""]
    rows, agree = [], []
    head = ("| Scenario | n | Legacy R3 coverage | Archived R3 coverage | Agreement with archive (R1, R3) | Hardened R3 | Outer R3 | "
            "Frozen R3 | Known-V R3 | Legacy R1 | Frozen R1 | Known-V R1 | Certificate unique, R3 (hi, lo) | Hull extended | v1.0.0 truncation |")
    L += [head, "|" + "---|" * head.count("|")[:-1] if False else "|" + "|".join(["---"] * (head.count("|") - 1)) + "|"]
    allp = {}
    for i, sc in enumerate(B.ALL_SCEN):
        fn = os.path.join(HERE, f"validation_reps_{i}.json")
        if not os.path.exists(fn):
            continue
        V = load(fn)["reps"]
        arch = main_reps(sc)
        n = len(V)
        cov = lambda key: np.mean([r[key]["cover"] for r in V])
        ag = []
        for rn in ("R1", "R3"):
            same = 0
            for r, a in zip(V, arch):
                x, y = r[f"legacy {rn}"], a[f"SIGNET {rn}"]
                ok = x["cover"] == y["cover"] and x.get("benefit", False) == y["benefit"] and x.get("harm", False) == y["harm"]
                if ok and y["width"] is not None and not y["empty"]:
                    w = (float(x["hi"]) - float(x["lo"]))
                    ok = abs(w - y["width"]) < 2e-3
                same += ok
            ag.append(same)
        agree.append((sc, ag, n))
        cu = [r["cert unique R3"] for r in V]
        uh = np.mean([c["hi"] is True for c in cu]); ul = np.mean([c["lo"] is True for c in cu])
        ext = sum(1 for r in V if r["hardened flags"]["hull_extended"])
        trunc = sum(1 for r in V if any((r["hardened flags"]["legacy_truncation"] or {}).values()))
        L.append(f"| {sc} | {n} | {cov('legacy R3'):.3f} | {np.mean([a['SIGNET R3']['cover'] for a in arch[:n]]):.3f} | "
                 f"{ag[0]}/{n}, {ag[1]}/{n} | {cov('hardened R3'):.3f} | {cov('outer R3'):.3f} | {cov('frozen R3'):.3f} | "
                 f"{cov('knownV R3'):.3f} | {cov('legacy R1'):.3f} | {cov('frozen R1'):.3f} | {cov('knownV R1'):.3f} | "
                 f"{uh:.2f}, {ul:.2f} | {ext} | {trunc} |")
        for var in ("hardened R3", "outer R3", "frozen R3", "knownV R3", "frozen R1", "knownV R1"):
            base = "legacy R1" if var.endswith("R1") else "legacy R3"
            for what in ("cover", "benefit"):
                r = contrast(f"{sc}: {var} vs {base} ({what})", [x[var] for x in V], [x[base] for x in V], what)
                allp.setdefault(var, []).append(r)
    L += ["", "## Paired differences from the legacy settings (coverage and certified benefit)", "",
          "| Contrast | A x/n | B x/n | Discordant (A only / B only) | Difference (95% CI) | p |", "|---|---|---|---|---|---|"]
    for var, lst in allp.items():
        for r in lst:
            L.append(fmt_row(r))
    widths = []
    for i, sc in enumerate(B.ALL_SCEN):
        fn = os.path.join(HERE, f"validation_reps_{i}.json")
        if not os.path.exists(fn):
            continue
        V = load(fn)["reps"]
        def w(key):
            v = [float(r[key]["hi"]) - float(r[key]["lo"]) for r in V if not r[key]["empty"]
                 and np.isfinite(float(r[key]["hi"])) and np.isfinite(float(r[key]["lo"]))]
            return np.median(v) if v else float("nan")
        widths.append(f"| {sc} | {w('legacy R3'):.3f} | {w('hardened R3'):.3f} | {w('outer R3'):.3f} | {w('frozen R3'):.3f} | {w('knownV R3'):.3f} |")
    L += ["", "## Median width of bounded rung-3 sets (log scale)", "",
          "| Scenario | Legacy | Hardened | Outer | Frozen | Known V |", "|---|---|---|---|---|---|"] + widths
    return "\n".join(L) + "\n", agree, allp


def main():
    rows = contrasts()
    head = ["| Contrast | Method A, x/n (proportion) | Method B, x/n (proportion) | Discordant literatures (A only / B only) | "
            "Difference A minus B (95% CI) | p |", "|---|---|---|---|---|---|"]
    txt = "\n".join(head + [fmt_row(r) for r in rows]) + "\n"
    rc = reproduction_checks()
    rep = ["", "## Reproduction of the archived (version 1.0.0) summaries by the re-runs", "",
           "The missing-controls stress scenario was regenerated after a generator fix (Supplementary S11.12); its rows compare the regenerated values with the version 1.0.0 archive and are not expected to agree.", "",
           "| Part | Scenario | Method | Re-run | Archived |", "|---|---|---|---|---|"]
    rep += [f"| {a} | {b} | {c} | {d:.4f} | {e:.4f} |" for a, b, c, d, e in rc]
    doc = __doc__.strip()
    open(os.path.join(HERE, "paired_contrasts.md"), "w").write(
        "# Method contrasts, paired (paired_contrasts.py)\n\n" + doc + "\n\n" + txt + "\n".join(rep) + "\n")
    json.dump(dict(contrasts=rows, reproduction=rc), open(os.path.join(HERE, "paired_contrasts.json"), "w"), indent=1)
    vs, agree, allp = validation_summary()
    open(os.path.join(HERE, "validation_summary.md"), "w").write(vs)
    json.dump(dict(agreement=agree, paired={k: v for k, v in allp.items()}),
              open(os.path.join(HERE, "validation_summary.json"), "w"), indent=1)
    print(txt)
    print("\n".join(rep))
    print(vs)


if __name__ == "__main__":
    main()
