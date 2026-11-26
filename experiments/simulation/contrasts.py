"""Method contrasts reported in the manuscript (Section 7). Chosen after the results were seen (post hoc).

Each contrast is a difference between two proportions (coverage or rate of certified benefit),
estimated from archived summaries. The methods were applied to the same simulated literatures
(identical seeds and generator), but the summaries do not keep the pairing, so the samples are
treated as independent. This is conservative for positively correlated outcomes.
Interval: Newcombe hybrid score interval (method 10, Newcombe 1998). P-value: two-sided Fisher exact test.
"""
import json
import os

from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
Z = 1.959963984540054


def wilson(x, n):
    p = x / n
    d = 1 + Z ** 2 / n
    c = p + Z ** 2 / (2 * n)
    h = Z * ((p * (1 - p) / n + Z ** 2 / (4 * n ** 2)) ** 0.5)
    return (c - h) / d, (c + h) / d


def newcombe(x1, n1, x2, n2):
    p1, p2 = x1 / n1, x2 / n2
    l1, u1 = wilson(x1, n1)
    l2, u2 = wilson(x2, n2)
    d = p1 - p2
    lo = d - ((p1 - l1) ** 2 + (u2 - p2) ** 2) ** 0.5
    hi = d + ((u1 - p1) ** 2 + (p2 - l2) ** 2) ** 0.5
    return d, lo, hi


def get(res, sc, m, what):
    v = res[sc][m]
    return round(v[what] * v["n"]) if "n" in v and v["n"] else round(v[what] * 400), v.get("n") or 400


def contrast(label, a, b):
    (x1, n1), (x2, n2) = a, b
    d, lo, hi = newcombe(x1, n1, x2, n2)
    p = stats.fisher_exact([[x1, n1 - x1], [x2, n2 - x2]])[1]
    ps = "< 0.0001" if p < 1e-4 else f"{p:.4f}"
    return (f"| {label} | {x1}/{n1} ({x1 / n1:.3f}) | {x2}/{n2} ({x2 / n2:.3f}) | "
            f"{d:+.3f} ({lo:+.3f} to {hi:+.3f}) | {ps} |")


def main():
    base = json.load(open("all_scenarios_summary.json"))
    sw = json.load(open("benchmark_v2_sweep.json"))["results"]
    het = json.load(open("benchmark_v2_hetallow.json"))["results"]
    rows = []
    K1, RE, BAS, PR = "SIGNET R3 k=1", "RE (all designs)", "BA-S (elicited directions and magnitudes)", "SIGNET protocol"

    def g(sc, m, w="coverage", src=base, n=400):
        v = src[sc][m]
        nn = v.get("n") or n
        return round(v[w] * nn), nn

    rows.append(contrast("Coverage, reference: SIGNET k=1 vs RE", g("reference", K1), g("reference", RE)))
    rows.append(contrast("Coverage, reference: SIGNET k=1 vs adaptive protocol", g("reference", K1), g("reference", PR)))
    rows.append(contrast("Coverage, 15% of signs wrong: SIGNET k=1 vs adaptive protocol",
                         g("15% of signs wrong", K1), g("15% of signs wrong", PR)))
    rows.append(contrast("Coverage, heterogeneous populations: SIGNET k=1 vs adaptive protocol",
                         g("heterogeneous populations", K1), g("heterogeneous populations", PR)))
    rows.append(contrast("Benefit certified, reference (true benefit): SIGNET k=1 vs BA-S",
                         g("reference", K1, "benefit"), g("reference", BAS, "benefit")))
    rows.append(contrast("False benefit, null effect: SIGNET k=1 vs RE",
                         g("null effect", K1, "benefit"), g("null effect", RE, "benefit")))
    sc2 = "one-signed literature, null effect, magnitudes doubled"
    rows.append(contrast("False benefit, one-signed literature, doubled biases: SIGNET k=1 vs BA-S",
                         g(sc2, K1, "benefit"), g(sc2, BAS, "benefit")))
    rows.append(contrast("False benefit, one-signed literature, doubled biases: SIGNET k=1 vs RE",
                         g(sc2, K1, "benefit"), g(sc2, RE, "benefit")))
    k4 = "mixed designs, true benefit | x4.0"
    rows.append(contrast("Coverage, biases x4 (mixed designs): SIGNET k=1 vs BA-S",
                         g(k4, K1, src=sw, n=200), g(k4, BAS, src=sw, n=200)))
    k3 = "one-signed literature, null effect | x3.0"
    rows.append(contrast("False benefit, biases x3 (one-signed): SIGNET k=1 vs BA-S",
                         g(k3, K1, "benefit", src=sw, n=200), g(k3, BAS, "benefit", src=sw, n=200)))
    rows.append(contrast("Coverage, contradictory evidence: SIGNET k=1 vs RE",
                         g("contradictory evidence", K1), g("contradictory evidence", RE)))
    rows.append(contrast("Coverage, heterogeneous populations: SIGNET R3 with vs without allowance",
                         g("heterogeneous populations", "SIGNET R3 + heterogeneity allowance", src=het),
                         g("heterogeneous populations", "SIGNET R3")))
    try:
        bnd = json.load(open("benchmark_v2_bounded.json"))["results"]
        rows.append(contrast("Benefit certified, reference: SIGNET-B k=1 vs SIGNET k=1",
                             g("reference", "SIGNET-B k=1", "benefit", src=bnd), g("reference", K1, "benefit")))
    except (FileNotFoundError, KeyError):
        pass
    try:
        cr = json.load(open("benchmark_v2_cellrobust.json"))["results"]
        m = "SIGNET R3 k=1 (group or cell)"
        rows.append(contrast("Coverage, contradictory evidence: group-or-cell tolerance vs group tolerance",
                             g("contradictory evidence", m, src=cr), g("contradictory evidence", K1)))
        rows.append(contrast("Coverage, contradictory evidence: group-or-cell tolerance vs RE",
                             g("contradictory evidence", m, src=cr), g("contradictory evidence", RE)))
        if "reference" in cr:
            rows.append(contrast("Benefit certified, reference: group-or-cell tolerance vs group tolerance",
                                 g("reference", m, "benefit", src=cr), g("reference", K1, "benefit")))
    except (FileNotFoundError, KeyError):
        pass
    head = ["| Contrast | Method A, x/n (proportion) | Method B, x/n (proportion) | Difference A minus B (95% CI) | p |",
            "|---|---|---|---|---|"]
    txt = "\n".join(head + rows) + "\n"
    open("contrasts.md", "w").write("# Method contrasts, post hoc (contrasts.py)\n\n" + __doc__.strip() + "\n\n" + txt)
    print(txt)


if __name__ == "__main__":
    main()
