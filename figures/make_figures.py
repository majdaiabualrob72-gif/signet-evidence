"""
Publication figures for the SIGNET manuscript (Research Synthesis Methods format).
Light background, greys plus one accent (blue for SIGNET; orange only where a second
method must be picked out); no tints. Figure 1 is a step-by-step guide to a SIGNET review
with the sigmoidoscopy question of Table 5 worked through each step; its numbers are read
from examples/sigmoidoscopy_distal.result.json when that file is present.
Run: python3 figures/make_figures.py [fig1 fig2 fig3 fig4 fig5]
Outputs: figures/out/AbuAlrobFig<k>.{pdf,png,tif}
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "figures", "out")
os.makedirs(OUT, exist_ok=True)

INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
BLUE, ORANGE = "#2a78d6", "#eb6834"
BAND = "#d9d8d3"
FONT = "Liberation Sans"   # Arial-compatible metrics
plt.rcParams.update({
    "font.family": FONT, "font.size": 8.5, "axes.edgecolor": INK2, "axes.linewidth": 0.6,
    "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "pdf.fonttype": 42, "ps.fonttype": 42})


def save(fig, k, dpi_tif=600):
    base = os.path.join(OUT, f"AbuAlrobFig{k}")
    # journal artwork rule: the figure number at the top of each figure file, placed just
    # above the drawn content (not at the canvas corner, which left an empty band)
    fig.canvas.draw()
    bb = fig.get_tightbbox(fig.canvas.get_renderer())      # inches
    W, H = fig.get_size_inches()
    fig.text(bb.x0 / W, (bb.y1 + 0.10) / H, f"Figure {k}", ha="left", va="bottom", fontsize=10,
             fontweight="bold", color="#1a1a1a", transform=fig.transFigure)
    fig.savefig(base + ".pdf", bbox_inches="tight")
    fig.savefig(base + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(base + ".tif", dpi=dpi_tif, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    # flatten to RGB on white: journal systems may reject TIFFs with an alpha channel
    from PIL import Image
    im = Image.open(base + ".tif"); dpi = im.info.get("dpi")
    if im.mode == "RGBA":
        bg = Image.new("RGB", im.size, (255, 255, 255)); bg.paste(im, mask=im.split()[-1])
        bg.save(base + ".tif", compression="tiff_lzw", dpi=dpi)
    print("saved", base)


# ----------------------------------------------------------------------
# Figure 1: how to run a SIGNET review, with a worked example at each step
# ----------------------------------------------------------------------
def _example_numbers():
    """Numbers for the worked example, read from the archived specification and CLI output."""
    ex = os.path.join(ROOT, "examples")
    spec = json.load(open(os.path.join(ex, "sigmoidoscopy_distal.json")))
    res = json.load(open(os.path.join(ex, "sigmoidoscopy_distal.result.json")))["targets"]["distal"]
    cells = {c["id"]: c for c in spec["cells"]}

    def rr(v):
        return "unbounded" if v in ("inf", "-inf", None) else f"{np.exp(float(v)):.2f}"

    def rset(s):
        lo = "0" if s["lo"] in ("-inf", None) else rr(s["lo"])
        return f"{lo} to {rr(s['hi'])}"
    d, p = cells["B_distal"], cells["B_proximal"]
    r1, r2 = res["rung1"], res["rung2"]
    return dict(
        distal=f"{d['estimate']:.2f} (95% CI {d['lower']:.2f} to {d['upper']:.2f})",
        proximal=f"{p['estimate']:.2f} ({p['lower']:.2f} to {p['upper']:.2f})",
        r1=rset(r1["set"]), r2=rset(r2["set"]), r2k1=rset(r2["tolerance_k1"]),
        p_ben=f"{r2['claims']['benefit']['p']:.3f}", p_harm=f"{r2['claims']['harm']['p']:.2f}",
        test_p=f"{r2['model_test']['p']:.2f}",
        default=("unbounded" if r2.get("default", {}).get("hi", "inf") == "inf" else rset(r2["default"])),
        cap1=r1["capacity"]["above"], cap2=r2["capacity"]["above"])


def fig1():
    from matplotlib.font_manager import FontProperties
    from matplotlib.patches import Circle, RegularPolygon, Rectangle, Arc
    from matplotlib.lines import Line2D
    E = _example_numbers()
    assert E["cap1"] is False and E["cap2"] is True, "capacity pattern of the worked example changed"
    W = 7.2
    X_SP, X_S0, X_S1, X_E0, X_E1 = 0.30, 0.52, 4.00, 4.16, 7.16
    PAD, GAP = 0.042, 0.046
    F_T, F_B, F_H = 7.5, 6.8, 8.2
    meas_fig = plt.figure(figsize=(1, 1), dpi=72)
    R = meas_fig.canvas.get_renderer()

    def tw(s, size, weight="normal", style="normal"):
        fp = FontProperties(family=FONT, size=size, weight=weight, style=style)
        return R.get_text_width_height_descent(s, fp, ismath=False)[0] / 72.0

    def wrap(text, width, size, weight="normal", style="normal"):
        out = []
        for para in text.split("\n"):
            line = ""
            for word in para.split(" "):
                cand = word if not line else line + " " + word
                if tw(cand, size, weight, style) <= width:
                    line = cand
                else:
                    out.append(line); line = word
            out.append(line)
        return out
    LH = lambda size: size * 1.24 / 72.0

    steps = [
        ("A", "Specification", "before any estimate is extracted"),
        (1, "Define the estimand", "Population, intervention, comparator, outcome and scale; an effect field and kernels if populations differ.",
         "Screening sigmoidoscopy versus none; death from distal colorectal cancer; log ratio scale.", "step"),
        (2, "List the evidence cells", "One per design, population and outcome, including negative and signed controls.",
         "Case-control cells: distal cancer death (measures the effect); proximal cancer death, beyond the scope's reach (kernel 0).", "step"),
        (3, "List the bias processes in each cell", "Sign each process with a documented direction; leave others unsigned.",
         "Healthy-screenee bias H in both cells, toward lower observed risk.", "step"),
        (4, "Encode named assumption groups", "Sign, zero, share, order, field-sign or monotone constraints, each with a rung and a source that predates the results.",
         "Rung 1 (sign): H lowers observed risk. Rung 2 (share): H is the same in the distal and proximal cells.", "step"),
        (5, "Check capacity for each claim", "Theorem 1, without data: can benefit or harm be certified at each rung?\nIf not: record the design or link that would make it certifiable.",
         "Rung 1: benefit uncertifiable, since a bias toward lower risk could explain any benefit. Rung 2: certifiable via the shared bias.", "check"),
        (6, "Compute the certificate weights", "From the design and the expected variances, never from the estimates.",
         "Benefit bound: distal minus proximal (weights +1 and −1), resting on the shared-bias group alone.", "step"),
        ("gate", "Register and hash", "the specification, capacity map and certificates before extraction.",
         "Protocol frozen and hashed before extraction (SHA-256 prefix d8c5c25b)."),
        ("B", "Analysis", "after extraction"),
        (7, "Extract and pool", "Convert estimates to the log scale; pool by random effects within each cell.",
         f"Distal OR {E['distal']}; proximal OR {E['proximal']}. One study per cell, so no pooling is needed.", "step"),
        (8, "Test the assumption set", "Conditional chi-squared test of all groups in force.\nIf rejected: list the minimal repairs.",
         f"Not rejected (p = {E['test_p']}): the data fit the assumptions exactly.", "check"),
        (9, "Compute SIGNET sets", "Each rung, tolerating 0, 1 or 2 failed groups or one broken cell.",
         f"Rung 1: {E['r1']}. Rung 2: {E['r2']}. Rung 2 tolerating one failed group: {E['r2k1']}.", "step"),
        (10, "Compute p-values for directional claims", "Benefit and harm: p = min(1, 4 p(certificate), 2 p(test)).",
         f"Benefit at rung 2: p = {E['p_ben']}; harm: p = {E['p_harm']}.", "step"),
        ("C", "Reporting", ""),
        (11, "Report the default as the primary result", "The set at the highest rung specified, tolerating one failed group or one broken cell, with its capacity; only this set supports a directional claim. Beside it, report the same rung without tolerance.",
         f"Default at rung 2: {E['default']}. It had no capacity for a claim: one cell measures the effect, and omitting it leaves nothing to bound the effect. Beside it, rung 2 without tolerance: {E['r2']}, conditional on the shared-bias assumption.", "primary"),
        (12, "Report what carries each bound", "Frontier of rungs and k, certificates, model test, refutation distance.",
         "Benefit at rung 2 rests on one named group, the shared bias: refutation distance 1.", "step"),
        (13, "Report structural gaps", "Each uncertifiable claim and the design or link that would close it.",
         "Gap: one design measures the effect. A trial would close it; PLCO later reported RR 0.50 (0.38 to 0.64).", "step"),
    ]
    sw = X_S1 - X_S0 - 2 * 0.09            # text width in step boxes
    ew = X_E1 - X_E0 - 0.12 - 0.04          # text width in example blocks
    # ---- layout pass
    rows = []
    for item in steps:
        if item[0] in ("A", "B", "C"):
            rows.append(dict(kind="stage", h=0.23, item=item)); continue
        if item[0] == "gate":
            tl = wrap(item[1] + " " + item[2], sw - 0.25, F_B)
            el = wrap(item[3], ew, F_B)
            h = 2 * PAD + max(len(tl), len(el)) * LH(F_B) + 0.04
            rows.append(dict(kind="gate", h=h, tl=tl, el=el, item=item)); continue
        k, title, body, ex, kind = item
        tl = wrap(title, sw, F_T, "bold")
        bl = wrap(body, sw, F_B)
        el = wrap(ex, ew, F_B)
        hs = 2 * PAD + len(tl) * LH(F_T) + len(bl) * LH(F_B)
        he = 2 * PAD + len(el) * LH(F_B)
        rows.append(dict(kind=kind, h=max(hs, he, 0.30), tl=tl, bl=bl, el=el, item=item))
    plt.close(meas_fig)
    sub = ["Screening sigmoidoscopy and distal colorectal cancer death,", "evidence of 1992"]
    head_h, key_h = 0.15 + LH(F_H) + len(sub) * LH(6.6), 0.30
    H = head_h + sum(r["h"] for r in rows) + GAP * (len(rows) - 1) + key_h + 0.10
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
    ax.set_aspect("equal")

    def text(x, y, s, size, weight="normal", style="normal", color=INK, **kw):
        ax.text(x, y, s, fontsize=size, fontweight=weight, fontstyle=style, color=color, va="top", ha="left", **kw)

    def box(x0, y0, x1, y1, ec=INK2, lw=0.8):
        ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle="round,pad=0,rounding_size=0.045",
                                    linewidth=lw, edgecolor=ec, facecolor="white", zorder=2))
    # ---- column headers
    y = H - 0.06
    text(X_S0, y, "Steps of a SIGNET review", F_H, "bold")
    text(X_E0, y, "Worked example (Table 5)", F_H, "bold")
    for i, ln in enumerate(sub):
        text(X_E0, y - LH(F_H) - 0.01 - i * LH(6.6), ln, 6.6, style="italic", color=INK2)
    y -= head_h
    ax.plot([X_S0, X_E1], [y + 0.07, y + 0.07], color=INK, lw=0.8, solid_capstyle="butt")
    markers = []
    for r in rows:
        top, bot = y, y - r["h"]
        it = r["item"]
        if r["kind"] == "stage":
            ax.plot([X_S0, X_E1], [top + 0.005, top + 0.005], color=MUTED, lw=0.5) if it[0] != "A" else None
            text(X_S0, top - 0.06, f"{it[0]}. {it[1]}", F_H, "bold")
            if it[2]:
                text(X_S0 + tw(f"{it[0]}. {it[1]}", F_H, "bold") + 0.10, top - 0.072, it[2], 7.0, style="italic", color=INK2)
        elif r["kind"] == "gate":
            box(X_S0, bot, X_E1, top, ec=BLUE, lw=1.2)
            yy = top - PAD - 0.02
            for i, ln in enumerate(r["tl"]):
                if i == 0:
                    lead = "Register and hash"
                    text(X_S0 + 0.09, yy, lead, F_B, "bold")
                    text(X_S0 + 0.09 + tw(lead + " ", F_B, "bold"), yy, ln[len(lead) + 1:], F_B)
                else:
                    text(X_S0 + 0.09, yy, ln, F_B)
                yy -= LH(F_B)
            yy = top - PAD - 0.02
            for ln in r["el"]:
                text(X_E0 + 0.12, yy, ln, F_B, color=INK); yy -= LH(F_B)
            ax.plot([X_E0, X_E0], [bot + 0.05, top - 0.05], color=MUTED, lw=0.5)
            markers.append(("lock", (top + bot) / 2))
        else:
            k = it[0]
            primary = r["kind"] == "primary"
            box(X_S0, bot, X_S1, top, ec=BLUE if primary else INK2, lw=1.6 if primary else 0.8)
            yy = top - PAD
            for ln in r["tl"]:
                text(X_S0 + 0.09, yy, ln, F_T, "bold"); yy -= LH(F_T)
            for ln in r["bl"]:
                bold_lead = ln.startswith("If not:") or ln.startswith("If rejected:")
                if bold_lead:
                    lead = ln.split(":")[0] + ":"
                    text(X_S0 + 0.09, yy, lead, F_B, "bold", color=INK)
                    text(X_S0 + 0.09 + tw(lead + " ", F_B, "bold"), yy, ln[len(lead) + 1:], F_B, color=INK2)
                else:
                    text(X_S0 + 0.09, yy, ln, F_B, color=INK2)
                yy -= LH(F_B)
            if primary:
                lab = "primary result"
                text(X_S1 - 0.08 - tw(lab, 6.4, "bold"), top - PAD - 0.01, lab, 6.4, "bold", color=BLUE)
            # example block: blue rule at left, text in ink
            ax.add_patch(Rectangle((X_E0, bot + 0.02), 0.028, r["h"] - 0.04, facecolor=BLUE, edgecolor="none", zorder=2))
            yy = top - PAD
            for ln in r["el"]:
                text(X_E0 + 0.12, yy, ln, F_B, color=INK); yy -= LH(F_B)
            markers.append((r["kind"], top - PAD - LH(F_T) / 2 + 0.01, k))
        y = bot - GAP
    # ---- spine with arrows between markers, then markers on top
    RAD = 0.105
    half = {"lock": 0.10, "check": 0.155}
    for m0, m1 in zip(markers[:-1], markers[1:]):
        a = m0[1] - half.get(m0[0], RAD) - 0.012
        b = m1[1] + half.get(m1[0], RAD) + 0.022
        ax.add_patch(FancyArrowPatch((X_SP, a), (X_SP, b), arrowstyle="-|>", mutation_scale=7,
                                     color=INK2, linewidth=0.9, zorder=1, shrinkA=0, shrinkB=0))
    for m in markers:
        if m[0] == "lock":
            yc = m[1]
            ax.add_patch(Rectangle((X_SP - 0.085, yc - 0.085), 0.17, 0.12, facecolor=BLUE, edgecolor=BLUE, zorder=4))
            ax.add_patch(Arc((X_SP, yc + 0.035), 0.11, 0.14, theta1=0, theta2=180, color=BLUE, linewidth=1.5, zorder=4))
            ax.add_patch(Rectangle((X_SP - 0.014, yc - 0.05), 0.028, 0.045, facecolor="white", edgecolor="none", zorder=5))
            continue
        kind, yc, k = m
        if kind == "check":
            ax.add_patch(RegularPolygon((X_SP, yc), numVertices=4, radius=0.15, orientation=0, facecolor=BLUE,
                                        edgecolor=BLUE, zorder=4))
        else:
            ax.add_patch(Circle((X_SP, yc), RAD, facecolor=INK, edgecolor=INK, zorder=4))
        ax.text(X_SP, yc - 0.004, str(k), fontsize=6.9 if k < 10 else 6.3, fontweight="bold", color="white",
                ha="center", va="center", zorder=5)
    # ---- key
    yk = y + GAP - 0.17
    xk = X_S0
    ax.add_patch(Circle((xk + 0.06, yk), 0.06, facecolor=INK, edgecolor=INK)); xk += 0.16
    ax.text(xk, yk, "step", fontsize=6.8, color=INK2, va="center", ha="left"); xk += tw("step", 6.8) + 0.25
    ax.add_patch(RegularPolygon((xk + 0.08, yk), numVertices=4, radius=0.085, orientation=0, facecolor=BLUE, edgecolor=BLUE)); xk += 0.24
    ax.text(xk, yk, "check with a branch", fontsize=6.8, color=INK2, va="center", ha="left"); xk += tw("check with a branch", 6.8) + 0.25
    ax.add_patch(Rectangle((xk, yk - 0.05), 0.12, 0.085, facecolor=BLUE, edgecolor=BLUE))
    ax.add_patch(Arc((xk + 0.06, yk + 0.035), 0.08, 0.10, theta1=0, theta2=180, color=BLUE, linewidth=1.1)); xk += 0.20
    ax.text(xk, yk, "registration point", fontsize=6.8, color=INK2, va="center", ha="left"); xk += tw("registration point", 6.8) + 0.25
    ax.add_patch(FancyBboxPatch((xk, yk - 0.06), 0.24, 0.12, boxstyle="round,pad=0,rounding_size=0.03", facecolor="white",
                                edgecolor=BLUE, linewidth=1.4)); xk += 0.32
    ax.text(xk, yk, "primary result", fontsize=6.8, color=INK2, va="center", ha="left"); xk += tw("primary result", 6.8) + 0.25
    ax.add_patch(Rectangle((xk, yk - 0.07), 0.028, 0.14, facecolor=BLUE, edgecolor="none")); xk += 0.10
    ax.text(xk, yk, "worked example", fontsize=6.8, color=INK2, va="center", ha="left")
    save(fig, 1, dpi_tif=1000)


# ----------------------------------------------------------------------
# Figure 5: real-data case series (function fig4)
# ----------------------------------------------------------------------
def _get(res, case):
    return [r for r in res if r["case"] == case][0]


def fig4():
    res = json.load(open(os.path.join(ROOT, "experiments", "realdata", "results_realdata.json")))
    aud = json.load(open(os.path.join(ROOT, "experiments", "audit_v12", "realdata_audit.json")))   # version 1.2.0 analyses
    rows_all = []
    A, B, C = _get(res, "A"), _get(res, "B"), _get(res, "C")
    z = 1.959964

    def ci(c):
        return c["m"] - z * c["se"], c["m"], c["m"] + z * c["se"]

    def sset(s):
        lo = -np.inf if s["lo"] in (None, "-inf") else float(s["lo"])
        hi = np.inf if s["hi"] in (None, "inf") else float(s["hi"])
        return lo, None, hi

    def k1(s):
        k = s.get("k1")
        lo = -np.inf if k["lo"] in (None, "-inf") else float(k["lo"])
        hi = np.inf if k["hi"] in (None, "inf") else float(k["hi"])
        return lo, None, hi

    def iv(lo, mid, hi):
        return float(lo), (None if mid is None else float(mid)), float(hi)

    NOCAP = (-np.inf, None, np.inf)          # the default: unbounded on both sides in every question
    panels = []
    ca, sa = A["comparators"]["cvd"], A["signet"]["cvd"]
    panels.append(("Hormone therapy and cardiovascular disease (evidence of 1994)",
                   [("Published best estimate", ci(ca["conventional"]), "c", "contradicted"),
                    ("Bias-adjusted, unsigned (BA-U)", ci(ca["bayes_unsigned"]), "c", "contradicted"),
                    ("Bias-adjusted, signed (BA-S)", ci(ca["bayes_signed"]), "c", "contradicted"),
                    ("Negative-control calibration", ci(ca["nco_calibration"]), "c", "contradicted"),
                    ("SIGNET rung 1", sset(sa["rung1"]), "s", "no claim"),
                    ("SIGNET rung 2 (signed control)", sset(sa["rung2"]), "s", "no claim"),
                    ("SIGNET rung 3 (null control)", sset(sa["rung3"]), "s", "contradicted"),
                    ("SIGNET rung 3, one failure", k1(sa["rung3"]), "s", "no claim"),
                    ("SIGNET default", NOCAP, "s", "no capacity")],
                   [("WHI E+P", 1.22, 1.09, 1.36), ("WHI E alone", 1.12, 1.01, 1.24)]))
    cb, sb = B["comparators"]["distal"], B["signet"]["distal"]
    panels.append(("Screening sigmoidoscopy and distal colorectal cancer death (evidence of 1992)",
                   [("Case-control estimate", ci(cb["conventional"]), "c", "confirmed"),
                    ("Bias-adjusted, unsigned (BA-U)", ci(cb["bayes_unsigned"]), "c", "confirmed"),
                    ("Bias-adjusted, signed (BA-S)", ci(cb["bayes_signed"]), "c", "confirmed"),
                    ("Negative-control calibration", ci(cb["nco_calibration"]), "c", "confirmed"),
                    ("SIGNET rung 1", sset(sb["rung1"]), "s", "missed"),
                    ("SIGNET rung 2 (anatomical control)", sset(sb["rung2"]), "s", "confirmed"),
                    ("SIGNET rung 2, one failure", k1(sb["rung2"]), "s", "missed"),
                    ("SIGNET default", NOCAP, "s", "no capacity")],
                   [("PLCO", 0.50, 0.38, 0.64)]))
    cc, sc = C["comparators"]["chd"], C["signet"]["chd"]
    pooled = C["pooled_reference"]
    rv, bs = aud["re_variants"]["vitamin_E"], aud["ba_s"]["vitamin_E"]["shared"]
    # the Wald row must be the frozen analysis: check that the audit reproduces it
    assert abs(rv["est"] - cc["conventional"]["m"]) < 1e-9 and abs(rv["wald"][0] - (cc["conventional"]["m"] - z * cc["conventional"]["se"])) < 1e-6
    panels.append(("Vitamin E supplements and coronary disease (evidence of 1993)",
                   [("Random-effects pool, Wald interval", ci(cc["conventional"]), "c", "not confirmed"),
                    ("Random-effects pool, Hartung-Knapp", iv(rv["hk"][0], rv["est"], rv["hk"][1]), "c", "not confirmed"),
                    ("Random-effects pool, modified Hartung-Knapp", iv(rv["mhk"][0], rv["est"], rv["mhk"][1]), "c", "no claim"),
                    ("Bias-adjusted, unsigned (BA-U)", ci(cc["bayes_unsigned"]), "c", "not confirmed"),
                    ("BA-S, independent bias terms", ci(cc["bayes_signed"]), "c", "not confirmed"),
                    ("BA-S, shared bias term", iv(bs["lo"], bs["m"], bs["hi"]), "c", "no claim"),
                    ("SIGNET rung 1", sset(sc["rung1"]), "s", "concordant null"),
                    ("SIGNET default", NOCAP, "s", "no capacity")],
                   [("WHS and PHS II (pooled MI)", pooled["point"], pooled["lower"], pooled["upper"])]))
    heights = [len(p[1]) + 1.2 for p in panels]
    fig, axes = plt.subplots(len(panels), 1, figsize=(6.3, 0.30 * sum(heights) + 0.9),
                             gridspec_kw=dict(height_ratios=heights, hspace=0.36), sharex=True)
    xmin, xmax = np.log(0.12), np.log(4.0)
    for ax, (title, rows, refs) in zip(axes, panels):
        n = len(rows)
        for j, (lab, p, lo, hi) in enumerate(refs):
            # no tints (journal artwork rule): trial estimate solid, its 95% CI dotted
            for edge in (lo, hi):
                ax.axvline(np.log(edge), color=INK2, lw=0.7, ls=(0, (1, 1.6)), zorder=1)
            ax.axvline(np.log(p), color=INK2, lw=1.0, ls="-", zorder=1)
            ax.text(np.log(p), n - 0.35 + 0.0 * j, lab, fontsize=7, color=INK2,
                    ha="left" if j == 0 else "right", va="bottom",
                    transform=ax.transData, zorder=5,
                    bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
        ax.axvline(0, color=MUTED, lw=0.7, ls=(0, (5, 3)), zorder=1)
        for i, (lab, (lo, mid, hi), kind, outcome) in enumerate(rows):
            y = n - 1 - i
            col = BLUE if kind == "s" else INK2
            l_ = max(lo, xmin); h_ = min(hi, xmax)
            ax.plot([l_, h_], [y, y], color=col, lw=2.0 if kind == "s" else 1.4, solid_capstyle="butt", zorder=3)
            if not np.isfinite(hi) or hi > xmax:
                ax.annotate("", xy=(xmax, y), xytext=(xmax - 0.25, y),
                            arrowprops=dict(arrowstyle="-|>", color=col, lw=1.4, mutation_scale=8), zorder=4)
            if not np.isfinite(lo) or lo < xmin:
                ax.annotate("", xy=(xmin, y), xytext=(xmin + 0.25, y),
                            arrowprops=dict(arrowstyle="-|>", color=col, lw=1.4, mutation_scale=8), zorder=4)
            if mid is not None:
                ax.plot([mid], [y], "o", ms=4, color=col, zorder=4)
            ax.text(xmax + 0.12, y, outcome, fontsize=7.2, color=INK2 if outcome not in ("contradicted",) else INK,
                    fontweight="bold" if outcome == "contradicted" else "normal", va="center")
        ax.set_yticks(range(n))
        ax.set_yticklabels([r[0] for r in rows][::-1], fontsize=7.6)
        ax.set_ylim(-0.7, n + 0.25)
        ax.set_xlim(xmin, xmax)
        ax.set_title(title, fontsize=8.6, loc="left", color=INK, pad=4)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="x", color=GRID, lw=0.5, zorder=0)
    ticks = [0.125, 0.25, 0.5, 1, 2, 4]
    axes[-1].set_xticks(np.log(ticks))
    axes[-1].set_xticklabels([str(t) for t in ticks])
    axes[-1].set_xlabel("Ratio (log scale). Solid line: trial estimate; dotted lines: its 95% CI;\ndashed line: no effect; arrows: unbounded limits", fontsize=7.8)
    axes[0].text(xmax + 0.12, len(panels[0][1]) - 0.35, "Scored against\nthe trial", fontsize=7, color=MUTED, va="bottom")
    save(fig, 5)   # real-data figure: Figure 5 (cited after the redesign figure)



# ----------------------------------------------------------------------
# Figure 2: coverage and width across 15 scenarios
# ----------------------------------------------------------------------
def fig2():
    sim = os.path.join(ROOT, "experiments", "simulation")
    allres = json.load(open(os.path.join(sim, "all_scenarios_summary.json")))
    order = ["no bias", "reference", "null effect", "dense small biases", "magnitudes doubled", "15% of signs wrong",
             "omitted bias", "no design opposing benefit", "one-signed literature, null effect",
             "one-signed literature, null effect, magnitudes doubled", "extreme outlier", "publication bias",
             "heterogeneous populations", "missing controls and variances", "contradictory evidence"]
    labels = {"one-signed literature, null effect": "one-signed, null effect",
              "one-signed literature, null effect, magnitudes doubled": "one-signed, null, doubled",
              "no design opposing benefit": "no design opposing benefit", "missing controls and variances": "missing controls, variances",
              "heterogeneous populations": "heterogeneous populations",
              "15% of signs wrong": "15% of signs reversed", "magnitudes doubled": "biases doubled"}
    meths = [("RE (all designs)", "Random-effects", INK2, "o", True),
             ("BA-S (elicited directions and magnitudes)", "Bias-adjusted, elicited magnitudes (BA-S)", ORANGE, "D", True),
             ("SIGNET R3", "SIGNET rung 3", BLUE, "o", False),
             ("SIGNET R3 k=1", "SIGNET rung 3, tolerant to one failed group", BLUE, "o", True)]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 5.6), sharey=True, gridspec_kw=dict(width_ratios=[1.35, 1], wspace=0.13))
    ys = np.arange(len(order))[::-1]
    for ax in axes:
        for y in ys:
            ax.axhline(y, color=GRID, lw=0.5, zorder=0)
    axes[0].axvline(0.95, color=INK2, lw=0.8, ls=(0, (3, 2)))
    axes[0].text(0.95, len(order) - 0.3, "nominal 0.95", fontsize=7, color=INK2, ha="right", va="bottom")
    offs = np.linspace(-0.27, 0.27, len(meths))
    for (key, lab, col, mk, filled), off in zip(meths, offs):
        cov, wid = [], []
        for sc in order:
            d = allres[sc].get(key, {})
            cov.append(d.get("coverage", np.nan))
            w = d.get("median_width")
            wid.append(np.nan if w is None else w)
        kw = dict(color=col, marker=mk, ms=4.5, ls="none", label=lab,
                  markerfacecolor=col if filled else "white", markeredgewidth=1.0)
        axes[0].plot(cov, ys + off, **kw)
        axes[1].plot(wid, ys + off, **kw)
        if key == "SIGNET R3 k=1":
            for y, sc in zip(ys, order):
                d = allres[sc].get(key, {})
                if d.get("unbounded", 0) and d.get("unbounded") >= 0.999:
                    axes[1].text(0.02, y + 0.18, "SIGNET sets unbounded", fontsize=6.6, color=BLUE, va="center")
    axes[0].set_yticks(ys)
    axes[0].set_yticklabels([labels.get(s, s) for s in order], fontsize=7.6)
    axes[0].set_xlim(0.15, 1.02)
    axes[0].set_xlabel("Coverage of the true effect (400 literatures each)")
    axes[1].set_xlim(0, 0.75)
    axes[1].set_xlabel("Median width of bounded sets (log scale)")
    for ax in axes:
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.tick_params(axis="y", length=0)
    for ax in axes:
        ax.axhline(4.5, color=MUTED, lw=0.6)
    axes[0].text(0.16, 4.38, "stress scenarios", fontsize=6.8, color=MUTED, va="top", ha="left")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, frameon=False, fontsize=7.4, bbox_to_anchor=(0.5, -0.04))
    save(fig, 2)


# ----------------------------------------------------------------------
# Figure 3: magnitude sweep
# ----------------------------------------------------------------------
def fig3():
    sim = os.path.join(ROOT, "experiments", "simulation")
    sw = json.load(open(os.path.join(sim, "benchmark_v2_sweep.json")))["results"]
    mults = [0.5, 1.0, 2.0, 3.0, 4.0]
    fams = [("one-signed literature, null effect", "benefit", "A. One-signed literature, null effect:\nshare of literatures claiming benefit"),
            ("mixed designs, true benefit", "coverage", "B. Mixed designs, true benefit:\ncoverage of the true effect")]
    meths = [("RE (all designs)", "Random-effects", INK2, "-", "o"),
             ("BA-U (unsigned bias priors)", "BA-U", MUTED, (0, (3, 2)), "s"),
             ("BA-S (elicited directions and magnitudes)", "BA-S", ORANGE, "-", "D"),
             ("SIGNET-B k=1", "SIGNET with bounds at 3\u00d7, tolerant", BLUE, (0, (3, 2)), "^"),
             ("SIGNET R3 k=1", "SIGNET, tolerant to one failure", BLUE, "-", "o")]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw=dict(wspace=0.3))
    for ax, (fam, metric, title) in zip(axes, fams):
        for key, lab, col, ls, mk in meths:
            vals = [sw[f"{fam} | x{m}"][key][metric] for m in mults]
            ax.plot(mults, vals, color=col, ls=ls, marker=mk, ms=4, lw=1.4, label=lab,
                    markerfacecolor=col if mk != "^" else "white")
        ax.set_xticks(mults)
        ax.set_xticklabels(["\u00d70.5", "\u00d71", "\u00d72", "\u00d73", "\u00d74"])
        ax.set_title(title, fontsize=8.2, loc="left", color=INK)
        ax.set_ylim(-0.03, 1.03)
        ax.grid(axis="y", color=GRID, lw=0.5)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[1].axhline(0.95, color=INK2, lw=0.8, ls=(0, (1.5, 1.5)))
    fig.supxlabel("True bias sizes relative to the elicited ones", fontsize=8.5, y=0.0)
    h, l = axes[0].get_legend_handles_labels()
    order = [0, 1, 2, 4, 3]       # established methods in the left column, SIGNET in the right
    fig.legend([h[i] for i in order], [l[i] for i in order], loc="upper center", ncol=2, frameon=False,
               fontsize=7.4, bbox_to_anchor=(0.5, -0.04), columnspacing=3.0)
    save(fig, 3)


# ----------------------------------------------------------------------
# Figure 4: redesigns evaluated after the stress tests (function fig5)
# ----------------------------------------------------------------------
def fig5():
    sim = os.path.join(ROOT, "experiments", "simulation")
    base = json.load(open(os.path.join(sim, "all_scenarios_summary.json")))
    het = json.load(open(os.path.join(sim, "benchmark_v2_hetallow.json")))["results"]
    bnd = json.load(open(os.path.join(sim, "benchmark_v2_bounded.json")))["results"]
    cr = json.load(open(os.path.join(sim, "benchmark_v2_cellrobust.json")))["results"]
    chp = os.path.join(sim, "benchmark_v2_cellhet.json")
    ch = json.load(open(chp))["results"] if os.path.exists(chp) else {}
    src = {"base": base, "het": het, "bnd": bnd, "cr": cr, "ch": ch}
    # label, source, key, colour, marker, filled
    rows = [("Random-effects", "base", "RE (all designs)", INK2, "o", True),
            ("BA-S", "base", "BA-S (elicited directions and magnitudes)", ORANGE, "D", True),
            ("SIGNET rung 3", "base", "SIGNET R3", BLUE, "o", False),
            ("Rung 3 + heterogeneity allowance", "het", "SIGNET R3 + heterogeneity allowance", BLUE, "s", False),
            ("Tolerant to one failed group", "base", "SIGNET R3 k=1", BLUE, "o", True),
            ("Tolerant + heterogeneity allowance", "het", "SIGNET R3 k=1 + heterogeneity allowance", BLUE, "s", True),
            ("Tolerant to one failed group or cell", "cr", "SIGNET R3 k=1 (group or cell)", BLUE, "^", True),
            ("Group or cell + heterogeneity allowance", "ch", "SIGNET R3 k=1 (group or cell) + heterogeneity allowance", BLUE, "v", True),
            ("Magnitude rung (3\u00d7), tolerant", "bnd", "SIGNET-B k=1", BLUE, "^", False)]
    panels = [("reference", "coverage", "A. Reference"),
              ("heterogeneous populations", "coverage", "B. Heterogeneous\npopulations"),
              ("contradictory evidence", "coverage", "C. Contradictory\nevidence"),
              ("reference", "benefit", "D. Reference")]
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 3.2), sharey=True,
                             gridspec_kw=dict(wspace=0.12, width_ratios=[1, 1, 1, 1]))
    ys = np.arange(len(rows))[::-1]
    for ax, (sc, metric, title) in zip(axes, panels):
        for y in ys:
            ax.axhline(y, color=GRID, lw=0.5, zorder=0)
        for y, (lab, s, key, col, mk, filled) in zip(ys, rows):
            d = src[s].get(sc, {}).get(key)
            if not isinstance(d, dict) or d.get(metric) is None:
                ax.text(0.33 if metric == "coverage" else 0.0, y, "not run", fontsize=6.4, color=MUTED,
                        ha="left", va="center")
                continue
            ax.plot([d[metric]], [y], color=col, marker=mk, ms=5, ls="none",
                    markerfacecolor=col if filled else "white", markeredgewidth=1.0, zorder=3)
        if metric == "coverage":
            ax.set_xlim(0.30, 1.03)
            ax.axvline(0.95, color=INK2, lw=0.8, ls=(0, (3, 2)))
            ax.set_xticks([0.4, 0.6, 0.8, 1.0])
            ax.set_xlabel("Coverage", fontsize=7.6)
        else:
            ax.set_xlim(-0.03, 1.03)
            ax.set_xticks([0, 0.5, 1.0])
            ax.set_xlabel("Certified benefit", fontsize=7.6)
        ax.set_title(title, fontsize=7.6, loc="left", color=INK)
        for s_ in ("top", "right", "left"):
            ax.spines[s_].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.tick_params(axis="x", labelsize=7)
    axes[0].set_yticks(ys)
    axes[0].set_yticklabels([r[0] for r in rows], fontsize=7.2)
    save(fig, 4)   # redesign figure: Figure 4 (cited first, in Section 7.4)


if __name__ == "__main__":
    which = sys.argv[1:] or ["fig1", "fig2", "fig3", "fig4", "fig5"]
    for w in which:
        globals()[w]()
