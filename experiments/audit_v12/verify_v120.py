"""Quick verification of signet 1.2.0 against the archived outputs of the forensic-audit revision.

Run from the repository root:   python3 experiments/audit_v12/verify_v120.py
It recomputes a sample of every new analysis and compares it with the archived output. Nothing in the archive is
overwritten: scripts that write next to themselves are run on copies of their outputs, which are then restored.
Time: about 15 minutes on two cores. Exit code 0 means that every comparison passed."""
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
import sim_audit as S                      # noqa: E402

FAIL = []


def report(name, ok, detail):
    print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}", flush=True)
    if not ok:
        FAIL.append(name)


def walk(a, b, path="", out=None, rtol=1e-6, atol=1e-9, skip=("signet_version", "seconds")):
    """Collect the paths at which two JSON objects differ."""
    out = [] if out is None else out
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k in skip:
                continue
            if k not in a or k not in b:
                out.append(f"{path}/{k}: missing on one side")
            else:
                walk(a[k], b[k], f"{path}/{k}", out, rtol, atol, skip)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append(f"{path}: lengths {len(a)} and {len(b)}")
        else:
            for i, (x, y) in enumerate(zip(a, b)):
                walk(x, y, f"{path}[{i}]", out, rtol, atol, skip)
    elif isinstance(a, bool) or isinstance(b, bool) or a is None or b is None or isinstance(a, str) or isinstance(b, str):
        if a != b:
            out.append(f"{path}: {a!r} against {b!r}")
    else:
        x, y = float(a), float(b)
        if math.isnan(x) and math.isnan(y):
            return out
        if math.isinf(x) or math.isinf(y):
            if x != y:
                out.append(f"{path}: {x} against {y}")
        elif abs(x - y) > atol + rtol * max(abs(x), abs(y)):
            out.append(f"{path}: {x!r} against {y!r}")
    return out


def rerun_and_compare(name, script, outputs, rtol=1e-6):
    """Run a script that writes its outputs next to itself, compare the JSON output with the archived copy, restore."""
    tmp = tempfile.mkdtemp()
    try:
        for f in outputs:
            shutil.copy2(f, os.path.join(tmp, os.path.basename(f)))
        env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
        r = subprocess.run([sys.executable, script], cwd=ROOT, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        if r.returncode != 0:
            report(name, False, "script failed: " + r.stderr.decode()[-300:])
            return
        j = [f for f in outputs if f.endswith(".json")][0]
        new = json.load(open(j)); old = json.load(open(os.path.join(tmp, os.path.basename(j))))
        d = walk(old, new, rtol=rtol)
        report(name, not d, f"{len(d)} differences from the archived output" + (": " + "; ".join(d[:3]) if d else ""))
    finally:
        for f in outputs:
            shutil.copy2(os.path.join(tmp, os.path.basename(f)), f)
        shutil.rmtree(tmp)


def same_record(a, b, tol=1e-6):
    for k in ("empty", "cover", "benefit", "harm", "bounded"):
        if a.get(k) != b.get(k):
            return False
    for k in ("lo", "hi"):
        x, y = a.get(k), b.get(k)
        if (x is None) != (y is None):
            return False
        if x is not None and abs(x - y) > tol:
            return False
    return True


def sample(kind, fun, per, start=0):
    n_ok = n_tot = 0
    for sc in (S.HELD_ALL if kind == "heldout" else S.SCEN15):
        fn = S.fname(kind, sc)
        if not os.path.exists(fn):
            report(f"{kind} | {sc}", False, "archived output missing")
            continue
        arch = json.load(open(fn))["records"]
        for r in range(start, start + per):
            new = fun(sc, r)
            for k, v in new.items():
                if isinstance(v, dict) and "cover" in v:
                    n_tot += 1
                    n_ok += same_record(v, arch[r][k])
    report(f"{kind}: sample of {per} literatures per scenario", n_ok == n_tot, f"{n_ok} of {n_tot} sets identical to the archive (limits within 1e-6)")


def main():
    import signet
    print("signet", signet.__version__, flush=True)
    rd = os.path.join(ROOT, "experiments", "realdata")
    rerun_and_compare("real-data checks of version 1.1.0 under this version", os.path.join(rd, "v11_checks.py"),
                      [os.path.join(rd, "v11_checks.json"), os.path.join(rd, "v11_checks.md")])
    rerun_and_compare("real-data analyses of the forensic-audit revision", os.path.join(HERE, "realdata_audit.py"),
                      [os.path.join(HERE, "realdata_audit.json"), os.path.join(HERE, "realdata_audit.md")])
    # comparators: first 25 literatures of each scenario
    arch = json.load(open(os.path.join(S.OUT, "comparators__15_scenarios.json")))["records"]
    n_ok = n_tot = 0
    for sc in S.SCEN15:
        for r in range(25):
            new = S.comparators_one(sc, S.seed_of(sc, r))
            for k, v in new.items():
                if isinstance(v, dict) and "cover" in v:
                    n_tot += 1
                    n_ok += same_record(v, arch[sc][r][k], tol=1e-9)
    report("comparators (random effects with three intervals, BA-S independent and shared): 25 literatures per scenario",
           n_ok == n_tot, f"{n_ok} of {n_tot} intervals identical to the archive")
    sample("ablation", S.job_ablation, 4)
    sample("groups", S.job_groups, 3)
    sample("heldout", S.job_heldout, 3)
    print("ALL PASS" if not FAIL else f"FAILED: {len(FAIL)} comparison(s): " + "; ".join(FAIL))
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
