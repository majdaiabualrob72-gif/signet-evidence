"""Command line interface.

    python -m signet validate spec.json
    python -m signet run spec.json --out results.json --report report.md
"""
import argparse
import json
import sys

from .spec import validate, run, report_md, SpecError


def main(argv=None):
    ap = argparse.ArgumentParser(prog="signet", description="Signed evidence triangulation")
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate", help="check a specification")
    v.add_argument("spec")
    r = sub.add_parser("run", help="run a specification")
    r.add_argument("spec")
    r.add_argument("--out", default=None)
    r.add_argument("--report", default=None)
    a = ap.parse_args(argv)
    spec = json.load(open(a.spec))
    try:
        validate(spec)
    except SpecError as e:
        print(f"invalid specification: {e}", file=sys.stderr)
        return 2
    if a.cmd == "validate":
        print("specification is valid")
        return 0
    res = run(spec)
    txt = json.dumps(res, indent=1)
    if a.out:
        open(a.out, "w").write(txt)
    else:
        print(txt)
    if a.report:
        open(a.report, "w").write(report_md(res))
    return 0


if __name__ == "__main__":
    sys.exit(main())
