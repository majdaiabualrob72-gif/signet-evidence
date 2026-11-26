#!/bin/sh
# Re-run the version 1.2.0 checks on your own computer and write a log.
# Usage, from the repository root:   sh scripts/reproduce_v120.sh
# Time: about 15 minutes (14.5 minutes in the development environment on 22 November 2026).
# Requires: Python 3.13 with requirements-lock.txt installed and `pip install -e .`
# Optional: R with the packages jsonlite, lpSolve and quadprog (for the R interface and the R reconstruction).
#
# What it covers: every unit check (versions 1.0.0, 1.1.0 and 1.2.0), the real-data analyses of the
# forensic-audit revision, a sample of every new simulation analysis compared with the archived output,
# the R interface and the R reconstruction of the worked example.
# What it does not cover: the full new simulation runs (about 4 core-hours), which are:
#   python3 experiments/audit_v12/sim_audit.py comparators
#   python3 experiments/audit_v12/sim_audit.py worker        (start it twice to use two cores)
#   python3 experiments/audit_v12/analyse_audit.py
#   python3 experiments/audit_v12/certificate_weights.py 400
#   python3 experiments/audit_v12/tolerance_sensitivity.py
PY=${PY:-python3}
LOG=REPRODUCTION_v1.2.0_$(date +%Y%m%d_%H%M).log
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
{
  echo "== SIGNET v1.2.0 reproduction, $(date)"
  echo "== machine: $(uname -a)"
  $PY -c "import sys, numpy, scipy, clarabel; import signet; print('python', sys.version.split()[0], 'numpy', numpy.__version__, 'scipy', scipy.__version__, 'clarabel', clarabel.__version__, 'signet', signet.__version__)"
  for t in tests/test_core.py tests/test_spec.py tests/test_cell_tolerance.py tests/test_v11.py tests/test_v12.py; do
    echo "== $t"; $PY $t; echo "exit code $?"
  done
  echo "== the four example specifications"
  for x in sigmoidoscopy_distal hormone_therapy_cvd vitamin_e_chd beta_carotene_chd_trial_validity_named; do
    $PY -m signet run examples/$x.json --out /tmp/signet_$x.json --report /tmp/signet_$x.md > /dev/null; echo "$x exit code $?"
    $PY - "$x" << 'PYEOF'
import json, sys
x = sys.argv[1]
a = json.load(open(f"examples/{x}.result.json")); b = json.load(open(f"/tmp/signet_{x}.json"))
def walk(u, v, out):
    if isinstance(u, dict):
        for k in u:
            if k in ("signet_version",): continue
            walk(u[k], v.get(k) if isinstance(v, dict) else None, out)
    elif isinstance(u, list):
        for i, w in enumerate(u):
            walk(w, v[i] if isinstance(v, list) and i < len(v) else None, out)
    elif isinstance(u, (int, float)) and not isinstance(u, bool):
        if not isinstance(v, (int, float)) or abs(u - v) > 1e-6 + 1e-6 * abs(u): out.append((u, v))
    elif u != v:
        out.append((u, v))
d = []; walk(a, b, d)
print(f"{x}: {len(d)} differences from the archived result")
PYEOF
  done
  echo "== comparison of a sample of every analysis of this revision with the archived output"
  $PY experiments/audit_v12/verify_v120.py; echo "exit code $?"
  echo "== quick check of the version 1.0.0 simulation under this version (20 literatures per base scenario)"
  $PY - << 'PYEOF'
import json, sys
sys.path.insert(0, "experiments/validation")
import validation_v11 as V
import benchmark as B
mis, tot = 0, 0
for i, sc in enumerate(B.ALL_SCEN):
    arch = json.load(open(f"experiments/simulation/benchmark_reps_{i}.json"))
    for r in range(20):
        x = V.one((sc, 7000 * (i + 1) + r))
        for rn in ("R1", "R3"):
            a, b = x[f"legacy {rn}"], arch[r][f"SIGNET {rn}"]
            same = a["cover"] == b["cover"] and a.get("benefit", False) == b["benefit"]
            mis += not same; tot += 1
print(f"legacy settings agree with the archived results in {tot - mis} of {tot} analyses")
PYEOF
  if command -v Rscript > /dev/null; then
    echo "== R interface"; Rscript r/tests/test_signet_r.R; echo "exit code $?"
    echo "== R reconstruction (second implementation, no call to the Python package)"
    Rscript r/reconstruction/signet_reconstruction.R; echo "exit code $?"
  else
    echo "== Rscript not found: the R interface and the R reconstruction were not run"
  fi
  echo "== done $(date)"
} 2>&1 | tee "$LOG"
echo "Log written to $LOG"
