#!/bin/sh
# Re-run the version 1.1.0 checks on your own computer and write a log.
# Usage, from the repository root:   sh scripts/reproduce_v110.sh
# Time: about 10 to 15 minutes on two cores (the quick validation uses 20 literatures per scenario).
# The full validation (400 per scenario, about 60 minutes) is:
#   python3 experiments/validation/validation_v11.py 400
# Requires: Python 3.13 with requirements-lock.txt installed and `pip install -e .`; optionally R with jsonlite.
PY=${PY:-python3}
LOG=REPRODUCTION_v1.1.0_$(date +%Y%m%d_%H%M).log
{
  echo "== SIGNET v1.1.0 reproduction, $(date)"
  echo "== machine: $(uname -a)"
  $PY -c "import sys, numpy, scipy, clarabel; import signet; print('python', sys.version.split()[0], 'numpy', numpy.__version__, 'scipy', scipy.__version__, 'clarabel', clarabel.__version__, 'signet', signet.__version__)"
  for t in tests/test_core.py tests/test_spec.py tests/test_cell_tolerance.py; do
    echo "== $t"; $PY $t; echo "exit code $?"
  done
  echo "== tests/test_v11.py"; $PY tests/test_v11.py; echo "exit code $?"
  echo "== real data under v1.1.0 against the archived v1.0.0 output"
  $PY experiments/realdata/v11_checks.py > /dev/null; echo "exit code $?"
  head -12 experiments/realdata/v11_checks.md
  echo "== quick validation (20 literatures per base scenario, compared with the archive)"
  $PY - << 'PYEOF'
import json, os, sys
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
  if command -v Rscript > /dev/null; then echo "== R interface"; Rscript r/tests/test_signet_r.R; echo "exit code $?"; fi
  echo "== done $(date)"
} 2>&1 | tee "$LOG"
echo "Log written to $LOG"
