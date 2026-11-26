# Reproduction record for version 1.2.0

Version 1.2.0 (22 November 2026) responds to a forensic audit of the submission package. This record lists what was run, where, and with what result.

**The runs in the next sections were executed by a language model (Claude Fable 5.2, Anthropic) in its development environment, a Linux container with two cores. The author re-ran the reproduction script on a second computer the same day (last section). No independent person has verified the work.**

## Environment

python 3.13.18 numpy 2.5.3 scipy 1.18.1 clarabel 0.11.1 signet 1.2.0; CVXPY 1.9.3; R 4.3.3 with jsonlite 1.8.8, lpSolve 5.6.20 and quadprog 1.5.8.

## Reproduction script

`sh scripts/reproduce_v120.sh`, run from the repository root on 22 November 2026.

- Log: `reproduction_logs/REPRODUCTION_v1.2.0_20261122_1550_development_environment.log` (SHA-256 ce8231ce86019452b16b1914dec6c14ec06eda88195e3f3bb6e90820edb64bb1).
- Unit checks: the five test scripts ended with ALL PASS (15, 6, 2, 14 and 13 checks; 50 in all).
- Exit codes recorded in the log: 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.
- Example specifications: 0 differences from the archived results for each of the four.
- `experiments/audit_v12/verify_v120.py`: 6 comparisons, all PASS (real-data checks of version 1.1.0 under this version; real-data analyses of this revision; comparators; samples of the ablation, group-reading and held-out analyses).
- Version 1.0.0 simulation, 20 literatures per base scenario: 400 of 400 analyses agree with the archive.
- R interface: 7 checks pass. Second implementation in R: 80 comparisons, 0 outside tolerance.

## Full analyses of Supplementary S12 (development environment)

| Command | Result |
|---|---|
| `python3 experiments/audit_v12/sim_audit.py comparators` | `out/comparators__15_scenarios.json`, `out/comparators__sweep.json` |
| `python3 experiments/audit_v12/sim_audit.py worker` (two processes) | 50 jobs: 15 ablation (400 literatures each), 15 group readings (200), 2 group readings (literatures 201 to 400 of two scenarios), 18 held-out (200). The two processes were stopped by the two-hour limit of the environment after 48 jobs; a third process ran the last two jobs from their beginning. Each job writes its file only when it is complete (logs: `experiments/audit_v12/out/worker_A.log` to `worker_D.log`) |
| `python3 experiments/audit_v12/analyse_audit.py` | `audit_numbers.json`, `audit_tables.md` |
| `python3 experiments/audit_v12/realdata_audit.py` | `realdata_audit.json`, `realdata_audit.md` |
| `python3 experiments/audit_v12/certificate_weights.py 400` | `certificate_weights.json`, `certificate_weights.md` |
| `python3 experiments/audit_v12/tolerance_sensitivity.py` | `tolerance_sensitivity.json`, `tolerance_sensitivity.md` |
| `Rscript r/reconstruction/signet_reconstruction.R` | `r/reconstruction/reconstruction_report.md` |
| `python3 figures/make_figures.py fig1 fig4` | Figures 1 and 5 |

## Agreement with the archive

| Quantity recomputed by the new runs | Literatures | Different |
|---|---|---|
| Rung 3: coverage, benefit claim, emptiness | 6,000 | 0 |
| Conditional test alone at rung 3 (base scenarios) | 4,000 | 0 |
| Group tolerance: coverage and benefit claim | 3,000 | 0 |
| Default: package function against the archived simulation code | 3,000 | 0 |
| Random-effects comparator, modified Hartung-Knapp interval | 6,000 | 0 |

The real-data checks of version 1.1.0 (`experiments/realdata/v11_checks.py`), run under version 1.2.0, returned the archived JSON output with no numerical difference; the archived files of version 1.1.0 were left in place.

## Not covered

- The complete simulations of versions 1.0.0 and 1.1.0 under version 1.2.0 (the records above and the sample of the reproduction script were recomputed).
- The full analyses of Supplementary S12 on a second computer (the reproduction script re-runs samples of each).
- Verification by an independent person.

## Author's re-run (22 November 2026)

Computer: MacBook Air (Apple M1), macOS; Python 3.13.7, NumPy 2.5.3, SciPy 1.18.1, Clarabel 0.11.1, signet 1.2.0; R 4.5.1 with lpSolve 5.6.23, quadprog 1.5.8 and jsonlite 2.0.0.

- `reproduction_logs/REPRODUCTION_v1.2.0_20261122_1812.log` (SHA-256 dcb501acf0c8f27642f93ddaa41d2cd993f0573392954d174b48636ac12123de): first attempt; stopped because NumPy was not installed. No check ran.
- `reproduction_logs/REPRODUCTION_v1.2.0_20261122_1813.log` (SHA-256 c1745167f8bc06184ae4e853548ab0db6b222269945f79d03582d115f694bdae): the five test scripts ended with ALL PASS (50 checks); the four example specifications showed 0 differences; the samples of the comparators, ablation, group-reading and held-out analyses and the real-data analyses of this revision were identical to the archive; the version 1.0.0 sample agreed in 400 of 400 analyses; the R interface passed 7 checks; the R implementation agreed in 80 of 80 comparisons.
- One comparison of `verify_v120.py` failed: the real-data checks of version 1.1.0 showed two differences. Both concern the lower limit of the group-tolerant rung-3 set for the control outcome of the hormone-therapy question, at correlations 0 and 0.75. The values lay within 10⁻⁶ of zero on the log scale and differed from the archive by 4.6 × 10⁻⁸ and 6.1 × 10⁻⁸, beyond the strict tolerance of the comparison (10⁻⁹ plus 10⁻⁶ relative) and far below the bisection tolerance of the method (10⁻⁴). No reported number changes. The tolerance was not changed after the run.
