# Reproduction record, signet version 1.1.0 (18 November 2026)

**Who ran this.** Every command below was executed by Claude Opus 5.5 (Anthropic) in its development environment, a Linux container with two cores, on 18 November 2026. REPRODUCTION_LOG.md is the author's re-run of version 1.0.0 on a second computer. On 19 November 2026 the author ran `sh scripts/reproduce_v110.sh` for version 1.1.0; that re-run and its scope are recorded at the end of this file.

## Environment

| Item | Version |
| --- | --- |
| Operating system | Linux 6.18.44-fc-v595 |
| Python | 3.13.18 |
| NumPy | 2.5.3 |
| SciPy (HiGHS) | 1.18.1 |
| Clarabel | 0.11.1 |
| CVXPY (cross-solver check only) | 1.9.3 |
| R (interface tests only) | R version 4.3.3 (2024-02-29) -- "Angel Food Cake" |
| signet | 1.1.0 |

## Tests

| Command | Result |
| --- | --- |
| `python tests/test_core.py` | ALL PASS (15 checks) |
| `python tests/test_spec.py` | ALL PASS (6 checks) |
| `python tests/test_cell_tolerance.py` | ALL PASS (2 checks) |
| `python -m pytest -q tests/test_v11.py` | 14 passed |
| `Rscript r/tests/test_signet_r.R` | ALL PASS (7 checks) |

## Analyses

| Command | Output | Result |
| --- | --- | --- |
| `python experiments/realdata/v11_checks.py` | experiments/realdata/v11_checks.md, .json | 6 numbers differ from the archived v1.0.0 results, all of them p-values (4 grid suprema corrected upward, one certificate p-value in the fifth decimal, 1 directional p-value); no set limit differs |
| `python -m signet run examples/<name>.json --out ... --report ...` (both examples) | examples/*.result.json, *.report.md | regenerated; reports identical to version 1.0.0; results differ only in two corrected grid suprema (p_cc) and in added fields |
| `python experiments/validation/validation_v11.py 400` | experiments/validation/validation_reps_0.json to _9.json | archived settings reproduced 8,000 of 8,000 rung-1 and rung-3 analyses (Supplementary S11) |
| `python experiments/validation/paired_reruns.py comparators hetallow bounded sweep cellrobust` | experiments/validation/reps_*.json | archived summaries reproduced in the 14 scenarios with a deterministic generator; largest absolute difference 0.0000 |
| `python experiments/validation/validation_default_v11.py 20` | experiments/validation/validation_default_reps.json | archived settings reproduced the archived default in 300 of 300 literatures; hardened settings changed coverage in 0 |
| `python experiments/simulation/benchmark_v2.py stress 400 3` and `python experiments/simulation/rerun_missing_controls_v110.py` | benchmark_v2_stress.json, benchmark_v2_bounded.json, benchmark_v2_cellrobust.json, all_scenarios_summary.json | the missing-controls scenario regenerated after the generator fix; largest change in a reported proportion 0.013 |
| `python experiments/simulation/make_tables.py`, `python experiments/simulation/paired_cell_vs_group.py`, `python ms/make_s5.py`, `python figures/make_figures.py` | tables, S5, figures | regenerated; Figure 2 changed slightly, Figures 1, 3, 4 and 5 are pixel-identical to version 8 |
| `python experiments/validation/paired_contrasts.py` | experiments/validation/paired_contrasts.md, validation_summary.md | Table 4 of the manuscript and Supplementary S11 |

Run times (two cores): validation_v11.py about 25 minutes for 4,000 literatures; paired_reruns.py about 40 minutes; validation_default_v11.py about 20 minutes.

The author script `sh scripts/reproduce_v110.sh` was itself run once in the development environment: every step passed (tests, real-data comparison, 400 of 400 quick-validation analyses agreeing with the archive, R interface) in about 11 minutes.

The 1.0.0 benchmark scripts were not re-run in full; they pass `outer_scan=0`, so that they reproduce version 1.0.0 under version 1.1.0, and the validation above re-ran their SIGNET computations on all 4,000 base literatures.

## Re-run by the author (19 November 2026)

The author ran `sh scripts/reproduce_v110.sh` on a second computer (the author's Mac: macOS, Apple M1) from the unzipped code archive (folder `signet-evidence`; the log reports signet 1.1.0). The code in `signet-evidence-v1.1.0.zip` as delivered on 18 November 2026 is byte-identical to this release. Set-up instructions (a virtual environment and `pip install numpy scipy clarabel cvxpy`) were given by Claude Sonnet 5.5 (Anthropic) in the claude.ai conversation; the commands were run by the author in Terminal. The log is archived as the author uploaded it.

| Item | Value |
| --- | --- |
| Log | `reproduction_logs/REPRODUCTION_v1.1.0_20261119_0940.log` (SHA-256 ed653e5225f50cbd2f0b283c0d0ca5aaf202ab9497bc9629a7e8c7af0f573302) |
| Run | 09:40:51 to 09:47:07, +03:00 |
| Python | 3.13.7 |
| NumPy | 2.5.5 (the lock file pins 2.5.3) |
| SciPy | 1.18.1 |
| Clarabel | 0.11.1 |
| signet | 1.1.0 |

| Step | Result on the author's computer |
| --- | --- |
| `tests/test_core.py` | ALL PASS (15 checks), exit code 0 |
| `tests/test_spec.py` | ALL PASS (6 checks), exit code 0 |
| `tests/test_cell_tolerance.py` | ALL PASS (2 checks), exit code 0 |
| `tests/test_v11.py` | 14 of 14 passed, exit code 0 |
| `experiments/realdata/v11_checks.py` | exit code 0; the first 12 lines of the regenerated report, including the four p-values printed to full precision, are identical to the archived report. The script prints only these lines, so the rest of the regenerated report was not compared. |
| Quick validation (20 literatures per base scenario, legacy settings) | 400 of 400 rung-1 and rung-3 analyses agree with the archived version 1.0.0 results |
| `r/tests/test_signet_r.R` (R interface) | ALL PASS (7 checks), exit code 0 |

Scope. This re-run covers the unit checks, the regeneration of the real-data checks, a 400-analysis sample of the base simulation and the R interface, which calls the Python engine. It does not cover the full validation (8,000 analyses), the paired re-runs, the regeneration of the missing-controls scenario, the tables or the figures; those were run only in the development environment. The R interface test is not a test of the independent R reconstruction, whose source code is not archived.

An earlier attempt at 09:38 stopped at the first import because NumPy was not installed in the system Python; no check ran. Its log is kept as `reproduction_logs/REPRODUCTION_v1.1.0_20261119_0938.log` (SHA-256 9ff508aa7545f2a39a3510091eeb28b91b8f6c3a9660baac4411ca7b13f405c8). Its real-data section shows the archived report, printed by `head`, not a regenerated one.
