# SIGNET Reproducibility Verification Log

**Verification record date:** 18 November 2026  
**Package:** `SIGNET_RSM_submission_package_v6`  
**Scope:** Commands and checks actually executed or evidenced in this conversation.

## Interpretation of statuses

- **PASS** — the specified check was executed and the observed result met its stated acceptance criterion.
- **FAIL** — the check was executed and the observed result did not meet its acceptance criterion.
- **NOT RUN** — the check did not complete, was blocked, or the available evidence does not establish that it ran to completion.
- Where a result was reported by the user from their Mac Terminal or shown in an uploaded log/screenshot, that provenance is stated. This log does not claim that those commands were executed on the assistant's runtime.

## Check register

| Check | Command | Expected result | Observed result | Status |
|---|---|---|---|---|
| ZIP integrity | ZIP CRC/integrity check (command not retained in conversation record) | Archive readable; no CRC errors | Archive integrity check reported OK | PASS |
| Package SHA-256 manifest | SHA-256 verification against `SHA256SUMS` (command not retained in conversation record) | All manifest entries match; no missing, mismatched, or unlisted files | 165 entries matched; zero mismatch, missing, unlisted, or parse errors | PASS |
| Real-data protocol hash | SHA-256 of `experiments/realdata/PROTOCOL_realdata_v1.md` | Match recorded protocol hash | Observed hash: `d8c5c25be363065e6b70412017bd25edf0663064facac9b21abf1032a6a7e399` | PASS |
| Initial runtime dependency check (assistant runtime) | Python/package-version inspection (exact command not retained) | Required package versions available | Python 3.13.5; NumPy 2.3.5; SciPy 1.17.0; CVXPY and Clarabel unavailable; environment did not match specified versions | NOT RUN |
| Initial specification test (assistant runtime) | `python tests/test_spec.py` | Test completes and passes | Exit code 1; blocked by `ModuleNotFoundError: No module named 'clarabel'` (also an artifact-tool warmup traceback) | NOT RUN |
| Initial real-data script (assistant runtime) | `python experiments/realdata/run_realdata.py` | Script completes and reproduces archived real-data results | Exit code 1; blocked by missing `clarabel`; comparison could not be performed in that run | NOT RUN |
| Simulation table rebuild — simulation tables | `python make_tables.py` | Regenerated `tables_simulation.md` identical to archived copy | Exit code 0; regenerated file identical to archived copy in temporary extraction | PASS |
| Simulation table rebuild — contrasts | `python contrasts.py` | Regenerated `contrasts.md` identical to archived copy | Exit code 0; regenerated file identical to archived copy in temporary extraction | PASS |
| Manuscript Table 3 audit (initial simulation rebuild) | `python audit/check_manuscript_numbers.py ../../01_manuscript/AbuAlrob_SIGNET_manuscript.md` | 90/90 Table 3 cells match regenerated outputs | Printed `Table 3 cells matching regenerated outputs: 90/90`; exit code 0 | PASS |
| Mac specification test | `python tests/test_spec.py` | Test completes; all assertions pass | User-reported Mac output included multiple `PASS` lines and `ALL PASS`; installed versions included Python 3.13.7, NumPy 2.5.3, SciPy 1.18.1, CVXPY 1.9.3, Clarabel 0.11.1. Exit code was blank in the uploaded log, so exit status itself was not confirmed there. | PASS |
| Mac cell-tolerance test | `python tests/test_cell_tolerance.py` | Test completes; all assertions pass | User-reported Mac output included two `PASS` lines and `ALL PASS`; exit code was blank in the uploaded log, so exit status itself was not confirmed there. | PASS |
| Mac cross-solver test | `python tests/verify_crosssolver.py` | Cross-solver checks complete and agree | User-reported output: 60/60 ledgers agree; max solver difference `8.58e-09`; certificate SD and estimate agree to shown precision; identified-set max difference `3.44e-09`; exact/fast rank agree at 15/15 points. Uploaded log did not show an `ALL PASS` line or exit code. | PASS |
| Mac core test | `python tests/test_core.py` | Exactly 15 lines beginning `PASS`; final line `ALL PASS` | User-reported Terminal output: 15 PASS lines; final line `ALL PASS`; exit code 0 | PASS |
| Mac real-data regeneration | `python experiments/realdata/run_realdata.py` | Exit code 0 and regenerated results identical to archived results | Uploaded log: exit code 0; line-by-line comparison reported `NO DIFFERENCES FOUND`; diff exit code 0 | PASS |
| Mac GENESIS VI default-set generation | `python experiments/realdata/default_sets_genesis6.py` | Required Default-column values appear in printed table | Uploaded log: exit code 0. Hormone-therapy cardiovascular rows, sigmoidoscopy, and vitamin E showed `0 to inf`; beta-carotene all-designs showed `0.97 to 1.07`; beta-carotene non-randomized showed `0.75 to inf`. | PASS |
| Simulation table rebuild — `tables_simulation.md` | `python make_tables.py` (Mac, from `experiments/simulation`) | New file identical to `_ARCHIVED.md` backup | User screenshot showed `IDENTICAL: tables_simulation.md`; script exit code 0 | PASS |
| Simulation table rebuild — `contrasts.md` | `python contrasts.py` (Mac, from `experiments/simulation`) | New file identical to `_ARCHIVED.md` backup | User screenshot showed `IDENTICAL: contrasts.md`; script exit code 0 | PASS |
| Simulation table rebuild — `paired_cell_vs_group.md` | `python paired_cell_vs_group.py` (Mac, from `experiments/simulation`) | New file identical to `_ARCHIVED.md` backup | User screenshot showed `IDENTICAL: paired_cell_vs_group.md`; script exit code 0 | PASS |
| Manuscript Table 3 audit (three-table rebuild) | `python audit/check_manuscript_numbers.py ../../01_manuscript/AbuAlrob_SIGNET_manuscript.md` | Print `Table 3 cells matching regenerated outputs: 90/90` | User screenshot showed `Table 3 cells matching regenerated outputs: 90/90` | PASS |
| Independent check A — primal LP vs Farkas dual | `python audit/independent_checks.py` | 4500/4500 agree | Uploaded log: `4500/4500 agree` | PASS |
| Independent check B — solver comparison | `python audit/independent_checks.py` | Maximum absolute difference approximately `2e-07` or smaller | Uploaded log: maximum absolute difference `2.01e-07`; maximum relative difference `2.28e-08` | PASS |
| Independent check C — four rejection rates | `python audit/independent_checks.py` | Each of four rejection rates between 0.035 and 0.060 | Uploaded log at alpha 0.05: 0.0375, 0.0495, 0.0460, 0.0405 | PASS |
| Independent check D — coverage and CC-set containment | `python audit/independent_checks.py` | SIGNET set coverage ≥0.95; `theta accepted but outside reported CC set: 0` | Uploaded log: SIGNET coverage 0.967; CC coverage 0.983; outside reported CC set 0. Additional dense scan found accepted values outside the set in 1/25 literatures. The explicitly specified criteria were met; the dense-scan diagnostic is noted for follow-up. | PASS |
| Independent check E — certificate residual | `python audit/independent_checks.py` | Largest residual below `1e-08` | Uploaded log: largest L1 residual `5.23e-09` | PASS |
| Independent checks script completion | `python audit/independent_checks.py` | Exit code 0 | Uploaded log: exit code 0 | PASS |
| Default equivalence — benefit claims | `python audit/default_equivalence.py` | 24/24 claims agree | Uploaded log: same benefit claim 24/24 | PASS |
| Default equivalence — coverage | `python audit/default_equivalence.py` | 24/24 coverage results agree | Uploaded log: same coverage 24/24 | PASS |
| Default equivalence script completion | `python audit/default_equivalence.py` | Exit code 0 | Uploaded log: exit code 0 | PASS |
| Trial benchmark 1 | `python benchmark.py 10` (run inside `experiments/simulation_TRIAL`) | Completes without errors; report duration | Uploaded log: exit code 0; duration 30 seconds | PASS |
| Trial benchmark 2 | `python benchmark_v2.py comparators 10` (run inside `experiments/simulation_TRIAL`) | Completes without errors; report duration | Uploaded log: exit code 0; duration 2 seconds | PASS |

## Notes and qualifications

1. The initial assistant-runtime attempts at `test_spec.py` and `run_realdata.py` were blocked by missing Clarabel. Later Mac runs supplied evidence of successful execution for the corresponding checks. The initial blocked attempts remain recorded separately as **NOT RUN**.
2. The uploaded three-test Mac log showed `ALL PASS` for `test_spec.py` and `test_cell_tolerance.py`, but did not show exit codes for those scripts. The cross-solver output showed the numerical checks but no final `ALL PASS` marker or exit code. Their PASS statuses above reflect the observed expected test outputs and the user's subsequent report that all three passed; the missing per-script exit-code evidence is disclosed.
3. The independent-check dense scan reported accepted values outside the reported set in 1/25 literatures, despite the separate requested count `theta accepted but outside reported CC set: 0`. The requested acceptance criterion was met, but the additional diagnostic should be investigated and transparently addressed before submission.
4. The small-repetition trial benchmarks are execution smoke tests only. Their numerical summaries are not expected to reproduce manuscript simulation estimates.
5. The GENESIS VI shell wrapper initially printed `RESULT: REVIEW TABLE EXPECTATIONS`; this was not a test failure. Manual inspection of the uploaded table confirmed all specified default values.
6. The simulation-table rebuild script initially had a zsh filename-array issue and stopped before running the checks. A corrected zsh-compatible procedure was subsequently used; all three files were reported identical and the 90/90 audit was observed.
7. No status in this log is based solely on an intended command. Checks that did not complete are marked **NOT RUN**.

## Overall record

All explicitly requested acceptance checks in the later Mac runs were reported as satisfied, including the four test scripts, real-data regeneration, GENESIS VI defaults, simulation-table comparisons, the 90/90 manuscript audit, independent numerical checks, default equivalence, and both trial benchmarks. The qualifications above remain part of the reproducibility record and should not be omitted from submission documentation.
