# SIGNET: signed evidence triangulation

Reference implementation and reproduction archive for the manuscript
"Signed evidence triangulation: what biased evidence can certify when the direction of bias is known but its size is not"
(Majd A. AbuAlrob; prepared for submission to *Research Synthesis Methods*).

SIGNET synthesizes evidence from designs with different biases when the direction of each bias is better known than its size.
Directions of bias, structural links between designs and design validity are encoded as named groups of linear constraints. The package then:

1. decides **capacity**: which directional claims (benefit, harm) the design of a body of evidence can certify at all, before any result is seen;
2. computes **certificates**: the weighted contrast of cells and the named assumptions behind each bound;
3. gives **finite-sample inference** by inverting a conditional chi-squared test for linear inequalities (valid in finite samples, with size at most the nominal level, for Gaussian estimates with known variances, and approximate when variances are estimated; the reported set is conservative);
4. reports **tolerance sets**, the **refutation distance** and the **default**, the primary result, which tolerates one failed assumption group or one broken evidence cell, with its capacity and, beside it, the set without tolerance, and with the minimal repairs when the assumption set is rejected.

Version 1.1.0 (18 November 2026) responds to an audit of version 1.0.0: it brackets the outer limits of every confidence set and reports the outer set, checks the hull densely by default, returns a validated upper bound for each directional p-value, makes the certificate choice invariant to how constraints are written, reports whether each certificate is unique, accepts correlations between cells and offers certificates frozen at a chosen covariance. See `CHANGELOG.md` and Supplementary S11 of the manuscript.

Version 1.2.0 (22 November 2026) responds to a forensic audit of the submission package. It changes no numerical result of versions 1.0.0 and 1.1.0. It reports no directional p-value for an empty set, computes the capacity of the default (`default_capacity`), prints the default with its capacity and the set without tolerance in the report, names the absolute numerical thresholds (`signet.core.TOL`) and warns about unusual scales. The archive adds the analyses of Supplementary S12 (`experiments/audit_v12/`) and a second implementation in R (`r/reconstruction/`). Version 1.2.0 has been run only in the development environment of the language model that wrote it; `scripts/reproduce_v120.sh` re-runs its checks. No independent person has verified the code.

Status: research software. It has been evaluated in simulation and in a small retrospective real-data series only. It is not validated for clinical or policy decisions.
The method, code, analyses and text were generated with a large language model under the author's direction; see `AI_USE.md`.

## Installation

Python 3.13 was used. Exact package versions are pinned.

```
python -m pip install -r requirements-lock.txt
python -m pip install -e .
```

On macOS, `python3` may be needed in place of `python`.

## Quick start: the worked example of Figure 1

```
python -m signet validate examples/sigmoidoscopy_distal.json
python -m signet run examples/sigmoidoscopy_distal.json --out sigmoidoscopy.json --report sigmoidoscopy.md
```

The report gives, for each rung, the capacity (whether each side of the effect can be bounded), the 95% SIGNET set, the set tolerant to one failed group, the p-values for benefit and harm, the model test and, in its last lines, the primary result (the default), the set at the same rung without tolerance and the capacity of the default. When one cell measures the target, the default cannot support a claim whatever the data, and the report says so.

## Running SIGNET from R

```r
source("r/signet.R")                 # needs the Python package installed, and the R package jsonlite
signet_capacity("examples/sigmoidoscopy_distal.json")
signet_summary(signet_run("examples/sigmoidoscopy_distal.json"))
```

See `r/README.md`; `Rscript r/tests/test_signet_r.R` reproduces the worked example of Figure 1. The R functions call the Python engine: they are an interface, not an independent implementation.

`r/reconstruction/signet_reconstruction.R` is different: it is a second implementation of the method in R, written from the mathematics, that does not call the Python package. `Rscript r/reconstruction/signet_reconstruction.R` (R packages jsonlite, lpSolve and quadprog) recomputes the four example specifications and compares 80 quantities with the results of the package. The same language model wrote both implementations, so the agreement is a consistency check, not independent verification.

## Running your own SIGNET review

Figure 1 of the manuscript shows the 13 steps and Box 1 shows how a review team can use them with ROBINS-I; this table shows where each step goes in the JSON specification and where its result appears. Supplementary S3 documents every field.

| Step | What to do | Where it goes |
|---|---|---|
| 1 | Define the estimand: population, intervention, comparator, outcome, scale | `name`, `field` (one context, or several if populations differ), `targets` |
| 2 | List the evidence cells: one per design, population and outcome, including controls | `cells`: `id`, `design`, `estimate`, `lower`, `upper` (ratio scale) or `scale: "log"` with `se`, `kernel` (0 for a control the exposure cannot affect), `source` |
| 3 | List the bias processes in each cell | `processes` in each cell |
| 4 | Encode named assumption groups, each with a rung and a dated source | `assumptions`: `type` (sign, zero, share, order, fieldsign, monotone, bound, linear), `process`, `direction` or `rows`, `rung`, `source` |
| 5 | Check capacity for each claim | `capacity` in the result: `above: false` means benefit cannot be certified at that rung |
| 6 | Fix the certificate weights and register | give the expected standard error of every cell in `certificate_se`, so that the weights are fixed before extraction; `claims.benefit.certificate` and `claims.harm.certificate`; the report prints the input SHA-256 |
| 7 | Extract and pool the estimates | enter the pooled estimate of each cell; if two cells draw on the same studies, give their correlation in `correlations` (or run the analysis over a range of values) |
| 8 | Test the assumption set | `model_test`; `minimal_repairs` when it is rejected |
| 9 | Compute the sets | `set` per rung, `tolerance_k1` (and `k2` if requested) |
| 10 | Compute claim p-values | `claims.benefit.p`, `claims.harm.p` |
| 11 | Report the default as the primary result, with its capacity and the set without tolerance beside it | `default`, `default_capacity` and `n_cells_measuring_target` at the highest rung; the last lines of the report |
| 12 | Report what carries each bound | certificates, `identified_at_estimates`, model test, tolerance sets |
| 13 | Report structural gaps | each claim with an unbounded side, with the design or link that would close it |

Capacity and certificate weights do not depend on the estimates. Steps 5 and 6 can therefore be completed before extraction by running the specification with placeholder estimates of 1.0 and the expected confidence limits; in the worked example this reproduced the capacity map and the certificate of the final analysis exactly. Register that specification and its hash before entering the real estimates.

## Layout

| Path | Contents |
|---|---|
| `signet/` | `core.py` (model, capacity, certificates, test, sets, tolerance, default, repairs), `spec.py` (JSON specification layer), `cli.py` |
| `tests/` | `test_core.py` (15 checks), `test_spec.py` (6 checks), `test_cell_tolerance.py` (2 checks of the default), `test_v11.py` (14 checks of version 1.1.0), `test_v12.py` (13 checks of version 1.2.0), `verify_crosssolver.py` (independent CVXPY verification) |
| `audit/` | Independent checks of the final audit and their outputs; the report of the R reconstruction cited in versions 8 to 10, whose source was never archived (kept for the record; no longer cited as evidence) |
| `examples/` | Four real-data specifications with their outputs (sigmoidoscopy, hormone therapy, vitamin E, and the method-agreement check with trial validity as a named assumption) |
| `data/realdata/` | Extracted estimates with access URL and location for each value; sources that could not be used |
| `experiments/realdata/` | Frozen protocol, its SHA-256 hash, amendments, analysis script, post hoc diagnosis, the default sets and outputs |
| `experiments/simulation/` | Simulation scripts, archived outputs, table and contrast generators |
| `figures/` | Figure script and outputs (`AbuAlrobFig1` to `AbuAlrobFig5`) |
| `prompts/` | The prompts that directed the work, verbatim |
| `r/` | R interface (`signet.R`) and its test; `r/reconstruction/`: second implementation in R with its comparison report |
| `experiments/validation/` | Version 1.1.0 validation, paired re-runs and paired contrasts (Supplementary S11) |
| `experiments/audit_v12/` | Analyses added after the forensic audit (Supplementary S12): comparator specifications, ablation of the Bonferroni split, three readings of an assumption group, held-out check of the default, real-data sensitivity analyses, certificate weights, threshold sensitivity, and `verify_v120.py` |
| `REPRODUCTION_LOG.md` | Record of the re-run of version 1.0.0 on a second computer (macOS), 18 November 2026 |
| `REPRODUCTION_RECORD_v1.1.0.md` | Commands, versions and results for version 1.1.0, in the development environment and on the author's computer |
| `reproduction_logs/` | Logs written by `scripts/reproduce_v110.sh` on the author's computer, 19 November 2026, as uploaded by the author |
| `scripts/reproduce_v110.sh` | Re-runs the version 1.1.0 checks and writes a log (about 10 minutes) |
| `scripts/reproduce_v120.sh` | Re-runs the version 1.2.0 checks and writes a log (about 15 minutes) |
| `REPRODUCTION_RECORD_v1.2.0.md` | Commands, versions and results for version 1.2.0, in the development environment only |

## Reproduction

Approximate times are for two cores.

```
python tests/test_core.py                               # 15 checks, about 10 min
python tests/test_spec.py                               # 6 checks, under 1 min
python tests/test_cell_tolerance.py                     # 2 checks of the default
python tests/verify_crosssolver.py                      # needs cvxpy
python experiments/realdata/run_realdata.py             # Table 5 and Figure 5 inputs, under 1 min
python experiments/realdata/default_sets_genesis6.py    # the default in the real-data series
python experiments/realdata/posthoc_contradiction.py    # post hoc diagnosis
python audit/independent_checks.py                      # capacity, test, size, inversion, certificates
python audit/default_equivalence.py                     # package default against the simulation code
cd experiments/simulation
python benchmark.py 400                                 # base scenarios, about 25 min
python benchmark_v2.py comparators 400                  # established comparators, 1 min
python benchmark_v2.py stress 400                       # stress scenarios, about 13 min
python benchmark_v2.py hetallow 400                     # heterogeneity allowance, 4 min
python benchmark_v2.py hetallow 400 0.05                # allowance at half the simulated SD, 4 min
python benchmark_v2.py sweep 200                        # magnitude sweep, about 25 min
python benchmark_v2.py bounded 400                      # magnitude rung, about 30 min
python benchmark_v2.py cellrobust 200                   # cell-level tolerance, about 35 min
python benchmark_v2.py cellhet 200                      # cell-level tolerance with the allowance, about 8 min
python make_tables.py                                   # tables_simulation.md
python contrasts.py                                     # contrasts.md (post hoc contrasts)
cd ../..
python figures/make_figures.py
```

All simulations use fixed seeds, so archived JSON outputs should be reproduced exactly with the pinned versions. `REPRODUCTION_LOG.md` records an independent re-run of the tests, the real-data analysis, the simulation tables and the audit checks on a second computer, with the expected and observed result of each command. For version 1.1.0, run `sh scripts/reproduce_v110.sh` from the repository root in a Python environment with NumPy, SciPy, Clarabel and CVXPY installed; the author's log of 19 November 2026 is in `reproduction_logs/`. For version 1.2.0, run `sh scripts/reproduce_v120.sh`; the full analyses of Supplementary S12 are listed at the top of that script and in Supplementary S9. The author has not yet re-run version 1.2.0.

## Use of generative AI

The method, proofs, code, analyses, figures and text were generated by Claude Opus 5.5 (Anthropic) in response to prompts written by the author, who is accountable for the content. The re-run of version 1.0.0 on a second computer was carried out by the author with assistance from ChatGPT (OpenAI); the re-run of version 1.1.0 was carried out by the author with set-up instructions from Claude Sonnet 5.5 (Anthropic). On 22 November 2026, after a forensic audit supplied by the author, Claude Fable 5.2 (Anthropic) re-ran the affected analyses, wrote the second implementation in R and revised the software to version 1.2.0. `AI_USE.md` gives the full record and `prompts/` the prompts, verbatim.

## Licence and citation

MIT licence (`LICENSE`). Citation metadata are in `CITATION.cff`; GitHub shows them through "Cite this repository". Each release is archived on Zenodo with its own DOI; please cite the DOI of the version used.
