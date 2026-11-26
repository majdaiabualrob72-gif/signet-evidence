# Changelog

## 1.2.0 (22 November 2026): response to the forensic audit of the submission package

No numerical result of versions 1.0.0 and 1.1.0 changes: under version 1.2.0 the real-data checks of version 1.1.0 return the archived output, and the new simulation runs reproduce every archived record that they recompute (Supplementary S12.13).

Changes to the inference code (`signet/core.py`):

- **Empty sets.** `claim_pvalue` takes `alpha` and `empty`; for an empty set it returns `p = None`, `p_upper = None`, `set_empty = True` and a note. When the test rejects every value of the effect, the formula min(1, 4 p_cert, 2 p_T) is small for benefit and for harm alike; the data then contradict the assumption set and support no direction.
- **Capacity of the default.** New function `default_capacity(model, include, cf, droppable=None)`: the default is a union of sets, so it is bounded on a side only if capacity is bounded on that side for the full assumption set, for the set without each droppable group and for each analysis that omits one cell. The function returns `above`, `below` and the omissions that leave each side unbounded. It uses no estimate.
- **Named thresholds.** The absolute thresholds of the test are collected in `TOL = dict(T_zero=1e-9, slack=1e-4, implicit=1e-7)`; their values are unchanged.

Changes to the specification layer (`signet/spec.py`):

- The result records `default_capacity` and `n_cells_measuring_target` at the highest rung, and `warnings`.
- The report prints the default, the set at the same rung without tolerance (conditional on every named assumption) and the capacity of the default; when one cell measures the target, it states that the default cannot support a directional claim and that its silence carries no information about the effect. It prints "not reported (empty set)" in place of a directional p-value for an empty set.
- `scale_warnings`: a warning when a cell's standard error is below 0.01 or above 100, because the bracketing step, the bisection tolerance and the thresholds above are absolute.

Version numbers: `signet.__version__` and `pyproject.toml` now agree (the latter still read 1.0.0 in version 1.1.0).

New tests: `tests/test_v12.py` (13 checks). The 37 earlier checks pass unchanged.

New examples: `examples/vitamin_e_chd.json` and `examples/beta_carotene_chd_trial_validity_named.json`, with their results and reports; the two earlier examples were re-run and are numerically identical, with the new fields added.

New analyses (`experiments/audit_v12/`, Supplementary S12): `sim_audit.py` (comparator specifications, ablation of the Bonferroni split, three readings of an assumption group, held-out check of the default), `analyse_audit.py`, `realdata_audit.py` and `realdata_audit_models.py`, `certificate_weights.py`, `tolerance_sensitivity.py` and `verify_v120.py`.

Second implementation in R: `r/reconstruction/signet_reconstruction.R` and its report. It does not call the Python package. The same language model wrote it, so it is not independent verification.

Figures: `figures/make_figures.py` regenerates Figure 1 (step 11 reports the capacity of the default and the set without tolerance) and Figure 5 (published estimates labelled as such, three intervals for the vitamin E pool, BA-S with shared terms, and the default with its lack of capacity). Figures 2 to 4 are unchanged.

`scripts/reproduce_v120.sh` re-runs the checks of this version and writes a log. Version 1.2.0 has been run only in the development environment.

Records: `REPRODUCTION_RECORD_v1.2.0.md` (what was run, with results), `reproduction_logs/REPRODUCTION_v1.2.0_20261122_1550_development_environment.log` and `REFERENCE_VERIFICATION.md` (78 references with the record each was checked against).

## 1.1.0, documentation update (19 November 2026)

No change to code or outputs. `AI_USE.md` also records the shortening of the final manuscript. Added `reproduction_logs/` with the author's run of `scripts/reproduce_v110.sh` and updated `REPRODUCTION_RECORD_v1.1.0.md`, `AI_USE.md` and `README.md` to report it.

## 1.1.0 (18 November 2026): response to the audit of version 8

Changes to the inference code (`signet/core.py`):

- **Outer set and bracketing.** `conf_set` now brackets each outer limit (where T(θ0) exceeds the n-degree quantile) up to `max_span = 1e4` from the estimate. Version 1.0.0 stopped the search at `span = 6`, which could truncate a very wide set without warning. A limit that cannot be bracketed is reported as infinite, with a warning in the result. The outer set `{θ0 : T(θ0) <= chi2_n(1 - alpha)}` is returned with every set (`outer`) and can be reported instead of the hull (`inversion="outer"`): because T is convex in θ0, it is an interval with exact end points, it contains every accepted value, and it is itself a valid confidence set.
- **Dense check by default.** `outer_scan` now defaults to 16: after the hull is refined, 16 points on each side between the hull and the outer limit are tested and the hull is extended to any accepted value (`hull_extended` in the result). `outer_scan=0` reproduces version 1.0.0; the archived benchmark scripts now pass it explicitly.
- **Validated directional p-value.** `claim_pvalue` adds the null value nearest the estimate to its grid and returns `p_cc_upper` and `p_upper`. Because the conditional p-value at θ0 never exceeds the n-degree tail probability at T(θ0), and T is convex, `p_cc_upper` bounds the supremum over the whole null region, so `p_upper = min(1, 4 p_cert, 2 p_cc_upper)` is a valid directional p-value.
- **Scale-invariant certificates.** Constraint rows are scaled to unit length when compiled (`Model.normalise_rows = True`), so the first stage of the certificate choice, which minimises the total multiplier on inequality assumptions, no longer depends on how an assumption is written. Set `normalise_rows = False` to reproduce version 1.0.0.
- **Frozen certificates and uniqueness.** `cert_bounds`, `signet_set` and `claim_pvalue` accept `cert_V`, a covariance used only to choose the certificate weights (for example within-study variances, which do not depend on the estimates); the standard deviation of the certified combination is still computed from V. `certificate_range` reports whether the certificate is unique, in which case its choice cannot depend on V or on the estimates. `check_certificate` verifies dual feasibility of a given weight vector and returns the offset of the non-homogeneous bound c'x <= λ'μ + ν'h + κ'g.
- **Docstrings** now state that certificate weights depend on the data through V when V is estimated from the same studies and the certificate is not unique.

Changes to the specification layer (`signet/spec.py`):

- New optional field `correlations`: a list of `{cells: [id, id], rho}`; the covariance of the cell estimates becomes S R S. Validation rejects unknown cells, |rho| >= 1 and correlation matrices that are not positive semidefinite.
- New optional field `certificate_se`: expected standard errors for every cell, recorded before extraction; certificate weights are then chosen with them (and with any correlations), so that the choice cannot depend on the estimates. Results record `certificate_selection`.
- New analysis option `inversion` (`"hull"` or `"outer"`).
- Results record `signet_version`, `inversion_checks` (outer set, bracketing, dense-check extension, warnings), `certificate_unique` and, for claims, `p_upper`.

New tests: `tests/test_v11.py` (14 checks). The 23 version 1.0.0 checks pass unchanged.

Simulation code: the generator of the stress scenario with missing controls and variances (`benchmark_v2.literature_v2`) looped over a Python set of design names when imputing variances, so its random draws depended on the per-process string-hash order and the scenario was not reproducible across runs. It now loops over the sorted names; the scenario was regenerated (`experiments/simulation/rerun_missing_controls_v110.py`, `benchmark_v2.py stress 400 3`), and the version 1.0.0 outputs are kept in `experiments/simulation/v1.0.0_missing_controls/`.

New analyses: `experiments/realdata/v11_checks.py`, `experiments/validation/validation_v11.py`, `experiments/validation/paired_reruns.py` and `experiments/validation/paired_contrasts.py` (see Supplementary S11).

## 1.0.0 (18 November 2026)

First archived release.
