# SIGNET report: Beta-carotene and coronary heart disease, published triangulation inputs, trial validity named (post hoc encoding, 22 November 2026)

Input SHA-256: 8d6bfd62a4bd13037f8f1d0b595b4e6832c14efb9738fa1e1876fe9e4355c367

Warning: cell D_rct: standard error 0.00769 on the log scale is below 0.01; set limits are computed to an absolute tolerance of 1e-4 (1.3% of this standard error)

## Target: chd

| Rung | Capacity (upper, lower) | 95% SIGNET set (ratio scale) | k = 1 set | p benefit | p harm | Model test p |
| --- | --- | --- | --- | --- | --- | --- |
| rung1 | unbounded, unbounded | 0 to inf | 0 to inf | 1.0000 | 1.0000 | 1.0 |
| rung2 | bounded, bounded | 0.970 to 1.073 | 0.749 to inf | 1.0000 | 0.7539 | 1.0 |
| rung3 | bounded, bounded | 0.986 to 1.019 | 0.970 to 1.073 | 1.0000 | 1.0000 | 0.4034 |

Default at rung3, tolerant to one failed assumption group or one broken evidence cell: 0.970 to 1.073
Beside it, rung3 without tolerance (conditional on every named assumption): 0.986 to 1.019
Capacity of the default: upper limit can be bounded; lower limit can be bounded.

