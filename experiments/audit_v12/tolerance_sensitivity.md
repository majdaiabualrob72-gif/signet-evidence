# Sensitivity of the results to the numerical thresholds (N10; signet 1.2.0)

Defaults: T_zero = 1e-09, slack = 0.0001, implicit = 1e-07 (signet.core.TOL). Each threshold was multiplied by 0.01 and by 100, one at a time.

| Setting | Real data: largest change in a set limit (log scale) | Real data: largest change in a p-value | Real data: sets changing between empty, bounded and unbounded | Simulation, archived rank: largest change in a limit | coverage or claim outcomes changed (of 300) | Simulation, exact rank: largest change in a limit | outcomes changed (of 300) |
|---|---|---|---|---|---|---|---|
| default | 0.0e+00 | 0.0e+00 | 0 | 0.0e+00 | 0 | 0.0e+00 | 0 |
| T_zero x 0.01 | 0.0e+00 | 0.0e+00 | 0 | 0.0e+00 | 0 | 0.0e+00 | 0 |
| T_zero x 100 | 0.0e+00 | 0.0e+00 | 0 | 0.0e+00 | 0 | 0.0e+00 | 0 |
| slack x 0.01 | 0.0e+00 | 0.0e+00 | 0 | 0.0e+00 | 0 | 0.0e+00 | 0 |
| slack x 100 | 0.0e+00 | 0.0e+00 | 0 | 0.0e+00 | 0 | 0.0e+00 | 0 |
| implicit x 0.01 | 0.0e+00 | 0.0e+00 | 0 | 2.2e-02 | 0 | 2.2e-02 | 0 |
| implicit x 100 | 0.0e+00 | 0.0e+00 | 0 | 1.4e-03 | 0 | 1.4e-03 | 0 |

Real-data analyses: 14 sets (every rung of the four questions, frozen and post hoc encodings of the method-agreement case). Simulation: 300 literatures (20 per scenario). At the default thresholds, exact and archived ranks differed by at most 0.0e+00 in a limit and in 0 outcomes.

## Rescaling the worked example

Estimates and standard errors of the sigmoidoscopy example (rung 2) were multiplied by c. A method without absolute settings would return limits
multiplied by c and unchanged p-values. The bracketing step (0.05), the bisection tolerance (1e-4) and the slack thresholds are absolute.

| c | Typical standard error | Lower limit / c | Upper limit / c | p for benefit | Largest relative change in a limit |
|---|---|---|---|---|---|
| 0.0001 | 2.6e-05 | -1.59167 | -0.07969 | 0.02791 | 6.0e-02 |
| 0.001 | 0.00026 | -1.59226 | -0.07704 | 0.02791 | 2.5e-02 |
| 0.01 | 0.0026 | -1.59344 | -0.07779 | 0.02791 | 3.5e-02 |
| 0.1 | 0.026 | -1.59356 | -0.07560 | 0.02791 | 5.8e-03 |
| 1 | 0.26 | -1.59356 | -0.07517 | 0.02791 | 0.0e+00 |
| 10 | 2.6 | -1.59356 | -0.07510 | 0.02791 | 8.9e-04 |
| 100 | 26 | -1.59356 | -0.07510 | 0.02791 | 9.3e-04 |
| 1000 | 2.6e+02 | -1.59356 | -0.07510 | 0.02791 | 9.4e-04 |
