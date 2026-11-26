# Certificate weights (finding N11)

Rule of the paper: least total weight on inequality assumptions first, least variance second.
Comparison: the certificate of least variance alone. Literatures with at least two cells measuring the effect.
Simulated literatures per scenario: 400 (archived seeds).

## Simulation, rung 3, upper limit (the side that supports a claim of benefit)

The last three columns compare the one-sided limits (level 0.0125) of the two rules on the log scale: mean of
(limit under the rule minus limit under least variance); a negative value means that the rule gave the tighter limit on average.
Limits closer than 1e-4 on the log scale are counted as the same.

| Scenario | Literatures with a certificate | Most precise cell given zero weight | Median SD ratio (rule / least variance) | 90th centile | Ratio above 1.05 | Mean difference in limit | Rule tighter | Same limit | Rule looser |
|---|---|---|---|---|---|---|---|---|---|
| reference | 389 | 185 (0.476) | 1.022 | 1.583 | 0.463 | +0.004 | 0.308 | 0.463 | 0.229 |
| null effect | 392 | 175 (0.446) | 1.000 | 1.474 | 0.434 | +0.000 | 0.332 | 0.487 | 0.181 |
| magnitudes doubled | 393 | 184 (0.468) | 1.000 | 1.385 | 0.425 | -0.015 | 0.557 | 0.326 | 0.117 |
| dense small biases | 392 | 168 (0.429) | 1.000 | 1.447 | 0.446 | +0.011 | 0.258 | 0.452 | 0.291 |
| 15% of signs wrong | 392 | 161 (0.411) | 1.000 | 1.355 | 0.452 | +0.004 | 0.291 | 0.467 | 0.242 |
| omitted bias | 393 | 186 (0.473) | 1.000 | 1.306 | 0.382 | -0.005 | 0.379 | 0.511 | 0.109 |
| no design opposing benefit | 242 | 125 (0.517) | 1.000 | 1.000 | 0.000 | -0.000 | 0.083 | 0.917 | 0.000 |
| no bias | 392 | 174 (0.444) | 1.000 | 1.432 | 0.444 | +0.016 | 0.115 | 0.513 | 0.372 |
| one-signed literature, null effect | 0 | not applicable | | | | | | | |
| one-signed literature, null effect, magnitudes doubled | 0 | not applicable | | | | | | | |
| extreme outlier | 395 | 195 (0.494) | 1.012 | 1.433 | 0.435 | -0.001 | 0.352 | 0.433 | 0.215 |
| publication bias | 392 | 176 (0.449) | 1.000 | 1.445 | 0.459 | +0.005 | 0.286 | 0.467 | 0.247 |
| heterogeneous populations | 389 | 158 (0.406) | 1.000 | 1.400 | 0.434 | -0.000 | 0.332 | 0.463 | 0.206 |
| missing controls and variances | 377 | 239 (0.634) | 1.028 | 1.586 | 0.469 | +0.002 | 0.371 | 0.435 | 0.194 |
| contradictory evidence | 400 | 172 (0.430) | 1.034 | 1.408 | 0.472 | +0.007 | 0.370 | 0.343 | 0.287 |

## Pooled over the 15 scenarios

| Rung and side | Literatures with a certificate | Most precise cell given zero weight | SD ratio above 1.05 | Mean difference in limit | Rule tighter | Same limit | Rule looser |
|---|---|---|---|---|---|---|---|
| Rung 1 upper | 2412 | 2152 (0.892) | 0.000 | -0.000 | 0.194 | 0.802 | 0.004 |
| Rung 1 lower | 5549 | 3813 (0.687) | 0.695 | -0.023 | 0.531 | 0.246 | 0.224 |
| Rung 2 upper | 4071 | 2889 (0.710) | 0.349 | +0.018 | 0.211 | 0.556 | 0.232 |
| Rung 2 lower | 5583 | 3651 (0.654) | 0.896 | -0.044 | 0.633 | 0.093 | 0.274 |
| Rung 3 upper | 4938 | 2298 (0.465) | 0.421 | +0.002 | 0.317 | 0.470 | 0.213 |
| Rung 3 lower | 5649 | 2637 (0.467) | 0.944 | -0.049 | 0.743 | 0.024 | 0.232 |

## Simulation, all rungs and sides

| Scenario | Rung and side | Literatures with a certificate | Zero weight on the most precise cell | Median SD ratio | Ratio above 1.05 | Mean difference in limit |
|---|---|---|---|---|---|---|
| reference | R1 hi | 205 | 175 (0.854) | 1.000 | 0.000 | -0.000 |
| reference | R1 lo | 386 | 287 (0.744) | 1.208 | 0.617 | -0.013 |
| reference | R2 hi | 322 | 233 (0.724) | 1.000 | 0.398 | +0.012 |
| reference | R2 lo | 390 | 271 (0.695) | 1.718 | 0.877 | -0.043 |
| reference | R3 hi | 389 | 185 (0.476) | 1.022 | 0.463 | +0.004 |
| reference | R3 lo | 398 | 188 (0.472) | 1.354 | 0.950 | -0.055 |
| null effect | R1 hi | 193 | 167 (0.865) | 1.000 | 0.000 | -0.000 |
| null effect | R1 lo | 391 | 290 (0.742) | 1.259 | 0.678 | -0.011 |
| null effect | R2 hi | 333 | 240 (0.721) | 1.000 | 0.372 | +0.016 |
| null effect | R2 lo | 392 | 275 (0.702) | 1.771 | 0.901 | -0.040 |
| null effect | R3 hi | 392 | 175 (0.446) | 1.000 | 0.434 | +0.000 |
| null effect | R3 lo | 395 | 174 (0.441) | 1.339 | 0.949 | -0.054 |
| magnitudes doubled | R1 hi | 209 | 193 (0.923) | 1.000 | 0.000 | -0.000 |
| magnitudes doubled | R1 lo | 390 | 281 (0.721) | 1.209 | 0.641 | -0.064 |
| magnitudes doubled | R2 hi | 326 | 240 (0.736) | 1.000 | 0.393 | -0.010 |
| magnitudes doubled | R2 lo | 396 | 274 (0.692) | 1.664 | 0.884 | -0.179 |
| magnitudes doubled | R3 hi | 393 | 184 (0.468) | 1.000 | 0.425 | -0.015 |
| magnitudes doubled | R3 lo | 400 | 191 (0.477) | 1.360 | 0.932 | -0.165 |
| dense small biases | R1 hi | 201 | 177 (0.881) | 1.000 | 0.000 | -0.000 |
| dense small biases | R1 lo | 387 | 278 (0.718) | 1.210 | 0.656 | -0.016 |
| dense small biases | R2 hi | 328 | 227 (0.692) | 1.000 | 0.393 | +0.023 |
| dense small biases | R2 lo | 389 | 260 (0.668) | 1.743 | 0.925 | -0.073 |
| dense small biases | R3 hi | 392 | 168 (0.429) | 1.000 | 0.446 | +0.011 |
| dense small biases | R3 lo | 398 | 166 (0.417) | 1.349 | 0.945 | -0.077 |
| 15% of signs wrong | R1 hi | 199 | 182 (0.915) | 1.000 | 0.000 | -0.000 |
| 15% of signs wrong | R1 lo | 391 | 271 (0.693) | 1.188 | 0.624 | +0.001 |
| 15% of signs wrong | R2 hi | 331 | 229 (0.692) | 1.000 | 0.366 | +0.021 |
| 15% of signs wrong | R2 lo | 395 | 257 (0.651) | 1.620 | 0.889 | -0.012 |
| 15% of signs wrong | R3 hi | 392 | 161 (0.411) | 1.000 | 0.452 | +0.004 |
| 15% of signs wrong | R3 lo | 396 | 158 (0.399) | 1.331 | 0.932 | -0.035 |
| omitted bias | R1 hi | 171 | 153 (0.895) | 1.000 | 0.000 | -0.000 |
| omitted bias | R1 lo | 388 | 266 (0.686) | 1.218 | 0.665 | -0.014 |
| omitted bias | R2 hi | 325 | 235 (0.723) | 1.000 | 0.351 | +0.011 |
| omitted bias | R2 lo | 391 | 262 (0.670) | 1.763 | 0.923 | -0.046 |
| omitted bias | R3 hi | 393 | 186 (0.473) | 1.000 | 0.382 | -0.005 |
| omitted bias | R3 lo | 398 | 187 (0.470) | 1.364 | 0.940 | -0.029 |
| no design opposing benefit | R1 lo | 346 | 198 (0.572) | 1.259 | 0.708 | -0.012 |
| no design opposing benefit | R2 hi | 242 | 125 (0.517) | 1.000 | 0.000 | -0.000 |
| no design opposing benefit | R2 lo | 346 | 183 (0.529) | 1.746 | 0.962 | -0.052 |
| no design opposing benefit | R3 hi | 242 | 125 (0.517) | 1.000 | 0.000 | -0.000 |
| no design opposing benefit | R3 lo | 346 | 183 (0.529) | 1.746 | 0.962 | -0.052 |
| no bias | R1 hi | 195 | 174 (0.892) | 1.000 | 0.000 | +0.000 |
| no bias | R1 lo | 385 | 267 (0.694) | 1.250 | 0.670 | +0.042 |
| no bias | R2 hi | 314 | 227 (0.723) | 1.000 | 0.420 | +0.032 |
| no bias | R2 lo | 388 | 259 (0.668) | 1.758 | 0.912 | +0.098 |
| no bias | R3 hi | 392 | 174 (0.444) | 1.000 | 0.444 | +0.016 |
| no bias | R3 lo | 397 | 175 (0.441) | 1.374 | 0.937 | +0.051 |
| one-signed literature, null effect | R1 lo | 265 | 127 (0.479) | 1.379 | 0.981 | -0.011 |
| one-signed literature, null effect | R2 lo | 265 | 127 (0.479) | 1.379 | 0.981 | -0.011 |
| one-signed literature, null effect | R3 lo | 265 | 127 (0.479) | 1.379 | 0.981 | -0.011 |
| one-signed literature, null effect, magnitudes doubled | R1 lo | 261 | 132 (0.506) | 1.420 | 0.981 | -0.091 |
| one-signed literature, null effect, magnitudes doubled | R2 lo | 261 | 132 (0.506) | 1.420 | 0.981 | -0.091 |
| one-signed literature, null effect, magnitudes doubled | R3 lo | 261 | 132 (0.506) | 1.420 | 0.981 | -0.091 |
| extreme outlier | R1 hi | 215 | 187 (0.870) | 1.000 | 0.000 | -0.000 |
| extreme outlier | R1 lo | 390 | 272 (0.697) | 1.166 | 0.603 | +0.007 |
| extreme outlier | R2 hi | 343 | 251 (0.732) | 1.000 | 0.394 | +0.020 |
| extreme outlier | R2 lo | 394 | 266 (0.675) | 1.735 | 0.878 | -0.025 |
| extreme outlier | R3 hi | 395 | 195 (0.494) | 1.012 | 0.435 | -0.001 |
| extreme outlier | R3 lo | 400 | 195 (0.487) | 1.374 | 0.917 | -0.048 |
| publication bias | R1 hi | 198 | 174 (0.879) | 1.000 | 0.000 | -0.000 |
| publication bias | R1 lo | 390 | 275 (0.705) | 1.198 | 0.636 | -0.009 |
| publication bias | R2 hi | 316 | 215 (0.680) | 1.000 | 0.392 | +0.021 |
| publication bias | R2 lo | 392 | 263 (0.671) | 1.711 | 0.893 | -0.037 |
| publication bias | R3 hi | 392 | 176 (0.449) | 1.000 | 0.459 | +0.005 |
| publication bias | R3 lo | 399 | 184 (0.461) | 1.350 | 0.920 | -0.049 |
| heterogeneous populations | R1 hi | 194 | 180 (0.928) | 1.000 | 0.000 | -0.000 |
| heterogeneous populations | R1 lo | 387 | 282 (0.729) | 1.208 | 0.633 | -0.006 |
| heterogeneous populations | R2 hi | 324 | 228 (0.704) | 1.000 | 0.386 | +0.010 |
| heterogeneous populations | R2 lo | 392 | 259 (0.661) | 1.696 | 0.901 | -0.038 |
| heterogeneous populations | R3 hi | 389 | 158 (0.406) | 1.000 | 0.434 | -0.000 |
| heterogeneous populations | R3 lo | 399 | 158 (0.396) | 1.338 | 0.930 | -0.050 |
| missing controls and variances | R1 hi | 221 | 197 (0.891) | 1.000 | 0.000 | -0.000 |
| missing controls and variances | R1 lo | 392 | 283 (0.722) | 1.214 | 0.630 | -0.003 |
| missing controls and variances | R2 hi | 221 | 197 (0.891) | 1.000 | 0.000 | -0.000 |
| missing controls and variances | R2 lo | 392 | 283 (0.722) | 1.214 | 0.630 | -0.003 |
| missing controls and variances | R3 hi | 377 | 239 (0.634) | 1.028 | 0.469 | +0.002 |
| missing controls and variances | R3 lo | 397 | 247 (0.622) | 1.420 | 0.929 | -0.041 |
| contradictory evidence | R1 hi | 211 | 193 (0.915) | 1.000 | 0.000 | -0.000 |
| contradictory evidence | R1 lo | 400 | 304 (0.760) | 1.411 | 0.885 | -0.154 |
| contradictory evidence | R2 hi | 346 | 242 (0.699) | 1.000 | 0.465 | +0.064 |
| contradictory evidence | R2 lo | 400 | 280 (0.700) | 1.822 | 0.968 | -0.114 |
| contradictory evidence | R3 hi | 400 | 172 (0.430) | 1.034 | 0.472 | +0.007 |
| contradictory evidence | R3 lo | 400 | 172 (0.430) | 1.371 | 0.980 | -0.023 |

## Real-data series

The control outcome of the hormone-therapy question is omitted: at rungs 2 and 3 it is bounded by assumption.

| Question | Rung | Side | Cells measuring the target | Most precise | Zero weight on it | Weights | SD under the rule | Least SD |
|---|---|---|---|---|---|---|---|---|
| hormone therapy, cardiovascular disease | R1 | lower | 1 | A_cvd | no | A_cvd: +1.000 | 0.0630 | 0.0630 |
| hormone therapy, cardiovascular disease | R2 | lower | 1 | A_cvd | no | A_cvd: +1.000 | 0.0630 | 0.0630 |
| hormone therapy, cardiovascular disease | R3 | upper | 1 | A_cvd | no | A_cvd: +1.000, A_cancer: -1.000 | 0.0994 | 0.0994 |
| hormone therapy, cardiovascular disease | R3 | lower | 1 | A_cvd | no | A_cvd: +1.000, A_cancer: -1.000 | 0.0994 | 0.0630 |
| sigmoidoscopy, distal cancer death | R1 | lower | 1 | B_distal | no | B_distal: +1.000 | 0.2590 | 0.2590 |
| sigmoidoscopy, distal cancer death | R2 | upper | 1 | B_distal | no | B_distal: +1.000, B_proximal: -1.000 | 0.3461 | 0.3461 |
| sigmoidoscopy, distal cancer death | R2 | lower | 1 | B_distal | no | B_distal: +1.000, B_proximal: -1.000 | 0.3461 | 0.2590 |
| vitamin E, coronary disease | R1 | lower | 2 | C_hpfs | no | C_nhs: +0.307, C_hpfs: +0.693 | 0.1234 | 0.1234 |
| beta-carotene, frozen encoding | R1 | upper | 3 | D_rct | no | D_rct: +1.000 | 0.0077 | 0.0077 |
| beta-carotene, frozen encoding | R1 | lower | 3 | D_rct | no | D_rct: +1.000 | 0.0077 | 0.0077 |
| beta-carotene, frozen encoding | R2 | upper | 3 | D_rct | no | D_rct: +0.895, D_mr: +0.105 | 0.0073 | 0.0073 |
| beta-carotene, frozen encoding | R2 | lower | 3 | D_rct | no | D_rct: +0.895, D_mr: +0.105 | 0.0073 | 0.0072 |
| beta-carotene, trial validity named | R2 | upper | 3 | D_rct | yes | D_mr: +1.000 | 0.0224 | 0.0224 |
| beta-carotene, trial validity named | R2 | lower | 3 | D_rct | yes | D_mr: +1.000 | 0.0224 | 0.0201 |
| beta-carotene, trial validity named | R3 | upper | 3 | D_rct | no | D_rct: +0.895, D_mr: +0.105 | 0.0073 | 0.0073 |
| beta-carotene, trial validity named | R3 | lower | 3 | D_rct | no | D_rct: +0.895, D_mr: +0.105 | 0.0073 | 0.0072 |
