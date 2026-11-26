# Method contrasts, post hoc (contrasts.py)

Method contrasts reported in the manuscript (Section 7). Chosen after the results were seen (post hoc).

Each contrast is a difference between two proportions (coverage or rate of certified benefit),
estimated from archived summaries. The methods were applied to the same simulated literatures
(identical seeds and generator), but the summaries do not keep the pairing, so the samples are
treated as independent. This is conservative for positively correlated outcomes.
Interval: Newcombe hybrid score interval (method 10, Newcombe 1998). P-value: two-sided Fisher exact test.

| Contrast | Method A, x/n (proportion) | Method B, x/n (proportion) | Difference A minus B (95% CI) | p |
|---|---|---|---|---|
| Coverage, reference: SIGNET k=1 vs RE | 399/400 (0.998) | 280/400 (0.700) | +0.298 (+0.253 to +0.344) | < 0.0001 |
| Coverage, reference: SIGNET k=1 vs adaptive protocol | 399/400 (0.998) | 377/400 (0.943) | +0.055 (+0.033 to +0.082) | < 0.0001 |
| Coverage, 15% of signs wrong: SIGNET k=1 vs adaptive protocol | 396/400 (0.990) | 353/400 (0.882) | +0.108 (+0.076 to +0.143) | < 0.0001 |
| Coverage, heterogeneous populations: SIGNET k=1 vs adaptive protocol | 375/400 (0.938) | 312/400 (0.780) | +0.157 (+0.110 to +0.205) | < 0.0001 |
| Benefit certified, reference (true benefit): SIGNET k=1 vs BA-S | 268/400 (0.670) | 376/400 (0.940) | -0.270 (-0.321 to -0.218) | < 0.0001 |
| False benefit, null effect: SIGNET k=1 vs RE | 0/400 (0.000) | 125/400 (0.312) | -0.312 (-0.360 to -0.268) | < 0.0001 |
| False benefit, one-signed literature, doubled biases: SIGNET k=1 vs BA-S | 0/400 (0.000) | 41/400 (0.102) | -0.102 (-0.136 to -0.075) | < 0.0001 |
| False benefit, one-signed literature, doubled biases: SIGNET k=1 vs RE | 0/400 (0.000) | 316/400 (0.790) | -0.790 (-0.827 to -0.746) | < 0.0001 |
| Coverage, biases x4 (mixed designs): SIGNET k=1 vs BA-S | 199/200 (0.995) | 120/200 (0.600) | +0.395 (+0.326 to +0.464) | < 0.0001 |
| False benefit, biases x3 (one-signed): SIGNET k=1 vs BA-S | 0/200 (0.000) | 132/200 (0.660) | -0.660 (-0.722 to -0.589) | < 0.0001 |
| Coverage, contradictory evidence: SIGNET k=1 vs RE | 297/400 (0.743) | 373/400 (0.932) | -0.190 (-0.240 to -0.140) | < 0.0001 |
| Coverage, heterogeneous populations: SIGNET R3 with vs without allowance | 375/400 (0.938) | 260/400 (0.650) | +0.287 (+0.234 to +0.339) | < 0.0001 |
| Benefit certified, reference: SIGNET-B k=1 vs SIGNET k=1 | 275/400 (0.688) | 268/400 (0.670) | +0.017 (-0.047 to +0.082) | 0.6497 |
| Coverage, contradictory evidence: group-or-cell tolerance vs group tolerance | 185/200 (0.925) | 297/400 (0.743) | +0.182 (+0.122 to +0.236) | < 0.0001 |
| Coverage, contradictory evidence: group-or-cell tolerance vs RE | 185/200 (0.925) | 373/400 (0.932) | -0.007 (-0.057 to +0.033) | 0.7364 |
| Benefit certified, reference: group-or-cell tolerance vs group tolerance | 138/200 (0.690) | 268/400 (0.670) | +0.020 (-0.060 to +0.097) | 0.6444 |
