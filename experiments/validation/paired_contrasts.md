# Method contrasts, paired (paired_contrasts.py)

paired_contrasts.py: paired re-analysis of the method contrasts (manuscript Table 4) and summary
of the v1.1.0 validation experiment.

Every contrast compares two methods on the same simulated literatures (same generator and seeds),
now using one record per literature: a = both succeed, b = only A succeeds, c = only B succeeds.
Difference in proportions (b - c) / n with the Newcombe paired-data interval (method 10, which
combines Wilson intervals for the two marginal proportions with a correlation correction) and
the exact two-sided McNemar test (binomial test on the b + c discordant literatures).

Inputs : experiments/simulation/benchmark_reps_{i}.json, benchmark_v2_stress_reps_{i}.json (archived)
         experiments/validation/reps_{part}.json (paired_reruns.py) and validation_reps_{i}.json
Outputs: experiments/validation/paired_contrasts.md, paired_contrasts.json, validation_summary.md

| Contrast | Method A, x/n (proportion) | Method B, x/n (proportion) | Discordant literatures (A only / B only) | Difference A minus B (95% CI) | p |
|---|---|---|---|---|---|
| Coverage, reference: SIGNET k=1 vs RE | 399/400 (0.998) | 280/400 (0.700) | 120 / 1 | +0.298 (+0.253 to +0.344) | < 0.0001 |
| Coverage, reference: SIGNET k=1 vs adaptive protocol | 399/400 (0.998) | 377/400 (0.943) | 23 / 1 | +0.055 (+0.033 to +0.082) | < 0.0001 |
| Coverage, 15% of signs wrong: SIGNET k=1 vs adaptive protocol | 396/400 (0.990) | 353/400 (0.882) | 43 / 0 | +0.108 (+0.079 to +0.142) | < 0.0001 |
| Coverage, heterogeneous populations: SIGNET k=1 vs adaptive protocol | 375/400 (0.938) | 312/400 (0.780) | 65 / 2 | +0.157 (+0.121 to +0.197) | < 0.0001 |
| Benefit certified, reference (true benefit): SIGNET k=1 vs BA-S | 268/400 (0.670) | 376/400 (0.940) | 2 / 110 | -0.270 (-0.316 to -0.225) | < 0.0001 |
| False benefit, null effect: SIGNET k=1 vs RE | 0/400 (0.000) | 125/400 (0.312) | 0 / 125 | -0.312 (-0.360 to -0.268) | < 0.0001 |
| False benefit, one-signed literature, doubled biases: SIGNET k=1 vs BA-S | 0/400 (0.000) | 41/400 (0.102) | 0 / 41 | -0.102 (-0.136 to -0.075) | < 0.0001 |
| False benefit, one-signed literature, doubled biases: SIGNET k=1 vs RE | 0/400 (0.000) | 316/400 (0.790) | 0 / 316 | -0.790 (-0.827 to -0.746) | < 0.0001 |
| Coverage, biases x4 (mixed designs): SIGNET k=1 vs BA-S | 199/200 (0.995) | 120/200 (0.600) | 79 / 0 | +0.395 (+0.326 to +0.464) | < 0.0001 |
| False benefit, biases x3 (one-signed): SIGNET k=1 vs BA-S | 0/200 (0.000) | 132/200 (0.660) | 0 / 132 | -0.660 (-0.722 to -0.589) | < 0.0001 |
| Coverage, contradictory evidence: SIGNET k=1 vs RE | 297/400 (0.743) | 373/400 (0.932) | 15 / 91 | -0.190 (-0.238 to -0.143) | < 0.0001 |
| Coverage, heterogeneous populations: SIGNET R3 with vs without allowance | 375/400 (0.938) | 260/400 (0.650) | 115 / 0 | +0.287 (+0.243 to +0.333) | < 0.0001 |
| Benefit certified, reference: SIGNET-B k=1 vs SIGNET k=1 | 275/400 (0.688) | 268/400 (0.670) | 9 / 2 | +0.017 (+0.000 to +0.035) | 0.0654 |
| Coverage, contradictory evidence: group-or-cell tolerance vs group tolerance | 185/200 (0.925) | 148/200 (0.740) | 37 / 0 | +0.185 (+0.132 to +0.243) | < 0.0001 |
| Coverage, contradictory evidence: group-or-cell tolerance vs RE | 185/200 (0.925) | 187/200 (0.935) | 8 / 10 | -0.010 (-0.056 to +0.035) | 0.8145 |
| Benefit certified, reference: group-or-cell tolerance vs group tolerance | 138/200 (0.690) | 138/200 (0.690) | 0 / 0 | +0.000 (-0.012 to +0.012) | 1.0000 |
| False benefit, null effect: group-or-cell tolerance vs RE | 0/200 (0.000) | 67/200 (0.335) | 0 / 67 | -0.335 (-0.403 to -0.270) | < 0.0001 |

## Reproduction of the archived (version 1.0.0) summaries by the re-runs

The missing-controls stress scenario was regenerated after a generator fix (Supplementary S11.12); its rows compare the regenerated values with the version 1.0.0 archive and are not expected to agree.

| Part | Scenario | Method | Re-run | Archived |
|---|---|---|---|---|
| comparators | reference | FE (all designs) | 0.3650 | 0.3650 |
| comparators | reference | BA-U (unsigned bias priors) | 0.9725 | 0.9725 |
| comparators | reference | BA-S (elicited directions and magnitudes) | 0.9925 | 0.9925 |
| comparators | reference | NCO-cal (calibrated, RE) | 0.9825 | 0.9825 |
| comparators | null effect | FE (all designs) | 0.3150 | 0.3150 |
| comparators | null effect | BA-U (unsigned bias priors) | 0.9700 | 0.9700 |
| comparators | null effect | BA-S (elicited directions and magnitudes) | 0.9900 | 0.9900 |
| comparators | null effect | NCO-cal (calibrated, RE) | 0.9650 | 0.9650 |
| comparators | magnitudes doubled | FE (all designs) | 0.1925 | 0.1925 |
| comparators | magnitudes doubled | BA-U (unsigned bias priors) | 0.7900 | 0.7900 |
| comparators | magnitudes doubled | BA-S (elicited directions and magnitudes) | 0.9525 | 0.9525 |
| comparators | magnitudes doubled | NCO-cal (calibrated, RE) | 0.9975 | 0.9975 |
| comparators | dense small biases | FE (all designs) | 0.2450 | 0.2450 |
| comparators | dense small biases | BA-U (unsigned bias priors) | 0.9550 | 0.9550 |
| comparators | dense small biases | BA-S (elicited directions and magnitudes) | 0.9825 | 0.9825 |
| comparators | dense small biases | NCO-cal (calibrated, RE) | 0.9825 | 0.9825 |
| comparators | 15% of signs wrong | FE (all designs) | 0.4400 | 0.4400 |
| comparators | 15% of signs wrong | BA-U (unsigned bias priors) | 0.9925 | 0.9925 |
| comparators | 15% of signs wrong | BA-S (elicited directions and magnitudes) | 0.9775 | 0.9775 |
| comparators | 15% of signs wrong | NCO-cal (calibrated, RE) | 0.9800 | 0.9800 |
| comparators | omitted bias | FE (all designs) | 0.6300 | 0.6300 |
| comparators | omitted bias | BA-U (unsigned bias priors) | 0.9925 | 0.9925 |
| comparators | omitted bias | BA-S (elicited directions and magnitudes) | 0.9725 | 0.9725 |
| comparators | omitted bias | NCO-cal (calibrated, RE) | 0.9800 | 0.9800 |
| comparators | no design opposing benefit | FE (all designs) | 0.2125 | 0.2125 |
| comparators | no design opposing benefit | BA-U (unsigned bias priors) | 0.9750 | 0.9750 |
| comparators | no design opposing benefit | BA-S (elicited directions and magnitudes) | 1.0000 | 1.0000 |
| comparators | no design opposing benefit | NCO-cal (calibrated, RE) | 0.9550 | 0.9550 |
| comparators | no bias | FE (all designs) | 0.9600 | 0.9600 |
| comparators | no bias | BA-U (unsigned bias priors) | 1.0000 | 1.0000 |
| comparators | no bias | BA-S (elicited directions and magnitudes) | 0.9875 | 0.9875 |
| comparators | no bias | NCO-cal (calibrated, RE) | 1.0000 | 1.0000 |
| comparators | one-signed literature, null effect | FE (all designs) | 0.0900 | 0.0900 |
| comparators | one-signed literature, null effect | BA-U (unsigned bias priors) | 0.9275 | 0.9275 |
| comparators | one-signed literature, null effect | BA-S (elicited directions and magnitudes) | 1.0000 | 1.0000 |
| comparators | one-signed literature, null effect | NCO-cal (calibrated, RE) | 0.7425 | 0.7425 |
| comparators | one-signed literature, null effect, magnitudes doubled | FE (all designs) | 0.0175 | 0.0175 |
| comparators | one-signed literature, null effect, magnitudes doubled | BA-U (unsigned bias priors) | 0.2700 | 0.2700 |
| comparators | one-signed literature, null effect, magnitudes doubled | BA-S (elicited directions and magnitudes) | 0.8975 | 0.8975 |
| comparators | one-signed literature, null effect, magnitudes doubled | NCO-cal (calibrated, RE) | 0.6675 | 0.6675 |
| cellrobust | extreme outlier | group or cell, coverage | 0.9900 | 0.9900 |
| cellrobust | extreme outlier | group or cell, benefit | 0.5900 | 0.5900 |
| cellrobust | publication bias | group or cell, coverage | 0.9950 | 0.9950 |
| cellrobust | publication bias | group or cell, benefit | 0.7200 | 0.7200 |
| cellrobust | heterogeneous populations | group or cell, coverage | 0.9450 | 0.9450 |
| cellrobust | heterogeneous populations | group or cell, benefit | 0.6200 | 0.6200 |
| cellrobust | missing controls and variances | group or cell, coverage | 1.0000 | 1.0000 |
| cellrobust | missing controls and variances | group or cell, benefit | 0.4500 | 0.4550 |
| cellrobust | contradictory evidence | group or cell, coverage | 0.9250 | 0.9250 |
| cellrobust | contradictory evidence | group or cell, benefit | 0.6050 | 0.6050 |
| cellrobust | reference | group or cell, coverage | 0.9950 | 0.9950 |
| cellrobust | reference | group or cell, benefit | 0.6900 | 0.6900 |
| cellrobust | null effect | group or cell, coverage | 0.9950 | 0.9950 |
| cellrobust | null effect | group or cell, benefit | 0.0000 | 0.0000 |
| cellrobust | magnitudes doubled | group or cell, coverage | 1.0000 | 1.0000 |
| cellrobust | magnitudes doubled | group or cell, benefit | 0.5900 | 0.5900 |
| cellrobust | dense small biases | group or cell, coverage | 1.0000 | 1.0000 |
| cellrobust | dense small biases | group or cell, benefit | 0.7200 | 0.7200 |
| cellrobust | 15% of signs wrong | group or cell, coverage | 0.9950 | 0.9950 |
| cellrobust | 15% of signs wrong | group or cell, benefit | 0.6150 | 0.6150 |
| cellrobust | omitted bias | group or cell, coverage | 1.0000 | 1.0000 |
| cellrobust | omitted bias | group or cell, benefit | 0.5400 | 0.5400 |
| cellrobust | no design opposing benefit | group or cell, coverage | 1.0000 | 1.0000 |
| cellrobust | no design opposing benefit | group or cell, benefit | 0.0600 | 0.0600 |
| cellrobust | no bias | group or cell, coverage | 1.0000 | 1.0000 |
| cellrobust | no bias | group or cell, benefit | 0.7500 | 0.7500 |
| cellrobust | one-signed literature, null effect | group or cell, coverage | 1.0000 | 1.0000 |
| cellrobust | one-signed literature, null effect | group or cell, benefit | 0.0000 | 0.0000 |
| cellrobust | one-signed literature, null effect, magnitudes doubled | group or cell, coverage | 1.0000 | 1.0000 |
| cellrobust | one-signed literature, null effect, magnitudes doubled | group or cell, benefit | 0.0000 | 0.0000 |
| sweep | mixed designs, true benefit | x4.0 | SIGNET R3 k=1 | 0.9950 | 0.9950 |
| sweep | mixed designs, true benefit | x4.0 | BA-S (elicited directions and magnitudes) | 0.6000 | 0.6000 |
| sweep | one-signed literature, null effect | x3.0 | SIGNET R3 k=1 | 1.0000 | 1.0000 |
| sweep | one-signed literature, null effect | x3.0 | BA-S (elicited directions and magnitudes) | 0.3400 | 0.3400 |
| hetallow | heterogeneous populations | R3 + allowance | 0.9375 | 0.9375 |
| bounded | reference | SIGNET-B k=1 benefit | 0.6875 | 0.6875 |
