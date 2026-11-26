# Reconstruction of SIGNET in R: comparison with the Python package

R 4.3.3; lpSolve 5.6.20; quadprog 1.5.8; jsonlite 1.8.8. Run on 2026-11-22 16:04 +03.

Comparisons: 80. Outside tolerance: 0. Largest absolute difference among finite values: 7.00e-05.
Tolerances: 0 for capacity, 5e-4 on the log scale for set limits (the bisection tolerances are 1e-4 in the package and 1e-5 here), 1e-6 for p-values.

Worked example (sigmoidoscopy, rung 2): test-based limits -1.59356 and -0.07510 on the log scale (RR 0.2032 to 0.9277); just inside the lower limit T = 7.373 with 2 degrees of freedom, p = 0.0251; just inside the upper limit T = 5.021 with 1 degree of freedom, p = 0.0250.

| Specification | Quantity | R | Python package | Absolute difference | Within tolerance |
|---|---|---|---|---|---|
| sigmoidoscopy_distal | rung 1 capacity above | 0 | 0 | 0.0e+00 | yes |
| sigmoidoscopy_distal | rung 1 capacity below | 1 | 1 | 0.0e+00 | yes |
| sigmoidoscopy_distal | rung 1 set, lower limit (log) | -1.47210 (RR 0.2294) | -1.47208 (RR 0.2294) | 1.8e-05 | yes |
| sigmoidoscopy_distal | rung 1 set, upper limit (log) | Inf (RR unbounded) | Inf (RR unbounded) | 0.0e+00 | yes |
| sigmoidoscopy_distal | rung 1 p for benefit | 1 | 1 | 0.0e+00 | yes |
| sigmoidoscopy_distal | rung 1 p for harm | 1 | 1 | 0.0e+00 | yes |
| sigmoidoscopy_distal | rung 1 group-tolerant set, lower (log) | -Inf (RR 0) | -Inf (RR 0) | 0.0e+00 | yes |
| sigmoidoscopy_distal | rung 1 group-tolerant set, upper (log) | Inf (RR unbounded) | Inf (RR unbounded) | 0.0e+00 | yes |
| sigmoidoscopy_distal | rung 2 capacity above | 1 | 1 | 0.0e+00 | yes |
| sigmoidoscopy_distal | rung 2 capacity below | 1 | 1 | 0.0e+00 | yes |
| sigmoidoscopy_distal | rung 2 set, lower limit (log) | -1.59356 (RR 0.2032) | -1.59356 (RR 0.2032) | 3.3e-06 | yes |
| sigmoidoscopy_distal | rung 2 set, upper limit (log) | -0.07510 (RR 0.9277) | -0.07517 (RR 0.9276) | 7.0e-05 | yes |
| sigmoidoscopy_distal | rung 2 p for benefit | 0.0279112 | 0.0279112 | 1.1e-16 | yes |
| sigmoidoscopy_distal | rung 2 p for harm | 1 | 1 | 0.0e+00 | yes |
| sigmoidoscopy_distal | rung 2 group-tolerant set, lower (log) | -1.62645 (RR 0.1966) | -1.62642 (RR 0.1966) | 2.8e-05 | yes |
| sigmoidoscopy_distal | rung 2 group-tolerant set, upper (log) | Inf (RR unbounded) | Inf (RR unbounded) | 0.0e+00 | yes |
| sigmoidoscopy_distal | rung 2 default, lower (log) | -Inf (RR 0) | -Inf (RR 0) | 0.0e+00 | yes |
| sigmoidoscopy_distal | rung 2 default, upper (log) | Inf (RR unbounded) | Inf (RR unbounded) | 0.0e+00 | yes |
| hormone_therapy_cvd | rung 1 capacity above | 0 | 0 | 0.0e+00 | yes |
| hormone_therapy_cvd | rung 1 capacity below | 1 | 1 | 0.0e+00 | yes |
| hormone_therapy_cvd | rung 1 set, lower limit (log) | -0.70327 (RR 0.4950) | -0.70325 (RR 0.4950) | 1.9e-05 | yes |
| hormone_therapy_cvd | rung 1 set, upper limit (log) | Inf (RR unbounded) | Inf (RR unbounded) | 0.0e+00 | yes |
| hormone_therapy_cvd | rung 1 p for benefit | 1 | 1 | 0.0e+00 | yes |
| hormone_therapy_cvd | rung 1 p for harm | 1 | 1 | 0.0e+00 | yes |
| hormone_therapy_cvd | rung 1 group-tolerant set, lower (log) | -Inf (RR 0) | -Inf (RR 0) | 0.0e+00 | yes |
| hormone_therapy_cvd | rung 1 group-tolerant set, upper (log) | Inf (RR unbounded) | Inf (RR unbounded) | 0.0e+00 | yes |
| hormone_therapy_cvd | rung 2 capacity above | 0 | 0 | 0.0e+00 | yes |
| hormone_therapy_cvd | rung 2 capacity below | 1 | 1 | 0.0e+00 | yes |
| hormone_therapy_cvd | rung 2 set, lower limit (log) | -0.59867 (RR 0.5495) | -0.59866 (RR 0.5495) | 9.6e-06 | yes |
| hormone_therapy_cvd | rung 2 set, upper limit (log) | Inf (RR unbounded) | Inf (RR unbounded) | 0.0e+00 | yes |
| hormone_therapy_cvd | rung 2 p for benefit | 1 | 1 | 0.0e+00 | yes |
| hormone_therapy_cvd | rung 2 p for harm | 1 | 1 | 0.0e+00 | yes |
| hormone_therapy_cvd | rung 2 group-tolerant set, lower (log) | -0.70327 (RR 0.4950) | -0.70327 (RR 0.4950) | 1.1e-06 | yes |
| hormone_therapy_cvd | rung 2 group-tolerant set, upper (log) | Inf (RR unbounded) | Inf (RR unbounded) | 0.0e+00 | yes |
| hormone_therapy_cvd | rung 3 capacity above | 1 | 1 | 0.0e+00 | yes |
| hormone_therapy_cvd | rung 3 capacity below | 1 | 1 | 0.0e+00 | yes |
| hormone_therapy_cvd | rung 3 set, lower limit (log) | -0.59867 (RR 0.5495) | -0.59867 (RR 0.5495) | 2.9e-06 | yes |
| hormone_therapy_cvd | rung 3 set, upper limit (log) | -0.15291 (RR 0.8582) | -0.15291 (RR 0.8582) | 2.9e-06 | yes |
| hormone_therapy_cvd | rung 3 p for benefit | 0.000314871 | 0.000314871 | 1.1e-19 | yes |
| hormone_therapy_cvd | rung 3 p for harm | 1 | 1 | 0.0e+00 | yes |
| hormone_therapy_cvd | rung 3 group-tolerant set, lower (log) | -0.70327 (RR 0.4950) | -0.70322 (RR 0.4950) | 4.6e-05 | yes |
| hormone_therapy_cvd | rung 3 group-tolerant set, upper (log) | Inf (RR unbounded) | Inf (RR unbounded) | 0.0e+00 | yes |
| hormone_therapy_cvd | rung 3 default, lower (log) | -Inf (RR 0) | -Inf (RR 0) | 0.0e+00 | yes |
| hormone_therapy_cvd | rung 3 default, upper (log) | Inf (RR unbounded) | Inf (RR unbounded) | 0.0e+00 | yes |
| vitamin_e_chd | rung 1 capacity above | 0 | 0 | 0.0e+00 | yes |
| vitamin_e_chd | rung 1 capacity below | 1 | 1 | 0.0e+00 | yes |
| vitamin_e_chd | rung 1 set, lower limit (log) | -0.75863 (RR 0.4683) | -0.75863 (RR 0.4683) | 5.4e-10 | yes |
| vitamin_e_chd | rung 1 set, upper limit (log) | Inf (RR unbounded) | Inf (RR unbounded) | 0.0e+00 | yes |
| vitamin_e_chd | rung 1 p for benefit | 1 | 1 | 0.0e+00 | yes |
| vitamin_e_chd | rung 1 p for harm | 1 | 1 | 0.0e+00 | yes |
| vitamin_e_chd | rung 1 group-tolerant set, lower (log) | -Inf (RR 0) | -Inf (RR 0) | 0.0e+00 | yes |
| vitamin_e_chd | rung 1 group-tolerant set, upper (log) | Inf (RR unbounded) | Inf (RR unbounded) | 0.0e+00 | yes |
| vitamin_e_chd | rung 1 default, lower (log) | -Inf (RR 0) | -Inf (RR 0) | 0.0e+00 | yes |
| vitamin_e_chd | rung 1 default, upper (log) | Inf (RR unbounded) | Inf (RR unbounded) | 0.0e+00 | yes |
| beta_carotene_chd_trial_validity_named | rung 1 capacity above | 0 | 0 | 0.0e+00 | yes |
| beta_carotene_chd_trial_validity_named | rung 1 capacity below | 0 | 0 | 0.0e+00 | yes |
| beta_carotene_chd_trial_validity_named | rung 1 set, lower limit (log) | -Inf (RR 0) | -Inf (RR 0) | 0.0e+00 | yes |
| beta_carotene_chd_trial_validity_named | rung 1 set, upper limit (log) | Inf (RR unbounded) | Inf (RR unbounded) | 0.0e+00 | yes |
| beta_carotene_chd_trial_validity_named | rung 1 p for benefit | 1 | 1 | 0.0e+00 | yes |
| beta_carotene_chd_trial_validity_named | rung 1 p for harm | 1 | 1 | 0.0e+00 | yes |
| beta_carotene_chd_trial_validity_named | rung 1 group-tolerant set, lower (log) | -Inf (RR 0) | -Inf (RR 0) | 0.0e+00 | yes |
| beta_carotene_chd_trial_validity_named | rung 1 group-tolerant set, upper (log) | Inf (RR unbounded) | Inf (RR unbounded) | 0.0e+00 | yes |
| beta_carotene_chd_trial_validity_named | rung 2 capacity above | 1 | 1 | 0.0e+00 | yes |
| beta_carotene_chd_trial_validity_named | rung 2 capacity below | 1 | 1 | 0.0e+00 | yes |
| beta_carotene_chd_trial_validity_named | rung 2 set, lower limit (log) | -0.03043 (RR 0.9700) | -0.03040 (RR 0.9701) | 3.0e-05 | yes |
| beta_carotene_chd_trial_validity_named | rung 2 set, upper limit (log) | 0.07004 (RR 1.0726) | 0.07001 (RR 1.0725) | 3.0e-05 | yes |
| beta_carotene_chd_trial_validity_named | rung 2 p for benefit | 1 | 1 | 0.0e+00 | yes |
| beta_carotene_chd_trial_validity_named | rung 2 p for harm | 0.75394 | 0.75394 | 5.6e-16 | yes |
| beta_carotene_chd_trial_validity_named | rung 2 group-tolerant set, lower (log) | -0.28932 (RR 0.7488) | -0.28927 (RR 0.7488) | 4.8e-05 | yes |
| beta_carotene_chd_trial_validity_named | rung 2 group-tolerant set, upper (log) | Inf (RR unbounded) | Inf (RR unbounded) | 0.0e+00 | yes |
| beta_carotene_chd_trial_validity_named | rung 3 capacity above | 1 | 1 | 0.0e+00 | yes |
| beta_carotene_chd_trial_validity_named | rung 3 capacity below | 1 | 1 | 0.0e+00 | yes |
| beta_carotene_chd_trial_validity_named | rung 3 set, lower limit (log) | -0.01422 (RR 0.9859) | -0.01422 (RR 0.9859) | 2.7e-09 | yes |
| beta_carotene_chd_trial_validity_named | rung 3 set, upper limit (log) | 0.01839 (RR 1.0186) | 0.01839 (RR 1.0186) | 2.7e-09 | yes |
| beta_carotene_chd_trial_validity_named | rung 3 p for benefit | 1 | 1 | 0.0e+00 | yes |
| beta_carotene_chd_trial_validity_named | rung 3 p for harm | 1 | 1 | 0.0e+00 | yes |
| beta_carotene_chd_trial_validity_named | rung 3 group-tolerant set, lower (log) | -0.03043 (RR 0.9700) | -0.03040 (RR 0.9701) | 3.0e-05 | yes |
| beta_carotene_chd_trial_validity_named | rung 3 group-tolerant set, upper (log) | 0.07004 (RR 1.0726) | 0.07001 (RR 1.0725) | 3.0e-05 | yes |
| beta_carotene_chd_trial_validity_named | rung 3 default, lower (log) | -0.03043 (RR 0.9700) | -0.03040 (RR 0.9701) | 3.0e-05 | yes |
| beta_carotene_chd_trial_validity_named | rung 3 default, upper (log) | 0.07004 (RR 1.0726) | 0.07001 (RR 1.0725) | 3.0e-05 | yes |
