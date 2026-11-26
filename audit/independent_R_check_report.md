# Independent Reviewer Report

**Manuscript:** SIGNET  
**Review type:** Independent methodological, statistical, and computational reproducibility assessment  
**Independent computational environment:** R  
**Reviewer role:** Independent reviewer

## 1. Overall Assessment

The manuscript presents a methodological framework that warrants careful evaluation of its mathematical formulation, statistical properties, numerical implementation, and reproducibility.

As part of this independent review, I conducted computational assessments using R to examine the implementation of key methodological components and to evaluate whether the numerical results reported in the supplementary worked example could be independently reproduced.

The assessment included examination of the linear programming formulations, certificate calculations, constrained optimization procedures, implicit-constraint identification, degrees-of-freedom determination, and confidence-set inversion.

The central numerical finding is that the supplementary sigmoidoscopy worked example (S3) was successfully reproduced following correction of an issue in the independent R implementation of the implicit-constraint identification procedure. The resulting confidence limits were consistent with the reported supplementary values.

Overall, the computational assessment provides evidence supporting the numerical reproducibility of this worked example. However, reproduction of an individual example should not be interpreted as comprehensive validation of every theoretical result, simulation, or application presented in the manuscript.

## 2. Independent Computational Assessment Using R

I used R to independently reconstruct and evaluate the relevant computational procedures.

The assessment focused on the following components:

### 2.1. Linear programming and capacity calculations

The capacity calculations were assessed using the linear programming formulation described in the supplementary methods.

The assessment examined whether the relevant optimization problems were feasible and whether the expected boundedness patterns were reproduced under the specified assumptions.

The capacity pattern for the sigmoidoscopy worked example was successfully reproduced.

### 2.2. Certificate calculations

The certificate optimization procedures were independently examined, including the feasibility of the certificate constraints and the corresponding certificate weights.

The certificate calculations were consistent with the expected algebraic structure of the worked example.

The certificate variance calculations were also examined as part of the numerical reconstruction.

### 2.3. Constrained optimization and test-statistic calculation

The constrained fitting procedure was evaluated across candidate values of the target parameter.

Particular attention was given to the identification of implicit inequality constraints and their contribution to the degrees-of-freedom calculation.

An important implementation issue was identified in the initial independent R reconstruction: the implicit-constraint identification procedure did not correctly classify constraints that could not be made strictly negative under the relevant equality and inequality restrictions.

The procedure was subsequently corrected using a linear programming minimization criterion. A constraint was classified as implicit when its minimum achievable value was not strictly negative, within the prespecified numerical tolerance.

This correction resulted in the appropriate identification of both implicit constraints at the lower confidence boundary of the worked example.

The degrees of freedom were then recalculated separately at each tested parameter value, consistent with the parameter-dependent formulation.

### 2.4. Confidence-set inversion

The confidence-set inversion was independently recalculated using R.

The analysis included:

- Evaluation of the test statistic over a 40-point parameter grid.
- Identification of accepted and rejected parameter values.
- Numerical bisection to determine the lower and upper confidence boundaries.
- Recalculation of the degrees of freedom at each tested parameter value.
- Local validation of the acceptance/rejection behavior around both confidence boundaries.

## 3. Reproduction of the Supplementary Sigmoidoscopy Worked Example

The corrected R implementation produced the following confidence limits:

| Parameter | Independent R result | Supplementary reported result |
|---|---:|---:|
| Lower risk-ratio limit | 0.20320 | Approximately 0.203 |
| Upper risk-ratio limit | 0.92763 | Approximately 0.928 |
| Lower log-risk-ratio limit | −1.59355 | Approximately −1.595 |
| Upper log-risk-ratio limit | −0.07512 | Approximately −0.075 |

The numerical agreement was satisfactory.

The boundary-specific degrees of freedom were:

| Confidence boundary | Implicit constraints | Degrees of freedom |
|---|---|---:|
| Lower | Both inequality constraints | 2 |
| Upper | None | 1 |

This finding illustrates the importance of correctly identifying implicit constraints at each tested parameter value rather than relying on the degrees of freedom calculated at the fitted center.

## 4. Boundary Validation

To assess whether the numerical inversion correctly identified the confidence-set boundaries, I evaluated six parameter values immediately around and at the estimated endpoints.

The acceptance threshold was \(p=0.025\).

| Evaluation point | P-value | Expected classification | Observed classification |
|---|---:|---|---|
| Just below lower limit | 0.024995 | Rejected | Rejected |
| At lower limit | 0.025001 | Accepted | Accepted |
| Just above lower limit | 0.025006 | Accepted | Accepted |
| Just below upper limit | 0.025006 | Accepted | Accepted |
| At upper limit | 0.025002 | Accepted | Accepted |
| Just above upper limit | 0.024998 | Rejected | Rejected |

All six evaluations demonstrated the expected local acceptance/rejection pattern.

The numerical confidence-set inversion therefore passed the local boundary validation.

## 5. Interpretation of the Computational Findings

The independent R assessment supports the following conclusions:

1. The capacity pattern for the worked example was reproduced.
2. The certificate calculations were consistent with the expected mathematical structure.
3. Correct identification of implicit inequality constraints was necessary to reproduce the reported lower confidence limit.
4. The degrees of freedom varied across the parameter space and were correctly recalculated in the corrected implementation.
5. The lower and upper confidence limits agreed with the supplementary worked example to the reported numerical precision.
6. The local boundary validation demonstrated the expected transition between accepted and rejected parameter values.

An important distinction is that the implementation issue described above was identified in my independent R reconstruction. It should not be interpreted as evidence of an error in the authors' software unless the corresponding procedure in the authors' original implementation is independently shown to contain the same issue.

## 6. Limitations of This Assessment

Although the numerical reproduction of the S3 worked example was successful, several limitations should be acknowledged.

First, reproducing a worked example does not independently establish all theoretical properties of the proposed method.

Second, the present R assessment does not, by itself, establish complete reproduction of every simulation scenario, sensitivity analysis, real-world application, or numerical experiment in the manuscript.

Third, numerical agreement should be distinguished from theoretical validity. The latter requires assessment of the assumptions, mathematical derivations, identification conditions, statistical guarantees, and interpretation of the proposed framework.

Finally, the local boundary checks establish numerical behavior around the two estimated endpoints but do not independently establish that the entire confidence set is connected or that no additional accepted regions exist outside the examined grid.

These considerations should be taken into account when interpreting the scope of the reproducibility findings.

## 7. Overall Reviewer Conclusion

The independent computational assessment using R successfully reproduced the confidence limits of the supplementary sigmoidoscopy worked example after correcting the implicit-constraint identification procedure in the independent implementation.

The final confidence set was approximately:

\[
\boxed{RR \in [0.2032,\;0.9276]}
\]

This agrees with the supplementary worked example, which reports rounded limits of approximately 0.203 to 0.928.

The findings provide encouraging evidence of numerical reproducibility for this specific example. They also demonstrate that correct implementation of the implicit-constraint and degrees-of-freedom procedures is essential to obtaining the reported confidence set.

I recommend that the authors ensure that the supplementary materials clearly explain the parameter-dependent identification of implicit constraints, the corresponding degrees-of-freedom calculation, and the numerical tolerances used in the implementation.

**Final assessment:** The S3 numerical worked example was independently reproduced using R, and all six local confidence-boundary validation checks passed. This conclusion is restricted to the computational components examined and should not be interpreted as comprehensive validation of the entire methodology.