# R interface to SIGNET

`signet.R` lets an analysis be run from R. The computation is done by the Python reference implementation in the parent folder; the R functions write the specification, call it and return R lists and data frames. Only the R package `jsonlite` is needed on the R side.

`signet.R` is an interface, not an independent implementation.

`reconstruction/signet_reconstruction.R` is a second implementation of the method in R that does not call the Python package (R packages jsonlite, lpSolve and quadprog). Run from the repository root, `Rscript r/reconstruction/signet_reconstruction.R` recomputes the four example specifications, compares 80 quantities with the results of the package and writes `reconstruction/reconstruction_report.md`. It was written on 22 November 2026 by the language model that revised the software, so it is a consistency check in a second language and not verification by an independent person (Supplementary S12.12).

The reconstruction of the worked example that versions 8 to 10 of the manuscript cited is not part of this archive: its report is `audit/independent_R_check_report.md`, its source code was never deposited, and it is no longer cited as evidence.

```r
# once, in a terminal at the repository root:
#   python3 -m pip install -r requirements-lock.txt && python3 -m pip install -e .
install.packages("jsonlite")
source("r/signet.R")                 # from the repository root

ex <- "examples/sigmoidoscopy_distal.json"
signet_validate(ex)                  # reports every error in the specification at once
signet_capacity(ex)                  # steps 5 and 6 of Figure 1, before extraction
res <- signet_run(ex, report = "sigmoidoscopy.md")
signet_summary(res)                  # sets, p-values, model test and the default
```

A specification can also be built as an R list with the fields documented in Supplementary S3 and passed to the same functions. If Python is not found, set `options(signet.python = "/path/to/python3")`; to use the repository without installing it, set `options(signet.repo = "/path/to/signet-evidence")`.

`Rscript r/tests/test_signet_r.R`, run from the repository root, reproduces the worked example of Figure 1 (7 checks).
