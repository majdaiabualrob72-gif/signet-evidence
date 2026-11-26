# Run from the repository root: Rscript r/tests/test_signet_r.R
source("r/signet.R")
options(signet.repo = ".")
ex <- "examples/sigmoidoscopy_distal.json"
ok <- TRUE
chk <- function(name, cond) { cat(if (cond) "PASS" else "FAIL", name, "\n"); ok <<- ok && cond }
chk("valid specification", isTRUE(signet_validate(ex)))
cap <- signet_capacity(ex)
chk("rung 1 cannot certify benefit", !cap$benefit_certifiable[cap$rung == "rung1"])
chk("rung 2 can certify benefit", cap$benefit_certifiable[cap$rung == "rung2"])
chk("certificate is distal minus proximal", grepl("B_distal \\+1.00; B_proximal -1.00", cap$benefit_certificate[cap$rung == "rung2"]))
s <- signet_summary(signet_run(ex))
r2 <- s[s$rung == "rung2", ]
chk("rung 2 set 0.203 to 0.928", r2$set == "0.203 to 0.928")
chk("benefit p about 0.028", abs(r2$p_benefit - 0.0279) < 1e-3)
chk("default unbounded", r2$default == "0 to inf")
cat(if (ok) "ALL PASS\n" else "SOME FAILED\n")
quit(status = if (ok) 0 else 1)
