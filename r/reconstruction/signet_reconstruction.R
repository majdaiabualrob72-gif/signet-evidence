# signet_reconstruction.R
#
# A second implementation of SIGNET in R, written from the mathematical description in the manuscript
# (Sections 3.1 to 3.6 and Supplementary S1 and S2). It does not call the Python package, shares no code
# with it and uses different solvers (lpSolve for linear programmes, quadprog for quadratic programmes)
# and a different inversion (a fine scan of the outer interval instead of a coarse grid with a dense check).
#
# Written and run on 22 November 2026 by a language model (Claude, Anthropic) in the same session that
# revised the software to version 1.2.0. It is a reconstruction in a second language with an archived
# source, not a verification by an independent person.
#
# Usage, from the repository root:   Rscript r/reconstruction/signet_reconstruction.R
# It reads the four example specifications in examples/, recomputes capacity, certificates, the test
# statistic with its degrees of freedom, the confidence sets, the group-tolerant set, the default and the
# directional p-values, and compares them with the results written by the Python package
# (examples/*.result.json). Exit status 0 means every comparison was within tolerance.

suppressMessages({library(jsonlite); library(lpSolve); library(quadprog)})
Z975 <- qnorm(0.975)
TOL_T <- 1e-9; TOL_SLACK <- 1e-4; TOL_IMPL <- 1e-7

# ---------------------------------------------------------------- model
build <- function(spec) {
  field <- unlist(spec$field); G <- length(field)
  cells <- spec$cells; n <- length(cells)
  ids <- vapply(cells, function(c) c$id, "")
  inc <- list()
  for (i in seq_len(n)) for (p in unlist(cells[[i]]$processes)) inc[[length(inc) + 1]] <- list(p = p, i = i)
  P <- G + length(inc)
  pidx <- function(p, i) G + which(vapply(inc, function(z) z$p == p && z$i == i, TRUE))
  rows_with <- function(p) unlist(lapply(inc, function(z) if (z$p == p) z$i else NULL))
  A <- matrix(0, n, P); y <- numeric(n); se <- numeric(n)
  for (i in seq_len(n)) {
    A[i, seq_len(G)] <- unlist(cells[[i]]$kernel)
    y[i] <- log(cells[[i]]$estimate)
    se[i] <- (log(cells[[i]]$upper) - log(cells[[i]]$lower)) / (2 * Z975)
  }
  for (m in seq_along(inc)) A[inc[[m]]$i, G + m] <- 1
  pnames <- c(field, vapply(inc, function(z) paste0(z$p, "@", ids[z$i]), ""))
  groups <- list()
  for (a in spec$assumptions) {
    cons <- list()
    add <- function(r, type) cons[[length(cons) + 1]] <<- list(row = r, type = type)
    rows <- if (!is.null(a$rows)) match(unlist(a$rows), ids) else if (!is.null(a$process)) rows_with(a$process) else NULL
    if (a$type == "sign") {
      s <- if (a$direction == "negative") -1 else 1
      for (i in rows) { r <- numeric(P); r[pidx(a$process, i)] <- -s; add(r, "le") }      # s * b >= 0
    } else if (a$type == "zero") {
      for (i in rows) { r <- numeric(P); r[pidx(a$process, i)] <- 1; add(r, "eq") }
    } else if (a$type == "share") {
      j0 <- pidx(a$process, rows[1])
      for (i in rows[-1]) { r <- numeric(P); r[pidx(a$process, i)] <- 1; r[j0] <- -1; add(r, "eq") }
    } else if (a$type == "order") {
      r <- numeric(P); r[pidx(a$process, match(a$row_lo, ids))] <- 1; r[pidx(a$process, match(a$row_hi, ids))] <- -1; add(r, "le")
    } else if (a$type == "fieldsign") {
      s <- if (a$direction == "positive") 1 else -1
      r <- numeric(P); r[match(a$field, field)] <- -s; add(r, "le")                       # s * f >= 0
    } else if (a$type == "linear") {
      stopifnot(is.null(a$rhs) || a$rhs == 0)
      r <- numeric(P); for (nm in names(a$coefs)) r[match(nm, pnames)] <- a$coefs[[nm]]
      sense <- if (is.null(a$sense)) "le" else a$sense
      if (sense == "ge") { r <- -r; sense <- "le" }
      add(r, sense)
    } else stop("assumption type not implemented in this reconstruction: ", a$type)
    groups[[a$name]] <- list(rung = a$rung, cons = cons)
  }
  list(A = A, y = y, V = diag(se^2, n), n = n, P = P, G = G, field = field, ids = ids, groups = groups, pnames = pnames)
}

compile <- function(M, include) {
  D <- matrix(0, 0, M$P); E <- matrix(0, 0, M$P)
  for (g in include) for (cn in M$groups[[g]]$cons) {
    r <- cn$row / sqrt(sum(cn$row^2))                                                    # unit-length rows
    if (cn$type == "le") D <- rbind(D, r) else E <- rbind(E, r)
  }
  list(D = D, E = E)
}

cvec <- function(M, cf) c(cf, numeric(M$P - M$G))

# ---------------------------------------------------------------- linear algebra and solvers
nullspace <- function(Mx, tol = 1e-10) {
  p <- ncol(Mx)
  if (nrow(Mx) == 0) return(diag(p))
  s <- svd(Mx, nu = 0, nv = p)
  r <- sum(s$d > tol * max(1, s$d[1]))
  if (r >= p) matrix(0, p, 0) else s$v[, (r + 1):p, drop = FALSE]
}
pinv_solve <- function(Mx, b) {
  if (nrow(Mx) == 0) return(numeric(ncol(Mx)))
  s <- svd(Mx); keep <- s$d > 1e-10 * max(1, s$d[1])
  as.numeric(s$v[, keep, drop = FALSE] %*% ((t(s$u[, keep, drop = FALSE]) %*% b) / s$d[keep]))
}
rank_of <- function(Mx) { if (length(Mx) == 0) return(0L); d <- svd(Mx, nu = 0, nv = 0)$d; sum(d > 1e-9 * max(1, d[1])) }

# LP with free variables: optimise obj'x subject to Aeq x = beq, Aub x <= bub (variables split into positive and negative parts)
lp_free <- function(obj, Aeq = NULL, beq = NULL, Aub = NULL, bub = NULL, maximize = FALSE, nonneg = integer(0)) {
  p <- length(obj)
  free <- setdiff(seq_len(p), nonneg)
  Tm <- cbind(diag(p), -diag(p)[, free, drop = FALSE])                                   # x = T z, z >= 0
  con <- rbind(if (!is.null(Aeq) && nrow(Aeq)) Aeq %*% Tm, if (!is.null(Aub) && nrow(Aub)) Aub %*% Tm)
  dir <- c(rep("=", if (is.null(Aeq)) 0 else nrow(Aeq)), rep("<=", if (is.null(Aub)) 0 else nrow(Aub)))
  rhs <- c(beq, bub)
  if (is.null(con)) { con <- matrix(0, 1, ncol(Tm)); dir <- "<="; rhs <- 0 }
  r <- lp(if (maximize) "max" else "min", as.numeric(obj %*% Tm), con, dir, rhs)
  list(status = r$status, value = r$objval, x = if (r$status == 0) as.numeric(Tm %*% r$solution) else NULL)
}

# ---------------------------------------------------------------- capacity (Theorem 1)
capacity <- function(M, include, cf, keep = seq_len(M$n)) {
  cp <- compile(M, include); cc <- cvec(M, cf); A <- M$A[keep, , drop = FALSE]
  out <- c(above = NA, below = NA)
  for (side in c("above", "below")) {
    s <- if (side == "above") 1 else -1
    r <- lp_free(s * cc, rbind(A, cp$E), numeric(nrow(A) + nrow(cp$E)), rbind(cp$D, s * cc), c(numeric(nrow(cp$D)), 1), maximize = TRUE)
    out[side] <- (r$status == 0 && r$value < 1e-7)
  }
  out
}

# ---------------------------------------------------------------- certificates (Theorem 2 and the choice rule)
certificate <- function(M, include, cf, side, V, keep = seq_len(M$n)) {
  cp <- compile(M, include); A <- M$A[keep, , drop = FALSE]; n <- nrow(A)
  cc <- cvec(M, cf) * (if (side == "hi") 1 else -1)
  mD <- nrow(cp$D); mE <- nrow(cp$E); k <- n + mD + mE
  Aeq <- cbind(t(A), t(cp$D), t(cp$E))                                                   # A'lam + D'nu + E'kap = c
  cost <- c(numeric(n), rep(1, mD), numeric(mE))
  r1 <- lp_free(cost, Aeq, cc, nonneg = if (mD) n + seq_len(mD) else integer(0))
  if (r1$status != 0) return(NULL)                                                       # no certificate: side unbounded
  w <- r1$value
  z0 <- pinv_solve(Aeq, cc); N <- nullspace(Aeq)                                         # all solutions: z0 + N u
  Q <- matrix(0, k, k); Q[seq_len(n), seq_len(n)] <- 2 * V
  if (ncol(N) == 0) z <- z0 else {
    Dm <- t(N) %*% Q %*% N; Dm <- Dm + 1e-10 * max(1, sum(diag(Dm)) / ncol(N)) * diag(ncol(N))
    dv <- -as.numeric(t(N) %*% Q %*% z0)
    Am <- NULL; bv <- NULL
    if (mD) {                                                                             # nu >= 0 and sum(nu) <= w
      Sel <- cbind(matrix(0, mD, n), diag(mD), matrix(0, mD, mE))
      Am <- rbind(Sel %*% N, -matrix(cost, 1) %*% N)
      bv <- c(-as.numeric(Sel %*% z0), -(w * (1 + 1e-9) + 1e-9) + sum(cost * z0))
    }
    u <- if (is.null(Am)) solve(Dm, dv) else solve.QP(Dm, dv, t(Am), bv)$solution
    z <- z0 + as.numeric(N %*% u)
  }
  lam <- z[seq_len(n)]
  list(lam = lam, sd = sqrt(max(0, as.numeric(t(lam) %*% V %*% lam))), slack_weight = w)
}

# ---------------------------------------------------------------- test statistic, active rank and p-value (Section 3.5)
fit <- function(M, cp, y, W, A, cc = NULL, theta0 = NULL) {
  Meq <- if (is.null(cc)) cp$E else rbind(cp$E, cc)
  beq <- if (is.null(cc)) numeric(nrow(cp$E)) else c(numeric(nrow(cp$E)), theta0)
  if (nrow(Meq)) {
    x0 <- pinv_solve(Meq, beq)
    if (max(abs(Meq %*% x0 - beq)) > 1e-8) return(list(ok = FALSE, T = Inf))             # theta0 excluded by the equalities
    B <- nullspace(Meq)
  } else { x0 <- numeric(M$P); B <- diag(M$P) }
  r0 <- y - as.numeric(A %*% x0)
  if (ncol(B) == 0) {
    x <- x0
    if (nrow(cp$D) && any(cp$D %*% x > 1e-9)) return(list(ok = FALSE, T = Inf))
  } else {
    AB <- A %*% B
    Dm <- 2 * t(AB) %*% W %*% AB
    Dm <- Dm + 1e-9 * max(1, sum(diag(Dm)) / ncol(B)) * diag(ncol(B))
    dv <- 2 * as.numeric(t(AB) %*% W %*% r0)
    if (nrow(cp$D)) {
      s <- tryCatch(solve.QP(Dm, dv, t(-cp$D %*% B), as.numeric(cp$D %*% x0)), error = function(e) NULL)
      if (is.null(s)) return(list(ok = FALSE, T = Inf))
      wv <- s$solution
    } else wv <- solve(Dm, dv)
    x <- x0 + as.numeric(B %*% wv)
  }
  r <- y - as.numeric(A %*% x)
  list(ok = TRUE, T = max(0, as.numeric(t(r) %*% W %*% r)), x = x, mu = as.numeric(A %*% x),
       slack = if (nrow(cp$D)) -as.numeric(cp$D %*% x) else numeric(0), Meq = Meq, beq = beq)
}

active_rank <- function(M, cp, A, fr) {
  cand <- which(fr$slack < TOL_SLACK); impl <- integer(0)
  for (j in cand) {                                                                       # can constraint j be slack on the optimal face?
    Aeq <- cbind(rbind(A, fr$Meq), 0); beq <- c(fr$mu, fr$beq)
    Aub <- rbind(cbind(cp$D, 0), c(cp$D[j, ], 1), c(numeric(M$P), 1)); bub <- c(numeric(nrow(cp$D)), 0, 1)
    r <- lp_free(c(numeric(M$P), 1), Aeq, beq, Aub, bub, maximize = TRUE)
    if (r$status != 0 || r$value < TOL_IMPL) impl <- c(impl, j)
  }
  N <- nullspace(rbind(fr$Meq, cp$D[impl, , drop = FALSE]))
  if (ncol(N) == 0) return(nrow(A))
  nrow(A) - rank_of(A %*% N)
}

pval <- function(M, cp, y, W, A, cc, theta0) {
  fr <- fit(M, cp, y, W, A, cc, theta0)
  if (!fr$ok) return(list(p = 0, T = Inf, r = NA))
  if (fr$T <= TOL_T) return(list(p = 1, T = fr$T, r = 0))
  r <- active_rank(M, cp, A, fr)
  list(p = if (r == 0) 0 else pchisq(fr$T, r, lower.tail = FALSE), T = fr$T, r = r)
}

# ---------------------------------------------------------------- confidence set by inversion (fine scan of the outer interval)
conf_set <- function(M, include, cf, a, y = M$y, V = M$V, keep = seq_len(M$n), nscan = 400, tol = 1e-5) {
  cp <- compile(M, include); cc <- cvec(M, cf); A <- M$A[keep, , drop = FALSE]; W <- solve(V); n <- nrow(A)
  base <- fit(M, cp, y, W, A)
  if (!base$ok) return(list(empty = TRUE))
  th0 <- sum(cc * base$x); crit_n <- qchisq(1 - a, n); cap <- capacity(M, include, cf, keep)
  Tof <- function(t) fit(M, cp, y, W, A, cc, t)$T
  acc <- function(t) pval(M, cp, y, W, A, cc, t)$p >= a
  lim <- c(lo = NA, hi = NA)
  for (side in c("lo", "hi")) {
    s <- if (side == "hi") 1 else -1
    if (!cap[if (side == "hi") "above" else "below"]) { lim[side] <- s * Inf; next }
    step <- 0.05; t_in <- th0; t_out <- th0 + s * step
    while (Tof(t_out) <= crit_n && step < 1e4) { t_in <- t_out; step <- 2 * step; t_out <- th0 + s * step }
    for (it in 1:80) { mid <- (t_in + t_out) / 2; if (Tof(mid) <= crit_n) t_in <- mid else t_out <- mid; if (abs(t_out - t_in) < tol) break }
    lim[side] <- t_out
  }
  lo_s <- if (is.finite(lim["lo"])) lim["lo"] else th0 - 6; hi_s <- if (is.finite(lim["hi"])) lim["hi"] else th0 + 6
  grid <- sort(unique(c(seq(lo_s, hi_s, length.out = nscan), th0)))
  ok <- vapply(grid, acc, TRUE)
  if (!any(ok)) return(list(empty = TRUE, outer = lim))
  i1 <- min(which(ok)); i2 <- max(which(ok))
  refine <- function(t_acc, t_rej) { for (it in 1:60) { mid <- (t_acc + t_rej) / 2; if (acc(mid)) t_acc <- mid else t_rej <- mid; if (abs(t_rej - t_acc) < tol) break }; t_acc }
  lo <- if (!is.finite(lim["lo"]) && i1 == 1) -Inf else if (i1 > 1) refine(grid[i1], grid[i1 - 1]) else grid[1]
  hi <- if (!is.finite(lim["hi"]) && i2 == length(grid)) Inf else if (i2 < length(grid)) refine(grid[i2], grid[i2 + 1]) else grid[length(grid)]
  list(empty = FALSE, lo = as.numeric(lo), hi = as.numeric(hi), outer = lim, th0 = th0)
}

signet_set <- function(M, include, cf, alpha = 0.05, y = M$y, V = M$V, keep = seq_len(M$n)) {
  z <- qnorm(1 - alpha / 4)
  ch <- certificate(M, include, cf, "hi", V, keep); cl <- certificate(M, include, cf, "lo", V, keep)
  U <- if (is.null(ch)) Inf else sum(ch$lam * y) + z * ch$sd
  L <- if (is.null(cl)) -Inf else -(sum(cl$lam * y) + z * cl$sd)
  cs <- conf_set(M, include, cf, alpha / 2, y, V, keep)
  if (cs$empty) return(list(empty = TRUE))
  lo <- max(L, cs$lo); hi <- min(U, cs$hi)
  if (lo > hi) return(list(empty = TRUE))
  list(empty = FALSE, lo = lo, hi = hi, cert_lo = L, cert_hi = U, test_lo = cs$lo, test_hi = cs$hi, ch = ch, cl = cl)
}

union_sets <- function(sets) {
  ok <- Filter(function(s) !s$empty, sets)
  if (!length(ok)) return(list(empty = TRUE))
  list(empty = FALSE, lo = min(vapply(ok, function(s) s$lo, 0)), hi = max(vapply(ok, function(s) s$hi, 0)))
}
tolerant_set <- function(M, include, cf, alpha = 0.05) {
  union_sets(c(list(signet_set(M, include, cf, alpha)), lapply(include, function(g) signet_set(M, setdiff(include, g), cf, alpha))))
}
cell_tolerant_set <- function(M, include, cf, alpha = 0.05) {
  sets <- list(signet_set(M, include, cf, alpha))
  for (i in seq_len(M$n)) {
    keep <- setdiff(seq_len(M$n), i)
    if (!length(keep)) return(list(empty = FALSE, lo = -Inf, hi = Inf))
    sets[[length(sets) + 1]] <- signet_set(M, include, cf, alpha, M$y[keep], M$V[keep, keep, drop = FALSE], keep)
  }
  union_sets(sets)
}
default_set <- function(M, include, cf, alpha = 0.05) union_sets(list(tolerant_set(M, include, cf, alpha), cell_tolerant_set(M, include, cf, alpha)))

claim_p <- function(M, include, cf, claim, set_empty) {
  if (set_empty) return(NA)                                                              # no directional p-value for an empty set
  side <- if (claim == "benefit") "hi" else "lo"
  ct <- certificate(M, include, cf, side, M$V)
  p_cert <- if (is.null(ct) || ct$sd <= 0) 1 else pnorm(sum(ct$lam * M$y) / ct$sd)
  cp <- compile(M, include); cc <- cvec(M, cf); W <- solve(M$V)
  base <- fit(M, cp, M$y, W, M$A); th0 <- sum(cc * base$x)
  t_star <- if (claim == "benefit") max(0, th0) else min(0, th0)
  grid <- sort(unique(c(if (claim == "benefit") seq(0, 3, length.out = 61) else seq(-3, 0, length.out = 61), t_star)))
  p_T <- max(vapply(grid, function(t) pval(M, cp, M$y, W, M$A, cc, t)$p, 0))
  min(1, 4 * p_cert, 2 * p_T)
}

# ---------------------------------------------------------------- comparison with the Python package
num <- function(v) if (is.null(v)) NA else if (is.character(v)) (if (v == "inf") Inf else if (v == "-inf") -Inf else as.numeric(v)) else as.numeric(v)
fmt <- function(v) if (is.na(v)) "NA" else if (!is.finite(v)) (if (v > 0) "unbounded" else "0") else sprintf("%.4f", exp(v))
rows <- list(); worst <- 0; fails <- 0
cmp <- function(spec_name, what, r_val, py_val, tol) {
  d <- if (is.na(r_val) && is.na(py_val)) 0 else if (is.infinite(r_val) || is.infinite(py_val)) (if (identical(r_val, py_val)) 0 else Inf) else abs(r_val - py_val)
  ok <- !is.na(d) && d <= tol
  if (is.finite(d)) worst <<- max(worst, d)
  if (!ok) fails <<- fails + 1
  rows[[length(rows) + 1]] <<- data.frame(spec = spec_name, quantity = what, R = r_val, python = py_val, abs_diff = d, ok = ok, stringsAsFactors = FALSE)
}
specs <- c("sigmoidoscopy_distal", "hormone_therapy_cvd", "vitamin_e_chd", "beta_carotene_chd_trial_validity_named")
for (sn in specs) {
  spec <- fromJSON(file.path("examples", paste0(sn, ".json")), simplifyVector = FALSE)
  res <- fromJSON(file.path("examples", paste0(sn, ".result.json")), simplifyVector = FALSE)
  M <- build(spec)
  for (tg in spec$targets) {
    cf <- vapply(M$field, function(f) if (is.null(tg$coefficients[[f]])) 0 else tg$coefficients[[f]], 0)
    rungs <- unlist(spec$analysis$rungs)
    for (rg in rungs) {
      include <- names(M$groups)[vapply(M$groups, function(g) g$rung <= rg, TRUE)]
      py <- res$targets[[tg$name]][[paste0("rung", rg)]]
      cap <- capacity(M, include, cf)
      cmp(sn, sprintf("rung %d capacity above", rg), as.numeric(cap["above"]), as.numeric(py$capacity$above), 0)
      cmp(sn, sprintf("rung %d capacity below", rg), as.numeric(cap["below"]), as.numeric(py$capacity$below), 0)
      s <- signet_set(M, include, cf)
      cmp(sn, sprintf("rung %d set, lower limit (log)", rg), if (s$empty) NA else s$lo, num(py$set$lo), 5e-4)
      cmp(sn, sprintf("rung %d set, upper limit (log)", rg), if (s$empty) NA else s$hi, num(py$set$hi), 5e-4)
      for (cl in c("benefit", "harm")) cmp(sn, sprintf("rung %d p for %s", rg, cl), claim_p(M, include, cf, cl, s$empty), num(py$claims[[cl]]$p), 1e-6)
      t1 <- tolerant_set(M, include, cf)
      cmp(sn, sprintf("rung %d group-tolerant set, lower (log)", rg), if (t1$empty) NA else t1$lo, num(py$tolerance_k1$lo), 5e-4)
      cmp(sn, sprintf("rung %d group-tolerant set, upper (log)", rg), if (t1$empty) NA else t1$hi, num(py$tolerance_k1$hi), 5e-4)
      if (rg == max(rungs)) {
        d <- default_set(M, include, cf)
        cmp(sn, sprintf("rung %d default, lower (log)", rg), if (d$empty) NA else d$lo, num(py$default$lo), 5e-4)
        cmp(sn, sprintf("rung %d default, upper (log)", rg), if (d$empty) NA else d$hi, num(py$default$hi), 5e-4)
      }
    }
  }
}
tab <- do.call(rbind, rows)
# degrees of freedom at the limits of the worked example (Supplementary S2)
spec <- fromJSON("examples/sigmoidoscopy_distal.json", simplifyVector = FALSE); M <- build(spec)
inc <- names(M$groups); cp <- compile(M, inc); W <- solve(M$V); s <- signet_set(M, inc, 1)
eps <- 2e-4
dfl <- pval(M, cp, M$y, W, M$A, cvec(M, 1), s$test_lo + eps); dfu <- pval(M, cp, M$y, W, M$A, cvec(M, 1), s$test_hi - eps)
out <- c("# Reconstruction of SIGNET in R: comparison with the Python package", "",
         sprintf("R %s; lpSolve %s; quadprog %s; jsonlite %s. Run on %s.", paste(R.version$major, R.version$minor, sep = "."),
                 as.character(packageVersion("lpSolve")), as.character(packageVersion("quadprog")), as.character(packageVersion("jsonlite")), format(Sys.time(), "%Y-%m-%d %H:%M %Z")), "",
         sprintf("Comparisons: %d. Outside tolerance: %d. Largest absolute difference among finite values: %.2e.", nrow(tab), fails, worst),
         "Tolerances: 0 for capacity, 5e-4 on the log scale for set limits (the bisection tolerances are 1e-4 in the package and 1e-5 here), 1e-6 for p-values.", "",
         sprintf("Worked example (sigmoidoscopy, rung 2): test-based limits %.5f and %.5f on the log scale (RR %.4f to %.4f); just inside the lower limit T = %.3f with %d degrees of freedom, p = %.4f; just inside the upper limit T = %.3f with %d degree of freedom, p = %.4f.",
                 s$test_lo, s$test_hi, exp(s$test_lo), exp(s$test_hi), dfl$T, dfl$r, dfl$p, dfu$T, dfu$r, dfu$p), "",
         "| Specification | Quantity | R | Python package | Absolute difference | Within tolerance |", "|---|---|---|---|---|---|",
         sprintf("| %s | %s | %s | %s | %s | %s |", tab$spec, tab$quantity,
                 ifelse(grepl("log", tab$quantity), sprintf("%.5f (RR %s)", tab$R, vapply(tab$R, fmt, "")), sprintf("%.6g", tab$R)),
                 ifelse(grepl("log", tab$quantity), sprintf("%.5f (RR %s)", tab$python, vapply(tab$python, fmt, "")), sprintf("%.6g", tab$python)),
                 sprintf("%.1e", tab$abs_diff), ifelse(tab$ok, "yes", "NO")))
writeLines(out, "r/reconstruction/reconstruction_report.md")
cat(out[1:9], sep = "\n"); cat("\n")
if (fails) { print(tab[!tab$ok, ]); quit(status = 1) }
cat("ALL COMPARISONS WITHIN TOLERANCE\n")
