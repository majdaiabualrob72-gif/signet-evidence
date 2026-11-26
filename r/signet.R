# R interface to the SIGNET reference implementation.
# The Python package does the computation; these functions write the specification,
# call the command-line tool and read its JSON result.

.signet_python <- function() {
  py <- getOption("signet.python", Sys.getenv("SIGNET_PYTHON", ""))
  if (!nzchar(py)) py <- Sys.which("python3")
  if (!nzchar(py)) py <- Sys.which("python")
  if (!nzchar(py)) stop("Python not found: set options(signet.python = '/path/to/python').")
  py
}

.signet_env <- function() {
  repo <- getOption("signet.repo", Sys.getenv("SIGNET_REPO", ""))
  if (nzchar(repo)) paste0("PYTHONPATH=", normalizePath(repo)) else character(0)
}

.spec_file <- function(spec) {
  if (is.character(spec) && length(spec) == 1 && file.exists(spec)) return(spec)
  f <- tempfile(fileext = ".json")
  jsonlite::write_json(spec, f, auto_unbox = TRUE, digits = NA, pretty = TRUE)
  f
}

.call <- function(args) {
  out <- suppressWarnings(system2(.signet_python(), c("-m", "signet", args),
                                  stdout = TRUE, stderr = TRUE, env = .signet_env()))
  status <- attr(out, "status")
  list(out = out, status = if (is.null(status)) 0L else status)
}

#' Validate a SIGNET specification (an R list or a JSON file); stops with every error found.
signet_validate <- function(spec) {
  r <- .call(c("validate", .spec_file(spec)))
  if (r$status != 0) stop(paste(r$out, collapse = "\n"))
  invisible(TRUE)
}

#' Run a specification. Returns the result as a list; optionally writes the Markdown report.
signet_run <- function(spec, report = NULL) {
  out <- tempfile(fileext = ".json")
  args <- c("run", .spec_file(spec), "--out", out)
  if (!is.null(report)) args <- c(args, "--report", report)
  r <- .call(args)
  if (r$status != 0) stop(paste(r$out, collapse = "\n"))
  jsonlite::fromJSON(out, simplifyVector = FALSE)
}

#' Capacity and certificate weights before extraction (steps 5 and 6 of Figure 1).
#' Estimates are replaced by 1.0 with the expected confidence limits, because neither
#' capacity nor the certificate weights depend on the estimates.
signet_capacity <- function(spec) {
  s <- if (is.character(spec)) jsonlite::fromJSON(spec, simplifyVector = FALSE) else spec
  for (i in seq_along(s$cells)) {
    c <- s$cells[[i]]
    if (is.null(c$scale) || c$scale == "ratio") {
      w <- if (!is.null(c$lower) && !is.null(c$upper)) sqrt(c$upper / c$lower) else exp(1.96 * c$se)
      s$cells[[i]]$estimate <- 1.0; s$cells[[i]]$lower <- 1 / w; s$cells[[i]]$upper <- w
    } else {
      s$cells[[i]]$estimate <- 0.0
    }
  }
  res <- signet_run(s)
  rows <- list()
  for (tg in names(res$targets)) for (rg in names(res$targets[[tg]])) {
    rec <- res$targets[[tg]][[rg]]
    cert <- rec$claims$benefit$certificate
    rows[[length(rows) + 1]] <- data.frame(
      target = tg, rung = rg,
      benefit_certifiable = isTRUE(rec$capacity$above),
      harm_certifiable = isTRUE(rec$capacity$below),
      benefit_certificate = if (is.null(cert)) NA_character_ else
        paste(sprintf("%s %+.2f", names(cert$weights), unlist(cert$weights)), collapse = "; "),
      assumptions_used = if (is.null(cert)) NA_character_ else paste(unlist(cert$assumptions), collapse = "; "),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

#' One row per target and rung: sets, p-values, model test and the default (primary result).
signet_summary <- function(result) {
  g <- function(x, k) if (is.null(x) || is.null(x[[k]])) NA else x[[k]]
  rows <- list()
  for (tg in names(result$targets)) for (rg in names(result$targets[[tg]])) {
    rec <- result$targets[[tg]][[rg]]
    rows[[length(rows) + 1]] <- data.frame(
      target = tg, rung = rg,
      set = g(rec$set, "ratio_scale"),
      tolerant_k1 = g(rec$tolerance_k1, "ratio_scale"),
      p_benefit = g(rec$claims$benefit, "p"), p_harm = g(rec$claims$harm, "p"),
      model_test_p = g(rec$model_test, "p"),
      default = g(rec$default, "ratio_scale"),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}
