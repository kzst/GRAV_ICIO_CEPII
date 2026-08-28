# panel_granger_robustness_fdr_placebo.R
#
# Robustness diagnostics for panel_granger_causality() (see
# panel_granger_causality.R), addressing Reviewer 1, Points 5 & 6 of
# The_World_Economy_Rev.md: with dozens of sectors/products x explanatory
# variables x up to 5 lags tested at a nominal alpha = 0.01, the RAW
# number of "significant" pairs is expected to include a non-trivial
# number of false positives from multiple testing alone, on top of the
# short-panel (26-28 annual observations) concern. This script adds two
# standard, complementary diagnostics that do NOT change the existing
# pipeline, only add reporting on top of it:
#
#   (1) panel_granger_causality_pvalues() -- repeats the SAME
#       restricted-vs-unrestricted Wald-F-test logic as
#       panel_granger_causality(), but returns, for every (variable,
#       sector/product) pair, the actual minimum p-value across the
#       tested lags together with a within-pair Bonferroni correction for
#       the number of lags tried (a standard "min-p" procedure), instead
#       of discarding the p-value after the first threshold-crossing lag.
#
#   (2) fdr_correct_panel_results() -- applies Benjamini-Hochberg
#       false-discovery-rate correction (Benjamini & Hochberg, 1995,
#       J. R. Stat. Soc. B, 57(1), 289-300, doi:10.1111/j.2517-6161.1995.tb02031.x)
#       across the resulting vector of (variable x sector/product)
#       p-values, and reports how many pairs remain significant at
#       q < 0.05 / q < 0.10, next to the raw count at alpha = 0.01.
#
#   (3) panel_granger_placebo_test() -- a circular-shift placebo/surrogate
#       test in the spirit of Theiler et al. (1992, Physica D, 58(1-4),
#       77-94, doi:10.1016/0167-2789(92)90102-S): every explanatory series
#       is independently shifted by a random number of years (with
#       wraparound), which destroys the genuine lead-lag alignment with
#       the dependent series while preserving each series' own
#       autocorrelation. Re-running the UNCHANGED panel_granger_causality()
#       on n_placebo such placebo datasets gives an empirical benchmark for
#       "how many significant pairs would we find by testing-artifact
#       alone"; comparing the OBSERVED count against this benchmark
#       quantifies how much of the reported structure exceeds chance.
#
# All functions use the SAME (Y, X) data structures as the existing
# panel_granger_causality(): Y is a (years x sectors/products) matrix;
# X is a (variables x sectors/products x years) array.
#
# RUNTIME NOTE: panel_granger_placebo_test() re-runs the full panel test
# n_placebo times. Each call is O(n_explanatory x n_sectors x max_lag) OLS
# fits, which is fast (milliseconds each), but with typical dimensions
# (e.g., ~20 variables x ~45 sectors x 5 lags) and n_placebo = 200 this can
# still take several minutes. We recommend starting with n_placebo = 100
# as a pilot, and/or running the placebo test only on the "universal
# predictor" subset flagged by all three methods before scaling up.
#
# This script has not been executed against your actual data (we do not
# have access to it). Please sanity-check on a small example first.
# ---------------------------------------------------------------------

## ---- (1) p-value-preserving version of the panel Granger test --------

panel_granger_causality_pvalues <- function(Y, X, max_lag = 5, impute = TRUE) {

  if (!requireNamespace("imputeTS", quietly = TRUE)) stop("Package 'imputeTS' is required.")
  if (!requireNamespace("lmtest", quietly = TRUE)) stop("Package 'lmtest' is required.")

  if (!is.matrix(Y)) stop("Y must be a matrix")
  if (!is.array(X) || length(dim(X)) != 3) stop("X must be a 3-dimensional array")

  n_years <- nrow(Y); n_sectors <- ncol(Y); n_explanatory <- dim(X)[1]
  if (dim(X)[2] != n_sectors) stop("dim(X)[2] must match ncol(Y)")
  if (dim(X)[3] != n_years)   stop("dim(X)[3] must match nrow(Y)")
  if (max_lag >= n_years)    stop("max_lag must be < n_years")

  if (impute) {
    for (j in 1:n_sectors) {
      if (any(is.na(Y[, j]))) Y[, j] <- imputeTS::na_interpolation(Y[, j], option = "linear")
    }
    for (i in 1:n_explanatory) for (j in 1:n_sectors) {
      if (any(is.na(X[i, j, ]))) X[i, j, ] <- imputeTS::na_interpolation(X[i, j, ], option = "linear")
    }
  }

  granger_test_lag <- function(y, x, lag) {
    tryCatch({
      n <- length(y)
      if (lag >= n - 1) return(NA_real_)
      data_df <- data.frame(y = y[(lag + 1):n])
      for (l in 1:lag) data_df[[paste0("y_lag", l)]] <- y[(lag + 1 - l):(n - l)]
      for (l in 1:lag) data_df[[paste0("x_lag", l)]] <- x[(lag + 1 - l):(n - l)]
      data_df <- na.omit(data_df)
      if (nrow(data_df) < lag + 2) return(NA_real_)
      f0 <- as.formula(paste("y ~", paste(paste0("y_lag", 1:lag), collapse = " + ")))
      f1 <- as.formula(paste("y ~", paste(c(paste0("y_lag", 1:lag), paste0("x_lag", 1:lag)), collapse = " + ")))
      m0 <- lm(f0, data = data_df); m1 <- lm(f1, data = data_df)
      lmtest::waldtest(m0, m1, test = "F")$`Pr(>F)`[2]
    }, error = function(e) NA_real_)
  }

  panel_granger_test <- function(Y_matrix, X_matrix, lag) {
    tryCatch({
      n_years <- nrow(Y_matrix); n_sectors <- ncol(Y_matrix)
      if (lag >= n_years - 1) return(NA_real_)
      panel_data <- data.frame()
      for (j in 1:n_sectors) for (t in (lag + 1):n_years) {
        row_data <- data.frame(sector = j, time = t, y = Y_matrix[t, j])
        for (l in 1:lag) row_data[[paste0("y_lag", l)]] <- Y_matrix[t - l, j]
        for (l in 1:lag) row_data[[paste0("x_lag", l)]] <- X_matrix[t - l, j]
        panel_data <- rbind(panel_data, row_data)
      }
      panel_data <- na.omit(panel_data)
      if (nrow(panel_data) < lag + 2) return(NA_real_)
      f0 <- as.formula(paste("y ~", paste(paste0("y_lag", 1:lag), collapse = " + "), "+ factor(sector)"))
      f1 <- as.formula(paste("y ~", paste(c(paste0("y_lag", 1:lag), paste0("x_lag", 1:lag)), collapse = " + "), "+ factor(sector)"))
      m0 <- lm(f0, data = panel_data); m1 <- lm(f1, data = panel_data)
      lmtest::waldtest(m0, m1, test = "F")$`Pr(>F)`[2]
    }, error = function(e) NA_real_)
  }

  results <- data.frame()
  for (j in 1:n_sectors) {
    for (i in 1:n_explanatory) {
      y_series <- Y[, j]; x_series <- X[i, j, ]
      pvals <- vapply(1:max_lag, function(lag) granger_test_lag(y_series, x_series, lag), numeric(1))
      p_min <- suppressWarnings(min(pvals, na.rm = TRUE)); if (!is.finite(p_min)) p_min <- NA_real_
      n_tested <- sum(is.finite(pvals))
      p_min_bonf <- if (is.finite(p_min) && n_tested > 0) min(1, p_min * n_tested) else NA_real_
      results <- rbind(results, data.frame(
        level = "individual", unit = j, variable = i,
        p_min_raw = p_min, n_lags_tested = n_tested, p_min_lag_bonf = p_min_bonf,
        stringsAsFactors = FALSE))
    }
  }
  for (i in 1:n_explanatory) {
    X_matrix <- t(X[i, , ])
    pvals <- vapply(1:max_lag, function(lag) panel_granger_test(Y, X_matrix, lag), numeric(1))
    p_min <- suppressWarnings(min(pvals, na.rm = TRUE)); if (!is.finite(p_min)) p_min <- NA_real_
    n_tested <- sum(is.finite(pvals))
    p_min_bonf <- if (is.finite(p_min) && n_tested > 0) min(1, p_min * n_tested) else NA_real_
    results <- rbind(results, data.frame(
      level = "pooled", unit = NA_integer_, variable = i,
      p_min_raw = p_min, n_lags_tested = n_tested, p_min_lag_bonf = p_min_bonf,
      stringsAsFactors = FALSE))
  }
  rownames(results) <- NULL
  if (!is.null(dimnames(X)[[1]])) results$variable_name <- dimnames(X)[[1]][results$variable]
  results
}

## ---- (2) Benjamini-Hochberg FDR correction ----------------------------

fdr_correct_panel_results <- function(pvalue_table, p_col = "p_min_lag_bonf",
                                       level_filter = "individual", alpha_raw = 0.01) {
  tab <- if (!is.null(level_filter)) pvalue_table[pvalue_table$level == level_filter, , drop = FALSE] else pvalue_table
  pv <- tab[[p_col]]
  q <- rep(NA_real_, length(pv))
  ok <- is.finite(pv)
  q[ok] <- p.adjust(pv[ok], method = "BH")
  tab$q_bh <- q
  tab$sig_raw_alpha01 <- is.finite(tab[[p_col]]) & tab[[p_col]] < alpha_raw
  tab$sig_fdr_q05     <- is.finite(tab$q_bh) & tab$q_bh < 0.05
  tab$sig_fdr_q10     <- is.finite(tab$q_bh) & tab$q_bh < 0.10

  summary_stats <- c(
    n_tests           = sum(ok),
    n_sig_raw_alpha01  = sum(tab$sig_raw_alpha01,  na.rm = TRUE),
    n_sig_fdr_q05      = sum(tab$sig_fdr_q05,      na.rm = TRUE),
    n_sig_fdr_q10      = sum(tab$sig_fdr_q10,      na.rm = TRUE)
  )
  attr(tab, "summary") <- summary_stats
  cat("FDR correction summary (level =", ifelse(is.null(level_filter), "all", level_filter), "):\n")
  print(summary_stats)
  tab
}

## ---- (3) Circular-shift placebo / surrogate test ----------------------

circular_shift <- function(v, shift) {
  n <- length(v); shift <- shift %% n
  if (shift == 0) return(v)
  c(v[(n - shift + 1):n], v[1:(n - shift)])
}

panel_granger_placebo_test <- function(Y, X, max_lag = 5, p_value = 0.01,
                                        n_placebo = 100, seed = 123, impute = TRUE) {

  if (!exists("panel_granger_causality")) {
    stop("panel_granger_causality() must be loaded first (source('panel_granger_causality.r')).")
  }

  set.seed(seed)
  n_years <- dim(X)[3]

  observed <- panel_granger_causality(Y, X, max_lag = max_lag, p_value = p_value, impute = impute)
  n_sig_observed <- sum(observed > 0)

  placebo_counts <- integer(n_placebo)
  for (b in seq_len(n_placebo)) {
    X_shift <- X
    for (i in 1:dim(X)[1]) {
      for (j in 1:dim(X)[2]) {
        shift_amt <- sample.int(n_years - 1, 1)
        X_shift[i, j, ] <- circular_shift(X[i, j, ], shift_amt)
      }
    }
    placebo_res <- panel_granger_causality(Y, X_shift, max_lag = max_lag, p_value = p_value, impute = impute)
    placebo_counts[b] <- sum(placebo_res > 0)
    if (b %% 10 == 0) cat(sprintf("  placebo iteration %d/%d done\n", b, n_placebo))
  }

  out <- list(
    n_sig_observed = n_sig_observed,
    placebo_counts = placebo_counts,
    placebo_mean   = mean(placebo_counts),
    placebo_sd     = sd(placebo_counts),
    placebo_p95    = as.numeric(quantile(placebo_counts, 0.95)),
    empirical_p    = mean(placebo_counts >= n_sig_observed),
    n_placebo      = n_placebo
  )
  cat(sprintf(
    "\nObserved significant pairs: %d\nPlacebo mean (SD): %.1f (%.1f)\nPlacebo 95th percentile: %.1f\nEmpirical p-value (share of placebos >= observed): %.3f\n",
    out$n_sig_observed, out$placebo_mean, out$placebo_sd, out$placebo_p95, out$empirical_p))
  out
}

# ---------------------------------------------------------------------
# Suggested usage:
#
#   pv_tab   <- panel_granger_causality_pvalues(Y, X, max_lag = 5)
#   fdr_tab  <- fdr_correct_panel_results(pv_tab, level_filter = "individual")
#   # -> report attr(fdr_tab, "summary") as a one-row robustness table
#
#   placebo  <- panel_granger_placebo_test(Y, X, max_lag = 5, p_value = 0.01,
#                                           n_placebo = 100)
#   # -> report placebo$n_sig_observed vs. placebo$placebo_mean /
#   #    placebo$placebo_p95 and placebo$empirical_p
#
# Run this once for the aggregate OECD-ICIO network-level Y/X (as used for
# Figures "causgravsect"/"causgravprod"), and, time permitting, repeat for
# the CEPII-BACI product-level Y/X (Figure "causgravprod"). If the
# empirical_p is small (e.g. < 0.05) and the FDR-corrected count is still a
# large share of the raw count, this is a strong robustness argument; if
# most of the raw "significant" count evaporates after FDR correction and
# is statistically indistinguishable from the placebo distribution, this
# should be reported transparently and the interpretation softened
# accordingly (see The_World_Economy_Rev.md, Points 2, 5 and 6).
# ---------------------------------------------------------------------
