# panel_fdr_reporting_helpers.R
#
# Convenience helpers for turning the LONG-FORMAT output of
# panel_granger_causality_pvalues() / fdr_correct_panel_results() /
# panel_granger_placebo_test() (see panel_granger_robustness_fdr_placebo.R)
# into two paper-ready things:
#
#   (1) build_sig_matrix()      -- a TRUE/FALSE matrix with the SAME shape
#       as TRAD_DIST/TRAD_GDPi/TRAD_GDPj (n_explanatory x (n_sectors+1),
#       last column = "Pooled"), so it can be overlaid cell-by-cell on the
#       EXISTING Heatmap(TRAD_*, ...) figures via cell_fun -- no new figure
#       needed, no 800+ row table needed.
#
#   (2) summarize_panel_fdr()   -- collapses fdr_correct_panel_results()
#       output for several Y-coefficients (Distance/Origin/Destination)
#       into ONE small data.frame (one row per coefficient x level),
#       ready for tokable().
#
#   (3) summarize_panel_placebo() -- same idea for panel_granger_placebo_test()
#       output: one row per coefficient.
#
# Does not modify panel_granger_causality.r, panel_granger_robustness_fdr_placebo.R,
# or the existing "causgravsect"/"causgravprod" chunks -- purely additive.
#
# This script has not been executed against your actual data (we do not
# have access to it). Please sanity-check on a small example first.
# ---------------------------------------------------------------------

## ---- (1) FDR-significance marker matrix, same shape as TRAD_*/BAYS_*/RFBS_* ----
## mat_template   : any of TRAD_DIST/TRAD_GDPi/TRAD_GDPj (used only for dim/dimnames)
## fdr_individual : fdr_correct_panel_results(pv_tab, level_filter = "individual")
## fdr_pooled     : fdr_correct_panel_results(pv_tab, level_filter = "pooled")
## sig_col        : "sig_fdr_q05" (default) or "sig_fdr_q10" or "sig_raw_alpha01"
## Returns a logical matrix, same dim/dimnames as mat_template.

build_sig_matrix <- function(mat_template, fdr_individual, fdr_pooled,
                              sig_col = "sig_fdr_q05") {
  n_var    <- nrow(mat_template)
  n_sector <- ncol(mat_template) - 1L   # last column of panel_granger_causality() output is "Pooled"
  m <- matrix(FALSE, nrow = n_var, ncol = n_sector + 1L,
              dimnames = dimnames(mat_template))

  for (k in seq_len(nrow(fdr_individual))) {
    i <- fdr_individual$variable[k]; j <- fdr_individual$unit[k]
    if (is.finite(i) && is.finite(j) && i <= n_var && j <= n_sector) {
      m[i, j] <- isTRUE(fdr_individual[[sig_col]][k])
    }
  }
  for (k in seq_len(nrow(fdr_pooled))) {
    i <- fdr_pooled$variable[k]
    if (is.finite(i) && i <= n_var) {
      m[i, n_sector + 1L] <- isTRUE(fdr_pooled[[sig_col]][k])
    }
  }
  m
}

## ---- (2) Compact FDR summary table across several Y-coefficients -----
## fdr_list        : named list, e.g. list(Distance = FDR_DIST, `Origin GDP` = FDR_GDPi, ...)
##                   each element = fdr_correct_panel_results(..., level_filter = "individual")
## fdr_pooled_list : same names, each element = fdr_correct_panel_results(..., level_filter = "pooled")
## Returns one data.frame: one row per (coefficient x level) -- a HANDFUL of rows, not 800+.

summarize_panel_fdr <- function(fdr_list, fdr_pooled_list) {
  rows <- list()
  for (nm in names(fdr_list)) {
    s <- attr(fdr_list[[nm]], "summary")
    rows[[length(rows) + 1]] <- data.frame(
      Coefficient = nm, Level = "Individual (sector-level)",
      N_tested    = unname(s["n_tests"]),
      N_sig_raw   = unname(s["n_sig_raw_alpha01"]),
      N_sig_FDR05 = unname(s["n_sig_fdr_q05"]),
      N_sig_FDR10 = unname(s["n_sig_fdr_q10"]),
      stringsAsFactors = FALSE)
  }
  for (nm in names(fdr_pooled_list)) {
    s <- attr(fdr_pooled_list[[nm]], "summary")
    rows[[length(rows) + 1]] <- data.frame(
      Coefficient = nm, Level = "Pooled (panel, sector FE)",
      N_tested    = unname(s["n_tests"]),
      N_sig_raw   = unname(s["n_sig_raw_alpha01"]),
      N_sig_FDR05 = unname(s["n_sig_fdr_q05"]),
      N_sig_FDR10 = unname(s["n_sig_fdr_q10"]),
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  out$Pct_retained_q05 <- ifelse(out$N_sig_raw > 0,
                                  round(100 * out$N_sig_FDR05 / out$N_sig_raw, 1), NA_real_)
  out
}

## ---- (3) Compact placebo summary table --------------------------------
## placebo_list : named list, e.g. list(Distance = PLACEBO_DIST, `Origin GDP` = PLACEBO_GDPi, ...)
##                each element = panel_granger_placebo_test() output

summarize_panel_placebo <- function(placebo_list) {
  rows <- lapply(names(placebo_list), function(nm) {
    p <- placebo_list[[nm]]
    data.frame(
      Coefficient    = nm,
      N_sig_observed = p$n_sig_observed,
      Placebo_mean   = round(p$placebo_mean, 1),
      Placebo_SD     = round(p$placebo_sd, 1),
      Placebo_p95    = round(p$placebo_p95, 1),
      Empirical_p    = round(p$empirical_p, 3),
      Robust_at_5pct = p$empirical_p < 0.05,
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------
# Suggested usage (after sourcing panel_granger_causality.r and
# panel_granger_robustness_fdr_placebo.R, and computing PV_*/FDR_*/
# FDR_*_POOLED/PLACEBO_* as shown in chat):
#
#   SIG_DIST <- build_sig_matrix(TRAD_DIST, FDR_DIST, FDR_DIST_POOLED)
#   # -> use SIG_DIST[i, j] inside the cell_fun of Heatmap(TRAD_DIST, ...)
#
#   fdr_summary <- summarize_panel_fdr(
#     list(Distance = FDR_DIST, `Origin GDP` = FDR_GDPi, `Destination GDP` = FDR_GDPj),
#     list(Distance = FDR_DIST_POOLED, `Origin GDP` = FDR_GDPi_POOLED, `Destination GDP` = FDR_GDPj_POOLED))
#
#   placebo_summary <- summarize_panel_placebo(
#     list(Distance = PLACEBO_DIST, `Origin GDP` = PLACEBO_GDPi, `Destination GDP` = PLACEBO_GDPj))
# ---------------------------------------------------------------------
