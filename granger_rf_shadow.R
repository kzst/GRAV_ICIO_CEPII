# granger_rf_shadow.R
#
# Shadow-feature ("Boruta-style") significance-corrected version of the
# AGGREGATE-level granger_rf() function (see granger_rf.R).
#
# WHY THIS SCRIPT EXISTS (see The_World_Economy_Rev.md, Point 6):
# The manuscript's Appendix (Sec. "Random Forest-Based Granger-Like
# Causality Analysis", boxed Eq. eq:rf_decision) documents the following
# decision rule for declaring that x_j Granger-precedes y:
#
#   VI_j > 0   AND   max_l VI(j,l)  >  max_{k,l} VI(shadow_{j,l}^(k))
#
# i.e. the variable's best-lag importance must exceed the importance of
# ITS OWN randomly permuted ("shadow") copies. This is exactly the logic
# already implemented for the PANEL RF test in rf_panel_granger_causality.R
# (function create_shadow_features() + the "> max_shadow_importance" check).
#
# However, the existing AGGREGATE/single-series function granger_rf()
# (granger_rf.R) does NOT implement this shadow comparison. Its internal
# "drop" rule is simply:
#
#   drop_flag <- isTRUE(top$importance <= drop_threshold)   # drop_threshold = 0 by default
#
# i.e. ANY strictly positive importance -- however small -- is treated as
# evidence of Granger precedence. Because permutation/impurity importances
# for continuous predictors are rarely exactly zero or negative, this bar
# is very easy to cross by chance, especially with only 26-28 annual
# observations. We believe this inconsistency (not a difference in
# modelling philosophy) is a leading explanation for why the RF approach
# "consistently identif[ies] more causal relationships" than the
# traditional and Bayesian tests (as the manuscript itself reports for
# Figure "rfnetcaus" and the "_rf" panels of Figures "causgravsect" /
# "causgravprod").
#
# This script does NOT modify granger_rf.R (so your existing, already
# validated pipeline stays intact); it adds a new function,
# granger_rf_shadow(), that you can run ALONGSIDE granger_rf() and compare
# against it. We recommend reporting both counts (pre- and post-shadow
# -correction) as a one-line robustness comparison in the Results text.
#
# IMPORTANT: this script has not been executed against your actual data
# (we do not have access to it). Please sanity-check it on one small,
# already-known example (e.g., the aggregate OECD-ICIO series used for
# Figure "netcaus") before relying on it for the full pipeline.
#
# ---------------------------------------------------------------------

granger_rf_shadow <- function(y, X,
                               max_lag = 5,
                               ntree = 500,
                               mtry = NULL,
                               min_node_size = 5,
                               importance_mode = c("permutation", "impurity"),
                               n_shadow = 20,
                               seed = 123) {

  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' is required. Install with install.packages('ranger').")
  }
  importance_mode <- match.arg(importance_mode)

  # ---- basic checks (mirrors granger_rf.R) ----
  y <- as.numeric(y)
  if (any(!is.finite(y))) stop("y must be numeric and finite.")
  if (is.matrix(X)) X <- as.data.frame(X)
  if (!is.data.frame(X)) stop("X must be a data.frame or matrix.")
  if (length(y) != nrow(X)) stop("y and X must have the same length/rows.")
  if (max_lag < 1) stop("max_lag must be >= 1.")
  if (n_shadow < 1) stop("n_shadow must be >= 1.")

  n <- length(y)
  orig_names <- colnames(X)
  if (is.null(orig_names)) orig_names <- paste0("X", seq_len(ncol(X)))
  safe_names <- make.names(orig_names, unique = TRUE)
  colnames(X) <- safe_names

  for (j in seq_along(safe_names)) {
    if (!is.numeric(X[[j]])) suppressWarnings(X[[j]] <- as.numeric(X[[j]]))
    if (any(!is.finite(X[[j]]))) {
      stop(sprintf("Column '%s' in X contains non-finite values after coercion.", orig_names[j]))
    }
  }

  build_lags <- function(vec, L) {
    nn <- length(vec)
    out <- matrix(NA_real_, nrow = nn, ncol = L)
    for (l in 1:L) out[(l + 1):nn, l] <- vec[1:(nn - l)]
    out
  }

  # ---- real (non-shadow) lagged predictors: y's own lags + every X lag ----
  df <- data.frame(y = y)
  y_lags <- build_lags(y, max_lag)
  colnames(y_lags) <- sprintf("y_L%d", 1:max_lag)
  df <- cbind(df, y_lags)

  for (j in seq_along(safe_names)) {
    lags_mat <- build_lags(X[[j]], max_lag)
    colnames(lags_mat) <- sprintf("%s_L%d", safe_names[j], 1:max_lag)
    df <- cbind(df, lags_mat)
  }

  df <- df[(max_lag + 1):n, , drop = FALSE]
  row.names(df) <- NULL
  complete_idx <- stats::complete.cases(df)
  df <- df[complete_idx, , drop = FALSE]
  if (nrow(df) < (max_lag + 5)) {
    stop("Too few complete observations after lagging/trimming for a reliable RF fit.")
  }

  set.seed(seed)

  # ---- shadow copies: for every X variable j, n_shadow permuted replicas
  #      of EACH of its L lag columns (matches the manuscript's per-variable
  #      shadow definition and rf_panel_granger_causality.R) ----
  shadow_df <- data.frame(row.names = seq_len(nrow(df)))
  for (j in seq_along(safe_names)) {
    for (l in 1:max_lag) {
      col_name <- sprintf("%s_L%d", safe_names[j], l)
      for (s in seq_len(n_shadow)) {
        shadow_name <- sprintf("%s_shadow%d", col_name, s)
        shadow_df[[shadow_name]] <- sample(df[[col_name]])
      }
    }
  }

  full_data <- cbind(df, shadow_df)
  p_full <- ncol(full_data) - 1L
  mtry_full <- if (is.null(mtry)) max(1L, floor(sqrt(p_full))) else max(1L, mtry)

  rf_final <- ranger::ranger(
    formula = y ~ .,
    data = full_data,
    num.trees = ntree,
    mtry = mtry_full,
    min.node.size = min_node_size,
    importance = importance_mode,
    seed = seed
  )
  varimp <- rf_final$variable.importance

  get_imp <- function(nm) {
    v <- varimp[[nm]]
    if (is.null(v) || !is.finite(v)) 0 else v
  }

  decision_rows <- vector("list", length(safe_names))
  lag_rows      <- vector("list", length(safe_names))

  for (j in seq_along(safe_names)) {
    lag_names_j <- sprintf("%s_L%d", safe_names[j], 1:max_lag)
    vi_j <- vapply(lag_names_j, get_imp, numeric(1))

    shadow_pattern <- sprintf("^%s_L[0-9]+_shadow[0-9]+$", safe_names[j])
    shadow_names_j <- grep(shadow_pattern, names(varimp), value = TRUE)
    shadow_vi_j <- vapply(shadow_names_j, get_imp, numeric(1))
    max_shadow_j <- if (length(shadow_vi_j)) max(shadow_vi_j, na.rm = TRUE) else 0

    per_lag_shadow_max <- vapply(1:max_lag, function(l) {
      pat_l <- sprintf("^%s_L%d_shadow[0-9]+$", safe_names[j], l)
      nm_l <- grep(pat_l, names(varimp), value = TRUE)
      if (length(nm_l)) max(vapply(nm_l, get_imp, numeric(1))) else NA_real_
    }, numeric(1))

    best_lag  <- which.max(vi_j)
    best_imp  <- vi_j[best_lag]
    vi_total  <- sum(vi_j)

    granger_precedes <- isTRUE(vi_total > 0) && isTRUE(best_imp > max_shadow_j)

    decision_rows[[j]] <- data.frame(
      variable               = orig_names[j],
      VI_total                = vi_total,
      best_lag                = best_lag,
      best_lag_importance     = best_imp,
      max_shadow_importance   = max_shadow_j,
      expected_lag            = if (granger_precedes) best_lag else 0L,
      granger_precedes_shadow = granger_precedes,
      granger_precedes_raw    = isTRUE(vi_total > 0) && isTRUE(best_imp > 0),  # old granger_rf.R rule, for comparison
      stringsAsFactors = FALSE
    )

    lag_rows[[j]] <- data.frame(
      variable           = orig_names[j],
      lag                = 1:max_lag,
      importance         = vi_j,
      max_shadow_this_lag = per_lag_shadow_max,
      stringsAsFactors = FALSE
    )
  }

  decision <- do.call(rbind, decision_rows)
  row.names(decision) <- NULL
  decision <- decision[order(-decision$granger_precedes_shadow, -decision$best_lag_importance), , drop = FALSE]

  lag_importance <- do.call(rbind, lag_rows)
  row.names(lag_importance) <- NULL

  vi_pos <- pmax(decision$VI_total, 0)
  variable_importance <- data.frame(
    variable = decision$variable,
    importance = if (sum(vi_pos) > 0) vi_pos / sum(vi_pos) else 0,
    granger_precedes_shadow = decision$granger_precedes_shadow,
    granger_precedes_raw    = decision$granger_precedes_raw,
    stringsAsFactors = FALSE
  )
  variable_importance <- variable_importance[order(-variable_importance$importance), , drop = FALSE]

  list(
    decision = decision,
    variable_importance = variable_importance,
    lag_importance = lag_importance,
    n_flagged_shadow = sum(decision$granger_precedes_shadow),
    n_flagged_raw    = sum(decision$granger_precedes_raw),
    details = list(
      n_effective = nrow(df),
      max_lag = max_lag,
      n_shadow = n_shadow,
      importance_mode = importance_mode,
      ntree = ntree,
      mtry = mtry_full,
      seed = seed
    )
  )
}

# ---------------------------------------------------------------------
# Suggested one-line robustness comparison for the Results text, e.g.:
#
#   res <- granger_rf_shadow(y = <your y series>, X = <your X data.frame>,
#                             max_lag = 5, n_shadow = 20)
#   cat("Variables flagged before shadow correction:", res$n_flagged_raw, "\n")
#   cat("Variables flagged after  shadow correction:", res$n_flagged_shadow, "\n")
#   print(res$decision)
#
# Re-run once per gravity coefficient series (distance, origin-GDP,
# destination-GDP; aggregate + each sector/product) exactly where
# granger_rf() is currently called, and report the "before/after" counts
# alongside your existing Figure "rfnetcaus" / "_rf" panels.
# ---------------------------------------------------------------------
