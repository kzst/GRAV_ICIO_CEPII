granger_rf <- function(y, X,
                       max_lag = 5,
                       ntree = 500,
                       mtry = NULL,
                       min_node_size = 5,
                       importance_mode = c("permutation", "impurity"),
                       oos_folds = 5,
                       include_y_lags = TRUE,
                       seed = 123,
                       return_models = FALSE,
                       drop_threshold = 0) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' is required. Install with install.packages('ranger').")
  }
  importance_mode <- match.arg(importance_mode)
  
  # Basic checks and preparation
  y <- as.numeric(y)
  if (any(!is.finite(y))) stop("y must be numeric and finite.")
  if (is.matrix(X)) X <- as.data.frame(X)
  if (!is.data.frame(X)) stop("X must be a data.frame or matrix.")
  if (length(y) != nrow(X)) stop("y and X must have the same length/rows.")
  if (max_lag < 1) stop("max_lag must be >= 1.")
  if (oos_folds < 2) stop("oos_folds must be >= 2.")
  
  n <- length(y)
  orig_names <- colnames(X)
  if (is.null(orig_names)) orig_names <- paste0("X", seq_len(ncol(X)))
  safe_names <- make.names(orig_names, unique = TRUE)
  colnames(X) <- safe_names
  
  # Ensure numeric columns
  for (j in seq_along(safe_names)) {
    if (!is.numeric(X[[j]])) suppressWarnings(X[[j]] <- as.numeric(X[[j]]))
    if (any(!is.finite(X[[j]]))) {
      stop(sprintf("Column '%s' in X contains non-finite values after coercion.", orig_names[j]))
    }
  }
  
  # Helper to build lag matrix
  build_lags <- function(vec, L) {
    n <- length(vec)
    out <- matrix(NA_real_, nrow = n, ncol = L)
    for (l in 1:L) {
      out[(l + 1):n, l] <- vec[1:(n - l)]
    }
    out
  }
  
  # Construct lagged dataset
  df <- data.frame(y = y)
  if (include_y_lags) {
    y_lags <- build_lags(y, max_lag)
    colnames(y_lags) <- sprintf("y_L%d", 1:max_lag)
    df <- cbind(df, y_lags)
  }
  for (j in seq_along(safe_names)) {
    lags_mat <- build_lags(X[[j]], max_lag)
    colnames(lags_mat) <- sprintf("%s_L%d", safe_names[j], 1:max_lag)
    df <- cbind(df, lags_mat)
  }
  
  # Drop initial rows and any remaining NA rows
  df <- df[(max_lag + 1):n, , drop = FALSE]
  row.names(df) <- NULL
  complete_idx <- stats::complete.cases(df)
  removed <- sum(!complete_idx)
  if (removed > 0) {
    df <- df[complete_idx, , drop = FALSE]
  }
  
  # Feature sets
  all_feature_names <- setdiff(names(df), "y")
  y_feature_names <- if (include_y_lags) grep("^y_L\\d+$", names(df), value = TRUE) else character(0)
  X_feature_names <- setdiff(all_feature_names, y_feature_names)
  
  p_full <- length(all_feature_names)
  p_base <- length(y_feature_names)
  if (is.null(mtry)) {
    mtry_full <- max(1L, floor(sqrt(max(1, p_full))))
    mtry_base <- max(1L, floor(sqrt(max(1, p_base))))
  } else {
    mtry_full <- max(1L, mtry)
    mtry_base <- max(1L, min(mtry, max(1, p_base)))
  }
  
  set.seed(seed)
  
  # Blocked (time-aware) cross-validation
  n_obs <- nrow(df)
  fold_id <- as.integer(cut(seq_len(n_obs), breaks = oos_folds, labels = FALSE))
  test_folds <- sort(unique(fold_id))[-1]  # folds 2..K
  
  sse_base <- 0
  sse_full <- 0
  n_test_total <- 0
  y_test_all <- numeric(0)
  
  for (f in test_folds) {
    train_idx <- which(fold_id < f)
    test_idx  <- which(fold_id == f)
    
    # Baseline: AR on y's lags (or mean if no y lags)
    train_df_base <- data.frame(y = df$y[train_idx])
    if (p_base > 0) train_df_base <- cbind(train_df_base, df[train_idx, y_feature_names, drop = FALSE])
    test_df_base  <- if (p_base > 0) df[test_idx, y_feature_names, drop = FALSE] else NULL
    
    if (p_base > 0) {
      rf_base <- ranger::ranger(
        formula = y ~ .,
        data = train_df_base,
        num.trees = ntree,
        mtry = mtry_base,
        min.node.size = min_node_size,
        importance = "none",
        seed = seed
      )
      pred_base <- predict(rf_base, data = test_df_base)$predictions
    } else {
      pred_base <- rep(mean(train_df_base$y), length(test_idx))
    }
    
    # Full model: y lags + X lags
    train_df_full <- data.frame(y = df$y[train_idx], df[train_idx, all_feature_names, drop = FALSE])
    test_df_full  <- df[test_idx, all_feature_names, drop = FALSE]
    rf_full_cv <- ranger::ranger(
      formula = y ~ .,
      data = train_df_full,
      num.trees = ntree,
      mtry = mtry_full,
      min.node.size = min_node_size,
      importance = "none",
      seed = seed
    )
    pred_full <- predict(rf_full_cv, data = test_df_full)$predictions
    
    # Accumulate errors
    y_test <- df$y[test_idx]
    sse_base <- sse_base + sum((y_test - pred_base)^2)
    sse_full <- sse_full + sum((y_test - pred_full)^2)
    n_test_total <- n_test_total + length(test_idx)
    y_test_all <- c(y_test_all, y_test)
  }
  
  mse_base <- sse_base / n_test_total
  mse_full <- sse_full / n_test_total
  rel_mse_reduction <- 1 - (mse_full / mse_base)
  var_ref <- stats::var(y_test_all)
  oos_r2_against_mean <- if (is.finite(var_ref) && var_ref > 0) 1 - (mse_full / var_ref) else NA_real_
  
  # Final full model on all data for importances
  final_data <- data.frame(y = df$y, df[, all_feature_names, drop = FALSE])
  rf_final <- ranger::ranger(
    formula = y ~ .,
    data = final_data,
    num.trees = ntree,
    mtry = mtry_full,
    min.node.size = min_node_size,
    importance = importance_mode,
    seed = seed
  )
  varimp <- rf_final$variable.importance
  feat_names <- names(varimp)
  variables <- sub("_L\\d+$", "", feat_names)
  lags <- as.integer(sub(".*_L", "", feat_names))
  
  import_df <- data.frame(
    feature = feat_names,
    variable = variables,
    lag = lags,
    importance = as.numeric(varimp),
    stringsAsFactors = FALSE
  )
  
  # Map safe variable names back to original (keep y as 'y')
  name_map <- setNames(orig_names, safe_names)  # safe -> original
  import_df$variable_orig <- import_df$variable
  idx_x <- import_df$variable != "y"
  import_df$variable_orig[idx_x] <- unname(name_map[import_df$variable[idx_x]])
  
  # X-only importance summaries
  x_imp_df <- import_df[import_df$variable != "y", c("variable_orig", "lag", "importance"), drop = FALSE]
  names(x_imp_df)[1] <- "variable"
  
  variable_imp <- aggregate(importance ~ variable, data = x_imp_df, sum)
  if (sum(variable_imp$importance) > 0) {
    variable_imp$importance <- variable_imp$importance / sum(variable_imp$importance)
  }
  variable_imp <- variable_imp[order(-variable_imp$importance), , drop = FALSE]
  variable_imp$rank <- seq_len(nrow(variable_imp))
  row.names(variable_imp) <- NULL
  
  lag_imp <- x_imp_df
  lag_imp$rank_within_variable <- NA_real_
  for (v in unique(lag_imp$variable)) {
    idx <- which(lag_imp$variable == v)
    lag_imp$rank_within_variable[idx] <- rank(-lag_imp$importance[idx], ties.method = "average")
  }
  lag_imp <- lag_imp[order(lag_imp$variable, -lag_imp$importance, lag_imp$lag), , drop = FALSE]
  row.names(lag_imp) <- NULL
  
  lag_overall <- aggregate(importance ~ lag, data = x_imp_df, sum)
  if (sum(lag_overall$importance) > 0) {
    lag_overall$importance <- lag_overall$importance / sum(lag_overall$importance)
  }
  lag_overall <- lag_overall[order(lag_overall$lag), , drop = FALSE]
  row.names(lag_overall) <- NULL
  
  y_imp_df <- import_df[import_df$variable == "y", c("lag", "importance"), drop = FALSE]
  y_imp_df <- y_imp_df[order(y_imp_df$lag), , drop = FALSE]
  row.names(y_imp_df) <- NULL
  
  # expected_lags: most likely lag per variable (argmax importance; tie -> smallest lag)
  # drop flag: TRUE if best lag importance <= drop_threshold; then expected_lag = 0
  expected_list <- lapply(split(x_imp_df, x_imp_df$variable), function(d) {
    d <- d[order(-d$importance, d$lag), , drop = FALSE]
    top <- d[1L, , drop = FALSE]
    total <- sum(d$importance)
    lag_share <- if (total > 0) top$importance / total else NA_real_
    drop_flag <- isTRUE(top$importance <= drop_threshold)
    data.frame(
      variable = top$variable,
      expected_lag = if (drop_flag) 0L else as.integer(top$lag),
      lag_importance = as.numeric(top$importance),
      lag_share = as.numeric(lag_share),
      drop = drop_flag,
      stringsAsFactors = FALSE
    )
  })
  expected_df <- do.call(rbind, expected_list)
  row.names(expected_df) <- NULL
  expected_df <- expected_df[order(match(expected_df$variable, variable_imp$variable)), , drop = FALSE]
  
  out <- list(
    variable_importance = variable_imp,
    lag_importance = lag_imp,
    lag_importance_overall = lag_overall,
    y_lag_importance = y_imp_df,
    expected_lags = expected_df,
    oos_metrics = list(
      mse_baseline = mse_base,
      mse_full = mse_full,
      relative_mse_reduction = rel_mse_reduction,
      oos_r2_against_mean = as.numeric(oos_r2_against_mean)
    ),
    details = list(
      n = n,
      n_effective = nrow(df),
      removed_rows_due_to_na = removed,
      max_lag = max_lag,
      importance_mode = importance_mode,
      oos_folds = oos_folds,
      ntree = ntree,
      mtry_full = mtry_full,
      mtry_base = mtry_base,
      min_node_size = min_node_size,
      seed = seed,
      drop_threshold = drop_threshold
    )
  )
  if (return_models) out$models <- list(final_full_model = rf_final)
  class(out) <- c("granger_rf_result", class(out))
  out
}