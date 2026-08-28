granger_vars_simple <- function(
    y, X,
    max_lag = 5,
    alpha = 0.01,
    BF10 = 5,
    lag_select = c("AIC","BIC","HQ","FPE","fixed"),
    include_y_lags = TRUE,
    impute = TRUE,
    impute_method = c("kalman","interpolation","locf"),
    seed = 123
) {
  lag_select   <- match.arg(lag_select)
  impute_method <- match.arg(impute_method)
  
  if (!requireNamespace("vars", quietly = TRUE)) {
    stop("Package 'vars' is required. Install with install.packages('vars').")
  }
  if (impute && !requireNamespace("imputeTS", quietly = TRUE)) {
    stop("Package 'imputeTS' is required. Install with install.packages('imputeTS').")
  }
  if (impute && !requireNamespace("biclust", quietly = TRUE)) {
    stop("Package 'biclust' is required. Install with install.packages('biclust').")
  }
  set.seed(seed)
  
  # -- Imputation helper --
  impute_series <- function(v) {
    v <- as.numeric(v); v[!is.finite(v)] <- NA_real_
    if (!impute) return(v)
    if (all(is.na(v))) return(rep(0, length(v)))
    out <- try({
      if (impute_method == "kalman") {
        imputeTS::na_kalman(v, model = "auto.arima", smooth = TRUE)
      } else if (impute_method == "interpolation") {
        imputeTS::na_interpolation(v, option = "linear")
      } else {
        imputeTS::na_locf(v, option = "locf", na_remaining = "nocb")
      }
    }, silent = TRUE)
    if (inherits(out, "try-error")) {
      out <- try(imputeTS::na_interpolation(v, option = "linear"), silent = TRUE)
      if (inherits(out, "try-error")) out <- imputeTS::na_locf(v, option = "locf", na_remaining = "nocb")
    }
    out
  }
  
  # -- Coerce and impute inputs --
  y <- impute_series(y)
  if (is.matrix(X)) X <- as.data.frame(X)
  if (!is.data.frame(X)) stop("X must be a data.frame or matrix.")
  if (length(y) != nrow(X)) stop("y and X must have the same number of rows.")
  
  orig_names <- colnames(X)
  if (is.null(orig_names)) orig_names <- paste0("X", seq_len(ncol(X)))
  safe_names <- make.names(orig_names, unique = TRUE)
  colnames(X) <- safe_names
  
  for (j in seq_along(safe_names)) {
    if (!is.numeric(X[[j]])) suppressWarnings(X[[j]] <- as.numeric(X[[j]]))
    X[[j]] <- impute_series(X[[j]])
  }
  n <- length(y)
  
  joint_df <- data.frame(y = y, X, check.names = FALSE)
  joint_df <- joint_df[stats::complete.cases(joint_df), , drop = FALSE]
  n_eff <- nrow(joint_df)
  if (n_eff < 10) stop("Not enough effective observations after trimming/imputation.")
  
  # Global lag order selection
  max_lag <- max(1, as.integer(max_lag))
  if (lag_select != "fixed") {
    p <- NA_integer_
    sel <- try(vars::VARselect(joint_df, lag.max = max_lag, type = "const"), silent = TRUE)
    if (!inherits(sel, "try-error")) {
      crit_map <- c(AIC="AIC(n)", BIC="SC(n)", HQ="HQ(n)", FPE="FPE(n)")
      key <- crit_map[[lag_select]]
      if (!is.null(sel$selection[[key]])) p <- as.integer(sel$selection[[key]])
    }
    if (!is.finite(p) || p < 1) p <- 1L
  } else {
    p <- max_lag
  }
  # Safe cap for short series (bivariate; conservative)
  p <- max(1L, min(p, max_lag, floor((n_eff - 5)/4)))
  if (!is.finite(p) || p < 1) p <- 1L
  
  # Build y-lag importance via AR(p) OLS
  build_lags <- function(vec, L) {
    nn <- length(vec)
    out <- matrix(NA_real_, nrow = nn, ncol = L)
    for (l in 1:L) out[(l+1):nn, l] <- vec[1:(nn - l)]
    colnames(out) <- paste0("L", 1:L)
    out
  }
  
  if (p > 0) {
    ylags <- build_lags(joint_df$y, p); colnames(ylags) <- paste0("y_L", 1:p)
    ar_df <- data.frame(y = joint_df$y, ylags, check.names = FALSE)
    ar_df <- ar_df[(p + 1):nrow(ar_df), , drop = FALSE]
    if (nrow(ar_df) >= (p + 2)) {
      ar_fit <- stats::lm(stats::as.formula(paste("y ~", paste(colnames(ylags), collapse=" + "))), data = ar_df)
      ar_coef <- summary(ar_fit)$coefficients
      y_imp <- data.frame(
        lag = 1:p,
        importance = sapply(paste0("y_L", 1:p), function(nm) {
          if (nm %in% rownames(ar_coef)) {
            val <- ar_coef[nm, "t value"]; if (is.finite(val)) abs(val) else 0
          } else 0
        }),
        stringsAsFactors = FALSE
      )
      tot_y <- sum(y_imp$importance, na.rm = TRUE)
      y_imp$importance <- if (tot_y > 0) y_imp$importance / tot_y else 0
    } else {
      y_imp <- data.frame(lag = 1:p, importance = 0)
    }
  } else {
    y_imp <- data.frame(lag = integer(0), importance = numeric(0))
  }
  
  # Containers
  lag_imp_rows <- list()
  results <- data.frame(
    variable = orig_names,
    p_value = NA_real_,
    expected_lag = 0L,
    lag_importance = 0,
    lag_share = 0,
    drop = TRUE,
    stringsAsFactors = FALSE
  )
  
  # Helper to create a non-empty lag-importance DF (length = p_used)
  make_zero_lag_imp <- function(var_name, p_used) {
    data.frame(
      variable = rep(var_name, p_used),
      lag = seq_len(p_used),
      importance = rep(0, p_used),
      rank_within_variable = rep(NA_real_, p_used),
      stringsAsFactors = FALSE
    )
  }
  
  # Per-variable (X_j -> y) bivariate VAR
  for (j in seq_along(safe_names)) {
    df_bi <- joint_df[, c("y", safe_names[j]), drop = FALSE]
    n_bi <- nrow(df_bi)
    p_pair <- max(1L, min(p, floor((n_bi - 5)/4)))
    if (!is.finite(p_pair) || p_pair < 1) p_pair <- 1L
    
    #fit <- NULL; p_try <- p_pair
    #while (p_try >= 1L && is.null(fit)) {
    #  fit_try <- try(vars::VAR(df_bi, p = p_try, type = "const"), silent = TRUE)
    #  if (!inherits(fit_try, "try-error")) fit <- fit_try else p_try <- p_try - 1L
    #}
    
    
    fit <- NULL; p_try <- 1
    BF<-0
    pval<-1
    LL<-TRUE
    while (LL) {
      fit_try <- vars::VAR(df_bi, p = p_try, type = "const")
      gc <-vars::causality(fit_try, cause = safe_names[j])
      pval <- as.numeric(gc$Granger$p.value)
      if (is.na(pval)) pval<-1
      BF<-bayesian_granger_causality(df_bi[,safe_names[j]],df_bi[,"y"],lag=p_try)$bayes_factor
      fit <- fit_try
      if (p_try>max_lag) LL<-FALSE
      if ((pval<alpha)&&(BF>BF10)) LL<-FALSE  
      if (LL) p_try<-p_try+1
    }
    
    if (p_try>max_lag) {
      lag_imp_rows[[length(lag_imp_rows) + 1]] <- make_zero_lag_imp(orig_names[j], 1L)
      results$p_value[j] <- NA_real_
      results$expected_lag[j] <- 0L
      results$lag_importance[j] <- 0
      results$lag_share[j] <- 0
      results$drop[j] <- TRUE
      next
    }else{
      p_try<-fit$p
    }
    p_used <- fit$p
    
    # Granger test x_j -> y
    #pval <- NA_real_
    #gc <- try(vars::causality(fit, cause = safe_names[j]), silent = TRUE)
    #if (!inherits(gc, "try-error") && !is.null(gc$Granger$p.value)) {
    #  pval <- as.numeric(gc$Granger$p.value)
    #  if (!is.finite(pval)) pval <- NA_real_
    #}
    
    # y-equation t-stats for x_j lags
    lag_imp <- make_zero_lag_imp(orig_names[j], p_used)
    if (!is.null(fit$varresult$y)) {
      summ <- summary(fit$varresult$y)
      coefs <- summ$coefficients
      imp_vals <- BF #numeric(p_used)
      for (L in seq_len(p_used)) {
        nm <- paste0(safe_names[j], ".l", L)
        imp_vals[L] <- if (nm %in% rownames(coefs)) {
          val <- bayesian_granger_causality(df_bi[,"y"],df_bi[,safe_names[j]],lag=L)$bayes_factor; if (is.finite(val)) abs(val) else 0
        } else 0
      }
      lag_imp$importance <- imp_vals
      lag_imp$rank_within_variable <- if (any(imp_vals > 0)) rank(-imp_vals, ties.method = "average") else NA_real_
    }
    lag_imp_rows[[length(lag_imp_rows) + 1]] <- lag_imp
    
    best_lag <- if (any(lag_imp$importance > 0)) p_try else 1L
    best_imp <- if (any(lag_imp$importance > 0)) BF else 0
    tot_imp  <- sum(lag_imp$importance, na.rm = TRUE)
    drop_flag <- !(is.finite(pval) && (pval <= alpha))
    
    results$p_value[j]       <- pval
    results$expected_lag[j]  <- if (drop_flag) 0L else as.integer(best_lag)
    results$lag_importance[j]<- if (drop_flag) 0 else best_imp
    results$lag_share[j]     <- if (!drop_flag && tot_imp > 0) best_imp / tot_imp else 0
    results$drop[j]          <- drop_flag
  }
  
  # Relative importance ONLY over significant variables; sum = 1
  vi <- data.frame(variable = results$variable, importance = 0, stringsAsFactors = FALSE)
  signif_mask <- is.finite(results$p_value) & (results$p_value <= alpha)
  scores <- ifelse(signif_mask & (results$p_value > 0),
                   -log10(pmax(results$p_value, .Machine$double.xmin)),
                   0)
  sum_scores <- sum(scores[signif_mask], na.rm = TRUE)
  vi$importance <- if (sum_scores > 0) {
    ifelse(signif_mask, scores / sum_scores, 0)
  } else {
    0
  }
  vi <- vi[order(-vi$importance, vi$variable), , drop = FALSE]
  vi$rank <- if (nrow(vi) > 0) seq_len(nrow(vi)) else integer(0)
  rownames(vi) <- NULL
  
  # Lag importance overall (sum over variables), normalized
  lag_imp_df <- if (length(lag_imp_rows)) do.call(rbind, lag_imp_rows) else
    data.frame(variable = character(0), lag = integer(0), importance = numeric(0), rank_within_variable = numeric(0))
  lag_overall <- if (nrow(lag_imp_df) > 0) aggregate(importance ~ lag, data = lag_imp_df, sum, na.rm = TRUE) else
    data.frame(lag = integer(0), importance = numeric(0))
  pos_sum_lag <- sum(lag_overall$importance[lag_overall$importance > 0], na.rm = TRUE)
  lag_overall$importance <- if (pos_sum_lag > 0) ifelse(lag_overall$importance > 0,
                                                        lag_overall$importance / pos_sum_lag, 0) else 0
  lag_overall <- lag_overall[order(lag_overall$lag), , drop = FALSE]
  rownames(lag_overall) <- NULL
  
  # X-only causal graph (nodes: X variables; edges: i -> j if i GC j, weight = estimated lag)
  # NOTE: no requireNamespace() check for igraph as requested.
  edges <- as.data.frame(matrix (NA, ncol=3,nrow = 0))
  colnames(edges)<-c("from","to","weight")
  p_global <- p  # reuse global p as starting point; reduce per pair as needed
  
  for (to_idx in seq_along(safe_names)) {
    for (from_idx in seq_along(safe_names)) {
      if (from_idx == to_idx) next
      a <- safe_names[from_idx]
      b <- safe_names[to_idx]
      df_pair <- data.frame(A = X[[a]], B = X[[b]])
      names(df_pair) <- c(a, b)
      df_pair <- df_pair[stats::complete.cases(df_pair), , drop = FALSE]
      n_pair <- nrow(df_pair)
      if (n_pair < 10) next
      
      p_pair <- max(1L, min(p_global, floor((n_pair - 5)/4)))
      fit <- NULL; p_try <- 1
      BF<-0
      pval_ab<-1
      LL<-TRUE
      while (LL) {
        fit_try <- vars::VAR(df_pair, p = p_try, type = "const")
        gc_ab <-vars::causality(fit_try, cause = a)
        pval_ab <- as.numeric(gc_ab$Granger$p.value)
        if (is.na(pval_ab)) pval_ab<-1
        BF<-bayesian_granger_causality(X[[a]],X[[b]],lag=p_try)$bayes_factor
        fit <- fit_try
        
        if (p_try>max_lag) LL<-FALSE
        if ((pval_ab<alpha)&&(BF>BF10)) LL<-FALSE  
        if (LL) p_try<-p_try+1
      }
      
      if (p_try>max_lag) {
        best_lag_ab <- 0
      }else{
        p_try<-fit$p
        best_lag_ab <- p_try
      }
      
      
      if (best_lag_ab > 0) {
        edges[nrow(edges)+1,1]<-orig_names[from_idx]
        edges[nrow(edges),2]<-orig_names[to_idx]
        edges[nrow(edges),3]<-best_lag_ab
      }
    }
  }
  
  explanatory_graph <- NULL
  if (nrow(edges) > 0) {
    explanatory_graph <- igraph::graph_from_data_frame(edges,
      directed = TRUE,
      vertices = data.frame(name = orig_names, stringsAsFactors = FALSE)
    )
  } else {
    explanatory_graph <- igraph::graph_from_data_frame(
      data.frame(from = character(0), to = character(0), weight = integer(0), stringsAsFactors = FALSE),
      directed = TRUE,
      vertices = data.frame(name = orig_names, stringsAsFactors = FALSE)
    )
  }
  
  out <- list(
    variable_importance = vi,
    lag_importance = lag_imp_df[order(lag_imp_df$variable, -lag_imp_df$importance, lag_imp_df$lag), , drop = FALSE],
    lag_importance_overall = lag_overall,
    y_lag_importance = y_imp,
    expected_lags = results,
    explanatory_graph = explanatory_graph,
    details = list(
      n = n,
      n_effective = n_eff,
      selected_p_global = p,
      alpha = alpha,
      lag_select = lag_select,
      include_y_lags = include_y_lags,
      impute = impute,
      impute_method = impute_method
    )
  )
  class(out) <- c("granger_vars_result", class(out))
  out
}
