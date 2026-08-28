bayesian_granger_causality <- function(x, y, lag = 1, prior_prob = 0.5) {
  
  # Input validation
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("Both x and y must be numeric vectors")
  }
  
  if (length(x) != length(y)) {
    stop("Time series x and y must have the same length")
  }
  
  if (lag < 1 || lag >= length(x) - 1) {
    stop("Lag must be between 1 and ", length(x) - 2)
  }
  
  if (prior_prob < 0 || prior_prob > 1) {
    stop("Prior probability must be between 0 and 1")
  }
  
  if (length(x) < (lag + 10)) {
    warning("Small sample size. Results may be unreliable.")
  }
  
  # Prepare data
  n <- length(x)
  effective_n <- n - lag
  
  # Create dependent variable (y values after lag period)
  y_dep <- y[(lag + 1):n]
  
  # Create lagged predictors
  # Y lags (own history)
  y_lags <- matrix(NA, nrow = effective_n, ncol = lag)
  for (i in 1:lag) {
    y_lags[, i] <- y[(lag + 1 - i):(n - i)]
  }
  colnames(y_lags) <- paste0("y_lag", 1:lag)
  
  # X lags (potential causes)  
  x_lags <- matrix(NA, nrow = effective_n, ncol = lag)
  for (i in 1:lag) {
    x_lags[, i] <- x[(lag + 1 - i):(n - i)]
  }
  colnames(x_lags) <- paste0("x_lag", 1:lag)
  
  # Create data frames for models
  restricted_data <- data.frame(y = y_dep, y_lags)
  unrestricted_data <- data.frame(y = y_dep, y_lags, x_lags)
  
  # Fit models using standard linear regression
  tryCatch({
    
    # Restricted model (H0: no causality - only y's own lags)
    if (lag == 1) {
      restricted_formula <- y ~ y_lag1
    } else {
      y_lag_terms <- paste0("y_lag", 1:lag, collapse = " + ")
      restricted_formula <- as.formula(paste("y ~", y_lag_terms))
    }
    
    model_restricted <- lm(restricted_formula, data = restricted_data)
    
    # Unrestricted model (H1: causality exists - includes x lags)
    if (lag == 1) {
      unrestricted_formula <- y ~ y_lag1 + x_lag1
    } else {
      y_lag_terms <- paste0("y_lag", 1:lag, collapse = " + ")
      x_lag_terms <- paste0("x_lag", 1:lag, collapse = " + ")
      unrestricted_formula <- as.formula(paste("y ~", y_lag_terms, "+", x_lag_terms))
    }
    
    model_unrestricted <- lm(unrestricted_formula, data = unrestricted_data)
    
    # Calculate Bayesian Information Criterion (BIC)
    bic_restricted <- BIC(model_restricted)
    bic_unrestricted <- BIC(model_unrestricted)
    
    # Calculate Bayes Factor using BIC approximation
    # BF10 = exp((BIC_H0 - BIC_H1)/2)
    log_bayes_factor <- (bic_restricted - bic_unrestricted) / 2
    bayes_factor <- exp(log_bayes_factor)
    
    # Calculate posterior probability using Bayes' theorem
    posterior_prob <- (bayes_factor * prior_prob) / (bayes_factor * prior_prob + (1 - prior_prob))
    
    # Perform classical F-test for comparison
    f_test <- anova(model_restricted, model_unrestricted)
    p_value_classical <- f_test$`Pr(>F)`[2]
    
    # Determine strength of evidence based on Jeffreys' scale
    if (bayes_factor > 100) {
      strength <- "Decisive evidence FOR causality"
      decision <- "Strong Accept"
      evidence_level <- "decisive_for"
    } else if (bayes_factor > 30) {
      strength <- "Very strong evidence FOR causality"
      decision <- "Accept"
      evidence_level <- "very_strong_for"
    } else if (bayes_factor > 10) {
      strength <- "Strong evidence FOR causality"
      decision <- "Accept"
      evidence_level <- "strong_for"
    } else if (bayes_factor > 3) {
      strength <- "Moderate evidence FOR causality"
      decision <- "Accept"
      evidence_level <- "moderate_for"
    } else if (bayes_factor > 1) {
      strength <- "Weak evidence FOR causality"
      decision <- "Weak evidence - collect more data"
      evidence_level <- "weak_for"
    } else if (bayes_factor > 1/3) {
      strength <- "Inconclusive evidence"
      decision <- "Inconclusive"
      evidence_level <- "inconclusive"
    } else if (bayes_factor > 1/10) {
      strength <- "Moderate evidence AGAINST causality"
      decision <- "Reject"
      evidence_level <- "moderate_against"
    } else if (bayes_factor > 1/30) {
      strength <- "Strong evidence AGAINST causality"
      decision <- "Reject"
      evidence_level <- "strong_against"
    } else {
      strength <- "Very strong evidence AGAINST causality"
      decision <- "Strong Reject"
      evidence_level <- "very_strong_against"
    }
    
    # Create detailed recommendation
    interpretation <- if (bayes_factor > 3) {
      "The data provides substantial evidence that x Granger-causes y. This suggests x contains information useful for predicting y beyond y's own past."
    } else if (bayes_factor > 1) {
      "The data provides some evidence for Granger causality, but it's not conclusive. Consider collecting more data or examining different lag structures."
    } else if (bayes_factor > 1/3) {
      "The evidence is inconclusive regarding Granger causality. The data doesn't strongly favor either hypothesis."
    } else {
      "The data provides evidence against Granger causality from x to y. x does not appear to contain useful information for predicting y."
    }
    
    # Model comparison details
    model_comparison <- data.frame(
      Model = c("Restricted (H0)", "Unrestricted (H1)"),
      BIC = c(bic_restricted, bic_unrestricted),
      AIC = c(AIC(model_restricted), AIC(model_unrestricted)),
      R_squared = c(summary(model_restricted)$r.squared, summary(model_unrestricted)$r.squared),
      Residual_SE = c(summary(model_restricted)$sigma, summary(model_unrestricted)$sigma)
    )
    
    # Generate full recommendation text
    recommendation_text <- paste0(
      "BAYESIAN GRANGER CAUSALITY TEST RESULTS\n",
      "========================================\n\n",
      "MAIN RESULTS:\n",
      "Bayes Factor: ", round(bayes_factor, 4), "\n",
      "Evidence Strength: ", strength, "\n",
      "Posterior Probability of Causality: ", round(posterior_prob, 4), "\n",
      "Recommendation: ", decision, "\n\n",
      "INTERPRETATION:\n",
      interpretation, "\n\n",
      "COMPARISON WITH CLASSICAL TEST:\n",
      "Classical Granger test p-value: ", round(p_value_classical, 6), "\n",
      "Classical test result: ", ifelse(p_value_classical < 0.05, "Significant (reject H0)", "Not significant (fail to reject H0)"), "\n\n",
      "MODEL COMPARISON:\n",
      "Restricted model BIC: ", round(bic_restricted, 2), "\n",
      "Unrestricted model BIC: ", round(bic_unrestricted, 2), "\n",
      "BIC difference: ", round(bic_restricted - bic_unrestricted, 2), " (positive favors causality)\n\n",
      "TECHNICAL DETAILS:\n",
      "Effective sample size: ", effective_n, "\n",
      "Number of lags: ", lag, "\n",
      "Log Bayes Factor: ", round(log_bayes_factor, 4), "\n"
    )
    
    # Return comprehensive results
    result <- list(
      # Main results
      bayes_factor = bayes_factor,
      log_bayes_factor = log_bayes_factor,
      posterior_prob = posterior_prob,
      decision = decision,
      strength = strength,
      evidence_level = evidence_level,
      
      # Additional information
      p_value_classical = p_value_classical,
      recommendation = recommendation_text,
      interpretation = interpretation,
      
      # Model details
      model_comparison = model_comparison,
      models = list(
        restricted = model_restricted,
        unrestricted = model_unrestricted
      ),
      
      # Data information
      sample_info = list(
        total_observations = n,
        effective_observations = effective_n,
        lag_order = lag,
        prior_prob = prior_prob
      )
    )
    
    return(result)
    
  }, error = function(e) {
    stop(paste("Error in model fitting:", e$message, "\nCheck your data for missing values or other issues."))
  })
}
