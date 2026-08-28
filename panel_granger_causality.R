panel_granger_causality <- function(Y, X, max_lag = 5, p_value = 0.01, impute = TRUE) {
  
  # Load required packages
  if (!requireNamespace("imputeTS", quietly = TRUE)) {
    stop("Package 'imputeTS' is required. Please install it using: install.packages('imputeTS')")
  }
  if (!requireNamespace("lmtest", quietly = TRUE)) {
    stop("Package 'lmtest' is required. Please install it using: install.packages('lmtest')")
  }
  
  # Validate inputs
  if (!is.matrix(Y)) {
    stop("Y must be a matrix")
  }
  if (!is.array(X) || length(dim(X)) != 3) {
    stop("X must be a 3-dimensional array")
  }
  
  n_years <- nrow(Y)
  n_sectors <- ncol(Y)
  n_explanatory <- dim(X)[1]
  
  # Check dimensions
  if (dim(X)[2] != n_sectors) {
    stop("The number of sectors in X (dimension 2) must match the number of columns in Y")
  }
  if (dim(X)[3] != n_years) {
    stop("The number of years in X (dimension 3) must match the number of rows in Y")
  }
  
  if (max_lag >= n_years) {
    stop("max_lag must be less than the number of years")
  }
  
  if (p_value <= 0 || p_value >= 1) {
    stop("p_value must be between 0 and 1")
  }
  
  # Perform imputation if requested
  if (impute) {
    # Impute Y (by column/sector)
    for (j in 1:n_sectors) {
      if (any(is.na(Y[, j]))) {
        Y[, j] <- imputeTS::na_interpolation(Y[, j], option = "linear")
      }
    }
    
    # Impute X (by explanatory variable and sector)
    for (i in 1:n_explanatory) {
      for (j in 1:n_sectors) {
        if (any(is.na(X[i, j, ]))) {
          X[i, j, ] <- imputeTS::na_interpolation(X[i, j, ], option = "linear")
        }
      }
    }
  }
  
  # Check for remaining NA values
  if (any(is.na(Y)) || any(is.na(X))) {
    warning("NA values remain after imputation. Results may be unreliable.")
  }
  
  # Initialize result matrix
  # Columns: sectors + 1 for pooled data
  # Rows: explanatory variables
  result_matrix <- matrix(0, nrow = n_explanatory, ncol = n_sectors + 1)
  colnames(result_matrix) <- c(paste0("Sector_", 1:n_sectors), "Pooled")
  rownames(result_matrix) <- paste0("Var_", 1:n_explanatory)
  
  # Function to perform Granger causality test for a given lag on a single time series
  granger_test_lag <- function(y, x, lag) {
    tryCatch({
      # Validate series length
      n <- length(y)
      if (lag >= n - 1) {
        return(NA)
      }
      
      # Create data frame with lagged values
      data_df <- data.frame(y = y[(lag + 1):n])
      
      # Add lags of y
      for (l in 1:lag) {
        data_df[[paste0("y_lag", l)]] <- y[(lag + 1 - l):(n - l)]
      }
      
      # Add lags of x
      for (l in 1:lag) {
        data_df[[paste0("x_lag", l)]] <- x[(lag + 1 - l):(n - l)]
      }
      
      # Remove any rows with NA
      data_df <- na.omit(data_df)
      
      # Check if we have enough observations
      if (nrow(data_df) < lag + 2) {
        return(NA)
      }
      
      # Restricted model (without x lags)
      formula_restricted <- as.formula(paste("y ~", paste(paste0("y_lag", 1:lag), collapse = " + ")))
      model_restricted <- lm(formula_restricted, data = data_df)
      
      # Unrestricted model (with x lags)
      formula_unrestricted <- as.formula(paste("y ~", 
                                               paste(c(paste0("y_lag", 1:lag), paste0("x_lag", 1:lag)), 
                                                     collapse = " + ")))
      model_unrestricted <- lm(formula_unrestricted, data = data_df)
      
      # Perform F-test (Wald test) to compare models
      test_result <- lmtest::waldtest(model_restricted, model_unrestricted, test = "F")
      p_val <- test_result$`Pr(>F)`[2]
      
      return(p_val)
    }, error = function(e) {
      warning(paste("Error in granger_test_lag:", e$message))
      return(NA)
    })
  }
  
  # Function to perform panel Granger test on stacked data
  panel_granger_test <- function(Y_matrix, X_matrix, lag) {
    tryCatch({
      n_years <- nrow(Y_matrix)
      n_sectors <- ncol(Y_matrix)
      
      if (lag >= n_years - 1) {
        return(NA)
      }
      
      # Create panel data structure
      panel_data <- data.frame()
      
      for (j in 1:n_sectors) {
        for (t in (lag + 1):n_years) {
          row_data <- data.frame(
            sector = j,
            time = t,
            y = Y_matrix[t, j]
          )
          
          # Add lags of y
          for (l in 1:lag) {
            row_data[[paste0("y_lag", l)]] <- Y_matrix[t - l, j]
          }
          
          # Add lags of x
          for (l in 1:lag) {
            row_data[[paste0("x_lag", l)]] <- X_matrix[t - l, j]
          }
          
          panel_data <- rbind(panel_data, row_data)
        }
      }
      
      # Remove any rows with NA
      panel_data <- na.omit(panel_data)
      
      # Check if we have enough observations
      if (nrow(panel_data) < lag + 2) {
        return(NA)
      }
      
      # Restricted model (without x lags) - includes sector fixed effects
      formula_restricted <- as.formula(paste("y ~", 
                                            paste(paste0("y_lag", 1:lag), collapse = " + "),
                                            "+ factor(sector)"))
      model_restricted <- lm(formula_restricted, data = panel_data)
      
      # Unrestricted model (with x lags) - includes sector fixed effects
      formula_unrestricted <- as.formula(paste("y ~", 
                                               paste(c(paste0("y_lag", 1:lag), paste0("x_lag", 1:lag)), 
                                                     collapse = " + "),
                                               "+ factor(sector)"))
      model_unrestricted <- lm(formula_unrestricted, data = panel_data)
      
      # Perform F-test
      test_result <- lmtest::waldtest(model_restricted, model_unrestricted, test = "F")
      p_val <- test_result$`Pr(>F)`[2]
      
      return(p_val)
    }, error = function(e) {
      warning(paste("Error in panel_granger_test:", e$message))
      return(NA)
    })
  }
  
  # Test for each sector individually
  #cat("Testing individual sectors...\n")
  for (j in 1:n_sectors) {
    for (i in 1:n_explanatory) {
      y_series <- Y[, j]
      x_series <- X[i, j, ]
      
      # Test different lags, starting from lag 1
      best_lag <- 0
      for (lag in 1:max_lag) {
        p_val <- granger_test_lag(y_series, x_series, lag)
        
        if (!is.na(p_val) && p_val < p_value) {
          best_lag <- lag
          break  # Found the minimum significant lag
        }
      }
      
      result_matrix[i, j] <- best_lag
    }
  }
  
  # Test for pooled data (panel approach)
  #cat("Testing pooled panel data...\n")
  for (i in 1:n_explanatory) {
    # Prepare X matrix for this explanatory variable
    # X[i, , ] has dimensions: sectors x years
    # We need years x sectors, so transpose
    X_matrix <- t(X[i, , ])
    
    # Test different lags
    best_lag <- 0
    for (lag in 1:max_lag) {
      p_val <- panel_granger_test(Y, X_matrix, lag)
      
      if (!is.na(p_val) && p_val < p_value) {
        best_lag <- lag
        break
      }
    }
    
    result_matrix[i, n_sectors + 1] <- best_lag
  }
  
  #cat("Analysis complete.\n")
  rownames(result_matrix)<-dimnames(X)[[1]]
  return(result_matrix)
}

