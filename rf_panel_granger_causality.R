rf_panel_granger_causality <- function(Y, X, max_lag = 5, n_shadow = 5, 
                                       impute = TRUE, ntree = 500,
                                       importance_type = "permutation",
                                       seed = NULL) {
  
  # Load required packages
  if (!requireNamespace("imputeTS", quietly = TRUE)) {
    stop("Package 'imputeTS' is required. Please install it using: install.packages('imputeTS')")
  }
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("Package 'randomForest' is required. Please install it using: install.packages('randomForest')")
  }
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
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
  
  if (n_shadow < 1) {
    stop("n_shadow must be at least 1")
  }
  
  if (ntree < 1) {
    stop("ntree must be at least 1")
  }
  
  if (!importance_type %in% c("permutation", "impurity")) {
    stop("importance_type must be either 'permutation' or 'impurity'")
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
  
  # Function to create shadow features (permuted versions)
  create_shadow_features <- function(data_df, n_shadow) {
    shadow_df <- data.frame(matrix(0, nrow = nrow(data_df), ncol = 0))
    
    for (col_name in names(data_df)) {
      if (col_name != "y") {  # Don't create shadows for the target variable
        for (s in 1:n_shadow) {
          shadow_col_name <- paste0(col_name, "_shadow", s)
          shadow_df[[shadow_col_name]] <- sample(data_df[[col_name]])
        }
      }
    }
    
    return(shadow_df)
  }
  
  # Function to perform Random Forest Granger test for a given lag
  rf_granger_test_lag <- function(y, x, lag, n_shadow, ntree, importance_type) {
    tryCatch({
      # Validate series length
      n <- length(y)
      if (lag >= n - 2) {
        return(list(x_importance = NA, max_shadow_importance = NA))
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
      if (nrow(data_df) < max(10, lag + 5)) {
        return(list(x_importance = NA, max_shadow_importance = NA))
      }
      
      # Create shadow features
      shadow_df <- create_shadow_features(data_df[, -1, drop = FALSE], n_shadow)
      
      # Combine original and shadow features
      full_data <- cbind(data_df, shadow_df)
      
      # Train Random Forest
      if (importance_type == "permutation") {
        rf_model <- randomForest::randomForest(
          y ~ ., 
          data = full_data, 
          ntree = ntree,
          importance = TRUE,
          mtry = max(1, floor(sqrt(ncol(full_data) - 1)))
        )
        # Use %IncMSE for permutation importance
        importance_scores <- randomForest::importance(rf_model)[, "%IncMSE"]
      } else {
        rf_model <- randomForest::randomForest(
          y ~ ., 
          data = full_data, 
          ntree = ntree,
          importance = TRUE,
          mtry = max(1, floor(sqrt(ncol(full_data) - 1)))
        )
        # Use IncNodePurity for impurity importance
        importance_scores <- randomForest::importance(rf_model)[, "IncNodePurity"]
      }
      
      # Extract importance of x lags
      x_lag_names <- paste0("x_lag", 1:lag)
      x_importance <- mean(importance_scores[x_lag_names], na.rm = TRUE)
      
      # Extract importance of shadow features
      shadow_names <- grep("_shadow", names(importance_scores), value = TRUE)
      shadow_importance <- importance_scores[shadow_names]
      max_shadow_importance <- max(shadow_importance, na.rm = TRUE)
      
      return(list(
        x_importance = x_importance,
        max_shadow_importance = max_shadow_importance
      ))
      
    }, error = function(e) {
      warning(paste("Error in rf_granger_test_lag:", e$message))
      return(list(x_importance = NA, max_shadow_importance = NA))
    })
  }
  
  # Function to perform Random Forest panel Granger test on stacked data
  rf_panel_granger_test <- function(Y_matrix, X_matrix, lag, n_shadow, ntree, importance_type) {
    tryCatch({
      n_years <- nrow(Y_matrix)
      n_sectors <- ncol(Y_matrix)
      
      if (lag >= n_years - 2) {
        return(list(x_importance = NA, max_shadow_importance = NA))
      }
      
      # Create panel data structure
      panel_data <- data.frame()
      
      for (j in 1:n_sectors) {
        for (t in (lag + 1):n_years) {
          row_data <- data.frame(
            sector = j,
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
      if (nrow(panel_data) < max(20, lag + 10)) {
        return(list(x_importance = NA, max_shadow_importance = NA))
      }
      
      # Convert sector to factor
      panel_data$sector <- as.factor(panel_data$sector)
      
      # Create shadow features (excluding y and sector)
      predictor_cols <- setdiff(names(panel_data), c("y", "sector"))
      shadow_df <- data.frame(matrix(0, nrow = nrow(panel_data), ncol = 0))
      
      for (col_name in predictor_cols) {
        for (s in 1:n_shadow) {
          shadow_col_name <- paste0(col_name, "_shadow", s)
          shadow_df[[shadow_col_name]] <- sample(panel_data[[col_name]])
        }
      }
      
      # Combine original and shadow features
      full_panel_data <- cbind(panel_data, shadow_df)
      
      # Train Random Forest
      if (importance_type == "permutation") {
        rf_model <- randomForest::randomForest(
          y ~ ., 
          data = full_panel_data, 
          ntree = ntree,
          importance = TRUE,
          mtry = max(1, floor(sqrt(ncol(full_panel_data) - 1)))
        )
        importance_scores <- randomForest::importance(rf_model)[, "%IncMSE"]
      } else {
        rf_model <- randomForest::randomForest(
          y ~ ., 
          data = full_panel_data, 
          ntree = ntree,
          importance = TRUE,
          mtry = max(1, floor(sqrt(ncol(full_panel_data) - 1)))
        )
        importance_scores <- randomForest::importance(rf_model)[, "IncNodePurity"]
      }
      
      # Extract importance of x lags
      x_lag_names <- paste0("x_lag", 1:lag)
      x_importance <- mean(importance_scores[x_lag_names], na.rm = TRUE)
      
      # Extract importance of shadow features
      shadow_names <- grep("_shadow", names(importance_scores), value = TRUE)
      shadow_importance <- importance_scores[shadow_names]
      max_shadow_importance <- max(shadow_importance, na.rm = TRUE)
      
      return(list(
        x_importance = x_importance,
        max_shadow_importance = max_shadow_importance
      ))
      
    }, error = function(e) {
      warning(paste("Error in rf_panel_granger_test:", e$message))
      return(list(x_importance = NA, max_shadow_importance = NA))
    })
  }
  
  # Test for each sector individually
  cat("Testing individual sectors using Random Forest...\n")
  cat(sprintf("Parameters: ntree=%d, n_shadow=%d, importance_type=%s\n", 
              ntree, n_shadow, importance_type))
  
  for (j in 1:n_sectors) {
    cat(sprintf("  Sector %d/%d\n", j, n_sectors))
    for (i in 1:n_explanatory) {
      y_series <- Y[, j]
      x_series <- X[i, j, ]
      
      # Store importance values for each lag
      importance_by_lag <- numeric(max_lag)
      max_shadow_by_lag <- numeric(max_lag)
      
      # Test all lags and track importance
      for (lag in 1:max_lag) {
        result <- rf_granger_test_lag(y_series, x_series, lag, n_shadow, ntree, importance_type)
        importance_by_lag[lag] <- ifelse(is.na(result$x_importance), -Inf, result$x_importance)
        max_shadow_by_lag[lag] <- ifelse(is.na(result$max_shadow_importance), Inf, result$max_shadow_importance)
      }
      
      # Find the lag with highest importance
      best_lag_idx <- which.max(importance_by_lag)
      best_importance <- importance_by_lag[best_lag_idx]
      
      # Check if best importance exceeds shadow features at that lag
      if (is.finite(best_importance) && 
          best_importance > max_shadow_by_lag[best_lag_idx] &&
          best_importance > 0) {
        result_matrix[i, j] <- best_lag_idx
        cat(sprintf("    Var_%d -> Sector_%d: Best lag = %d (Importance = %.4f > Shadow = %.4f)\n", 
                    i, j, best_lag_idx, best_importance, max_shadow_by_lag[best_lag_idx]))
      } else {
        result_matrix[i, j] <- 0
        cat(sprintf("    Var_%d -> Sector_%d: No significant lag (Importance <= Shadow)\n", i, j))
      }
    }
  }
  
  # Test for pooled data (panel approach)
  cat("\nTesting pooled panel data using Random Forest...\n")
  for (i in 1:n_explanatory) {
    # Prepare X matrix for this explanatory variable
    X_matrix <- t(X[i, , ])
    
    # Store importance values for each lag
    importance_by_lag <- numeric(max_lag)
    max_shadow_by_lag <- numeric(max_lag)
    
    # Test all lags
    for (lag in 1:max_lag) {
      result <- rf_panel_granger_test(Y, X_matrix, lag, n_shadow, ntree, importance_type)
      importance_by_lag[lag] <- ifelse(is.na(result$x_importance), -Inf, result$x_importance)
      max_shadow_by_lag[lag] <- ifelse(is.na(result$max_shadow_importance), Inf, result$max_shadow_importance)
    }
    
    # Find the lag with highest importance
    best_lag_idx <- which.max(importance_by_lag)
    best_importance <- importance_by_lag[best_lag_idx]
    
    # Check if best importance exceeds shadow features
    if (is.finite(best_importance) && 
        best_importance > max_shadow_by_lag[best_lag_idx] &&
        best_importance > 0) {
      result_matrix[i, n_sectors + 1] <- best_lag_idx
      cat(sprintf("  Var_%d -> Pooled: Best lag = %d (Importance = %.4f > Shadow = %.4f)\n", 
                  i, best_lag_idx, best_importance, max_shadow_by_lag[best_lag_idx]))
    } else {
      result_matrix[i, n_sectors + 1] <- 0
      cat(sprintf("  Var_%d -> Pooled: No significant lag (Importance <= Shadow)\n", i))
    }
  }
  
  cat("\nRandom Forest analysis complete.\n")
  cat("Note: Non-zero values indicate the lag with highest importance that exceeds shadow features.\n")
  cat("      Zero values indicate the variable does not exceed random importance at any lag.\n")
  rownames(result_matrix)<-dimnames(X)[[1]]
  return(result_matrix)
}

