# forecast_backtest.R
#
# Out-of-sample backtest for the ARIMA / Bayesian-ARIMA forecasting
# exercise (see autoarima.R, forecasts.R), addressing Reviewer 1, Point 9
# of The_World_Economy_Rev.md ("the paper does not provide an evaluation
# of forecast accuracy").
#
# Holds out the final `h_holdout` years of each coefficient series, refits
# autoarima()'s underlying model on the remaining (truncated) years only,
# forecasts forward h_holdout steps, and compares the forecasts to the
# actually-observed held-out values using MAE and RMSE (both in original
# units and normalized by the in-sample standard deviation, so that series
# with very different scales remain comparable in one table).
#
# This reuses forecast::auto.arima() / bayesforecast::auto.sarima() exactly
# as already called by autoarima.R; it only wraps them in a train/test
# split. It mirrors the blocked, expanding-window logic already used for
# the Random-Forest causality tests (see granger_rf.R's oos_folds), so no
# new forecasting methodology is introduced -- only a standard backtest on
# top of the existing one.
#
# This script has not been executed against your actual data (we do not
# have access to it). Please sanity-check on one series first, and note
# the caveat below about bayesforecast's return object field names.
# ---------------------------------------------------------------------

forecast_backtest <- function(data, start = 1995, end = 2020,
                               h_holdout = 5, bayesian = FALSE) {

  if (!requireNamespace("forecast", quietly = TRUE)) stop("Package 'forecast' is required.")
  if (bayesian && !requireNamespace("bayesforecast", quietly = TRUE)) {
    stop("Package 'bayesforecast' is required for bayesian = TRUE.")
  }
  if (is.null(dim(data))) data <- matrix(data, ncol = 1)

  n_years <- end - start + 1
  if (nrow(data) != n_years) {
    stop(sprintf("nrow(data) = %d does not match end - start + 1 = %d.", nrow(data), n_years))
  }
  if (h_holdout >= n_years - 3) {
    stop("h_holdout is too large relative to the series length; leave at least 3 in-sample years.")
  }
  train_end <- end - h_holdout

  results <- list()
  for (i in 1:ncol(data)) {
    series_name <- colnames(data)[i]
    if (is.null(series_name)) series_name <- paste0("V", i)

    train_ts       <- ts(data[1:(train_end - start + 1), i], start = start, end = train_end)
    actual_holdout <- data[(train_end - start + 2):n_years, i]

    fc_mean <- tryCatch({
      if (bayesian) {
        fit <- bayesforecast::auto.sarima(train_ts, seasonal = FALSE)
        fc  <- bayesforecast::forecast(fit, h = h_holdout)
        # CAVEAT: bayesforecast forecast objects may not always expose a
        # "$mean" field under this exact name across package versions.
        # Run `str(fc)` / `names(fc)` once on a single series and adjust
        # the extraction below (e.g., fc$median, fc$forecast$mean, ...)
        # if this does not work out of the box.
        as.numeric(fc$mean)
      } else {
        fit <- forecast::auto.arima(train_ts)
        fc  <- forecast::forecast(fit, h = h_holdout)
        as.numeric(fc$mean)
      }
    }, error = function(e) {
      warning(sprintf("Series %s: forecasting failed (%s); skipping.", series_name, conditionMessage(e)))
      rep(NA_real_, h_holdout)
    })

    if (length(fc_mean) != length(actual_holdout) || all(is.na(fc_mean))) {
      warning(sprintf("Series %s: forecast length/content mismatch; skipping.", series_name))
      next
    }

    err <- actual_holdout - fc_mean
    mae  <- mean(abs(err), na.rm = TRUE)
    rmse <- sqrt(mean(err^2, na.rm = TRUE))
    in_sample_sd <- sd(train_ts)

    results[[series_name]] <- data.frame(
      series    = series_name,
      mae       = mae,
      rmse      = rmse,
      mae_norm  = if (in_sample_sd > 0) mae  / in_sample_sd else NA_real_,
      rmse_norm = if (in_sample_sd > 0) rmse / in_sample_sd else NA_real_,
      h_holdout = h_holdout,
      bayesian  = bayesian,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}

compare_backtests <- function(bt_traditional, bt_bayesian) {
  merged <- merge(bt_traditional, bt_bayesian, by = "series",
                   suffixes = c("_arima", "_bayes"))
  merged$rmse_norm_improvement_pct <- 100 * (merged$rmse_norm_arima - merged$rmse_norm_bayes) / merged$rmse_norm_arima
  merged[order(-merged$rmse_norm_improvement_pct), , drop = FALSE]
}

# ---------------------------------------------------------------------
# Suggested usage, mirroring the existing pipeline's variable names
# (replace `your_coefficient_matrix` with whichever already-loaded object
# holds the year-by-series coefficient series in a (years x series)
# layout, e.g. the same object currently passed into autoarima()):
#
#   bt_arima      <- forecast_backtest(your_coefficient_matrix, start = 1995,
#                                       end = 2020, h_holdout = 5, bayesian = FALSE)
#   bt_bayesarima <- forecast_backtest(your_coefficient_matrix, start = 1995,
#                                       end = 2020, h_holdout = 5, bayesian = TRUE)
#   comparison    <- compare_backtests(bt_arima, bt_bayesarima)
#   print(comparison)
#
# Repeat with end = 2022, h_holdout = 5 (i.e. train through 2017) for the
# CEPII-BACI series. Report the resulting normalized RMSE (and the
# traditional-vs-Bayesian improvement) as the new "Forecast accuracy"
# paragraph recommended in The_World_Economy_Rev.md, Point 9 -- this is
# the only genuinely new *computation* requested for that point.
# ---------------------------------------------------------------------
