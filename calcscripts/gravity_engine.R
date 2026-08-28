## =====================================================================
##  gravity_engine.R
##  Shared estimation engine for the yearly gravity-model scripts.
##
##  Two model variants (chosen by the caller via MODEL = "OLS" / "PPML"):
##    * OLS  : log-log gravity. Dependent = log(Vol); the continuous
##             regressors (distance, GDP_i, GDP_j) are log-transformed.
##    * PPML : Poisson pseudo-maximum-likelihood gravity
##             (Santos Silva & Tenreyro, 2006). Dependent = Vol in LEVELS
##             (handles zeros natively); the continuous regressors are
##             log-transformed, so coefficients are elasticities and are
##             directly comparable with the OLS estimates.
##
##  Inference : heteroskedasticity-robust (HC1) standard errors.
##  Backends  : uses 'fixest' (feols / fepois) when installed for speed and
##              robust SEs; otherwise falls back to base R (lm / glm) with
##              'sandwich'::vcovHC, and finally to classical SEs. The scripts
##              therefore run in any environment.
##  Robustness: every estimation is wrapped in tryCatch and returns NA on any
##              failure (too few observations, constant / collinear dummies,
##              non-convergence, ...). No Inf / NaN ever reaches the output.
##
##  NOTE ON FIXED EFFECTS: per the agreed setup NO country fixed effects are
##  used. In a *single-year* regression, exporter / importer dummies would be
##  collinear with GDP_i / GDP_j and would absorb the GDP effect. If FE are
##  ever required, they must be estimated separately (and GDP dropped); the
##  results should then be stored in a separate object/file.
## =====================================================================

## Output coefficient order (12 slopes; "adjreg" = model fit appended last).
COEF_NAMES <- c("dist", "GDPi", "GDPj", "comlang_off", "comlang_ethno",
                "contig", "col45", "col_dep_ever", "comleg_pretrans",
                "fta_wto", "eu_o", "eu_d")

## Control (dummy) variables: used as-is, NEVER log-transformed.
DUMMY_NAMES <- c("comlang_off", "comlang_ethno", "contig", "col45",
                 "col_dep_ever", "comleg_pretrans", "fta_wto", "eu_o", "eu_d")

## Full column set of a finished result vector (coefficients + fit).
OUT_COLS <- c(COEF_NAMES, "adjreg")

## Optional backends (detected once; the engine still runs if absent).
HAVE_FIXEST   <- requireNamespace("fixest",   quietly = TRUE)
HAVE_SANDWICH <- requireNamespace("sandwich", quietly = TRUE)

## Coerce a column of any type (numeric / integer / logical / factor /
## character) to numeric without throwing an error.
to_num <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  if (is.logical(x)) return(as.numeric(x))
  if (is.factor(x))  return(suppressWarnings(as.numeric(as.character(x))))
  suppressWarnings(as.numeric(x))
}

## Positivity-safe logarithm implementing the requested zero handling:
##   x  > 0            -> log(x)
##   x <= 0 (e.g. 0)   -> treated as 1  ->  log(1) = 0
##   NA / +-Inf / NaN  -> NA            (never returns +-Inf)
safe_log <- function(x) {
  x   <- to_num(x)
  out <- rep(NA_real_, length(x))
  fin <- is.finite(x)
  out[fin & x >  0] <- log(x[fin & x > 0])
  out[fin & x <= 0] <- 0
  out
}

## Replace any non-finite value (Inf, -Inf, NaN) with NA (no Inf in output).
na_if_nonfinite <- function(v) { v[!is.finite(v)] <- NA_real_; v }

## All-NA result placeholder.
empty_result <- function() {
  b <- stats::setNames(rep(NA_real_, length(COEF_NAMES)), COEF_NAMES)
  list(beta = b, p = b, fit = NA_real_, model_p = NA_real_)
}

## Load a single named object from an .RData file (robust to a differing
## object name: falls back to the first object stored in the file).
load_one <- function(path, varname = NULL) {
  e <- new.env()
  load(path, envir = e)
  nms <- ls(e)
  if (!is.null(varname) && varname %in% nms) return(get(varname, envir = e))
  if (length(nms) >= 1L) return(get(nms[1L], envir = e))
  stop("No object found in ", path)
}

## Build a model frame whose column names already equal the OUTPUT names.
build_model_frame <- function(d) {
  v <- to_num(d$Vol)
  v[!is.finite(v) | v < 0] <- NA_real_          # PPML response must be >= 0
  mf <- data.frame(
    y_ols  = safe_log(d$Vol),                    # OLS  response: log(Vol)
    y_ppml = v,                                  # PPML response: Vol (levels)
    dist   = safe_log(d$Dist_ij),
    GDPi   = safe_log(d$GDP_i),
    GDPj   = safe_log(d$GDP_j),
    stringsAsFactors = FALSE
  )
  for (nm in DUMMY_NAMES) mf[[nm]] <- to_num(d[[nm]])
  mf
}

## Predictors that actually vary in this cell (constant ones are unusable).
.usable_predictors <- function(mf, preds) {
  preds[vapply(preds, function(p) {
    vv <- stats::var(mf[[p]]); is.finite(vv) && vv > 0
  }, logical(1))]
}

## Robust vcov extraction for a fixest model, tolerant to package version.
.fixest_robust_vcov <- function(m) {
  V <- tryCatch(stats::vcov(m, vcov = "hetero"), error = function(e) NULL)
  if (is.null(V)) V <- tryCatch(stats::vcov(m, se = "hetero"), error = function(e) NULL)
  if (is.null(V)) V <- stats::vcov(m)
  V
}

## Low-level fit; returns a standardized list or NULL on failure.
.fit_model <- function(fml, mf, model) {
  yv <- all.vars(fml)[1]
  if (model == "OLS") {
    if (HAVE_FIXEST) {
      m <- fixest::feols(fml, data = mf, warn = FALSE, notes = FALSE)
      list(co = stats::coef(m), V = .fixest_robust_vcov(m),
           fitted = as.numeric(stats::predict(m)), y = mf[[yv]], n = nrow(mf))
    } else {
      m <- stats::lm(fml, data = mf)
      V <- tryCatch(if (HAVE_SANDWICH) sandwich::vcovHC(m, type = "HC1")
                    else stats::vcov(m), error = function(e) stats::vcov(m))
      list(co = stats::coef(m), V = V,
           fitted = as.numeric(stats::fitted(m)), y = mf[[yv]], n = nrow(mf))
    }
  } else {
    if (HAVE_FIXEST) {
      m <- fixest::fepois(fml, data = mf, warn = FALSE, notes = FALSE)
      list(co = stats::coef(m), V = .fixest_robust_vcov(m),
           fitted = as.numeric(stats::predict(m)), y = mf[[yv]], n = nrow(mf))
    } else {
      m <- suppressWarnings(stats::glm(
             fml, data = mf, family = stats::poisson(link = "log"),
             control = stats::glm.control(maxit = 100)))
      if (!isTRUE(m$converged)) return(NULL)
      V <- tryCatch(if (HAVE_SANDWICH) sandwich::vcovHC(m, type = "HC1")
                    else stats::vcov(m), error = function(e) stats::vcov(m))
      list(co = stats::coef(m), V = V,
           fitted = as.numeric(stats::fitted(m)), y = mf[[yv]], n = nrow(mf))
    }
  }
}

## Turn a fitted model into the required coefficient / p-value / fit outputs.
.finalize_fit <- function(fit, model) {
  co <- fit$co; V <- fit$V; n <- fit$n
  co <- co[is.finite(co)]
  if (length(co) == 0L) return(NULL)
  keep_v <- rownames(V) %in% names(co)
  V <- V[keep_v, keep_v, drop = FALSE]
  slopes <- setdiff(names(co), "(Intercept)")
  slopes <- slopes[slopes %in% rownames(V)]
  if (length(slopes) == 0L) return(NULL)

  b <- co[slopes]
  se_all <- sqrt(diag(V)); names(se_all) <- rownames(V)
  se <- se_all[slopes]
  ok <- is.finite(b) & is.finite(se) & se > 0

  pv <- stats::setNames(rep(NA_real_, length(slopes)), slopes)
  if (model == "OLS") {
    dfres <- n - length(co)
    if (dfres > 0 && any(ok))
      pv[ok] <- 2 * stats::pt(-abs(b[ok] / se[ok]), df = dfres)
  } else {
    if (any(ok))
      pv[ok] <- 2 * stats::pnorm(-abs(b[ok] / se[ok]))
  }

  if (model == "OLS") {
    ss_res <- sum((fit$y - fit$fitted)^2)
    ss_tot <- sum((fit$y - mean(fit$y))^2)
    r2 <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
    k  <- length(slopes)
    fitv <- if (is.finite(r2) && (n - k - 1) > 0)
              1 - (1 - r2) * (n - 1) / (n - k - 1) else NA_real_
  } else {
    fitv <- tryCatch(suppressWarnings(
              stats::cor(fit$y, fit$fitted, use = "complete.obs"))^2,
              error = function(e) NA_real_)
  }

  Vs <- V[slopes, slopes, drop = FALSE]
  W  <- tryCatch(as.numeric(t(b) %*% solve(Vs) %*% b), error = function(e) NA_real_)
  model_p <- if (is.finite(W) && W >= 0)
               stats::pchisq(W, df = length(slopes), lower.tail = FALSE) else NA_real_

  list(coef = b, pval = pv, fit = fitv, model_p = model_p)
}

## Core estimation for one data subset ("cell").
.estimate_core <- function(d, model) {
  if (is.null(d) || nrow(d) == 0L) return(empty_result())
  mf <- build_model_frame(d)
  yvar  <- if (model == "OLS") "y_ols" else "y_ppml"
  preds <- c("dist", "GDPi", "GDPj", DUMMY_NAMES)

  bad  <- rowSums(vapply(mf[preds], function(z) !is.finite(z),
                         logical(nrow(mf))))
  keep <- is.finite(mf[[yvar]]) & (bad == 0)
  mf <- mf[keep, , drop = FALSE]
  n  <- nrow(mf)
  if (n < 3L) return(empty_result())

  usable <- .usable_predictors(mf, preds)
  if (length(usable) == 0L)       return(empty_result())
  if (n <= (length(usable) + 1L)) return(empty_result())

  fml <- stats::as.formula(paste(yvar, "~", paste(usable, collapse = " + ")))
  fit <- .fit_model(fml, mf, model)
  if (is.null(fit)) return(empty_result())
  est <- .finalize_fit(fit, model)
  if (is.null(est)) return(empty_result())

  out <- empty_result()
  out$beta[names(est$coef)] <- est$coef
  out$p[names(est$pval)]    <- est$pval
  out$fit     <- est$fit
  out$model_p <- est$model_p
  out
}

## Public entry point: fully error-tolerant; guarantees NA (never Inf).
estimate_cell <- function(d, model = c("OLS", "PPML")) {
  model <- match.arg(model)
  res <- tryCatch(.estimate_core(d, model), error = function(e) empty_result())
  res$beta    <- na_if_nonfinite(res$beta)
  res$p       <- na_if_nonfinite(res$p)
  res$fit     <- na_if_nonfinite(res$fit)
  res$model_p <- na_if_nonfinite(res$model_p)
  res
}

## Convenience: build the length-13 output vectors (coefficients + fit).
cell_beta_vec <- function(r) c(r$beta[COEF_NAMES], adjreg = unname(r$fit))
cell_pval_vec <- function(r) c(r$p[COEF_NAMES],    adjreg = unname(r$model_p))
