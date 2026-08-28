## =====================================================================
##  CEPII_PROD_PPML.R                                     (script 10 of 10)
##  PPML (Poisson pseudo-maximum-likelihood) gravity regressions on the RAW
##  (non-aggregated) 6-digit CEPII BACI data. One regression per
##  (product, year). Dependent = Vol in LEVELS; continuous regressors logged.
##
##  The raw data are too large for a single .RData, so they are stored one
##  file per year in a sub-directory (DYADIC_6DIGIT_YEARS/):
##    DYADIC_BACI_6DIGIT_1995.RData ... DYADIC_BACI_6DIGIT_2023.RData
##  each holding one object 'DYADIC_BACI_6DIGIT_YEAR'.
##  Files are read ONE YEAR AT A TIME to keep the memory footprint low.
##
##  Results: two 3D arrays saved to '3D_grav_products_time_ALL_PPML.RData':
##    var_products_time_ALL   : beta values
##    var_products_time_ALL_p : p-values
##  Dimensions of each array:
##    [[1]] coefficients (12) + "adjreg"   (fit / model p-value)
##    [[2]] products (6-digit; union across all processed years)
##    [[3]] years
## =====================================================================

## ---- USER CONFIGURATION (adjust the paths to your environment) -------
MODEL       <- "PPML"
input_dir   <- "DYADIC_6DIGIT_YEARS"                 # folder with yearly files
file_tmpl   <- "DYADIC_BACI_6DIGIT_%d.RData"          # %d = year
input_var   <- "DYADIC_BACI_6DIGIT_YEAR"              # object name inside each
year_range  <- 1995:2023                              # years to look for
output_file <- "3D_grav_products_time_ALL_PPML.RData"

## ---- Load the shared estimation engine -------------------------------
.here <- "gravity_engine.R"
if (!file.exists(.here)) .here <- file.path("calcscripts", "gravity_engine.R")
if (!file.exists(.here)) stop("gravity_engine.R not found; keep it next to this script.")
source(.here)

## ---- Pass 1: estimate year by year, keep per-year coef/p matrices ----
beta_by_year <- list()
p_by_year    <- list()
years_done   <- character(0)

for (yr in year_range) {
  f <- file.path(input_dir, sprintf(file_tmpl, yr))
  if (!file.exists(f)) { message("Missing file, skipped: ", f); next }

  DAT <- as.data.frame(load_one(f, input_var), stringsAsFactors = FALSE)
  prods <- sort(unique(as.character(DAT$product)))

  bm <- matrix(NA_real_, length(OUT_COLS), length(prods), dimnames = list(OUT_COLS, prods))
  pm <- matrix(NA_real_, length(OUT_COLS), length(prods), dimnames = list(OUT_COLS, prods))

  idx_p <- split(seq_len(nrow(DAT)), factor(as.character(DAT$product), levels = prods))
  for (p in prods) {
    rp <- idx_p[[p]]
    if (length(rp) == 0L) next
    r <- estimate_cell(DAT[rp, , drop = FALSE], MODEL)
    bm[, p] <- cell_beta_vec(r)
    pm[, p] <- cell_pval_vec(r)
  }

  beta_by_year[[as.character(yr)]] <- bm
  p_by_year[[as.character(yr)]]    <- pm
  years_done <- c(years_done, as.character(yr))
  message("Year ", yr, " done (", length(prods), " products).")
  rm(DAT, bm, pm, idx_p); gc(verbose = FALSE)
}

if (length(years_done) == 0L) stop("No yearly files found in ", input_dir)

## ---- Pass 2: union of products, assemble the 3D array ---------------
all_products <- sort(unique(unlist(lapply(beta_by_year, colnames), use.names = FALSE)))
years <- years_done

dims <- c(length(OUT_COLS), length(all_products), length(years))
dn   <- list(OUT_COLS, all_products, years)
var_products_time_ALL   <- array(NA_real_, dim = dims, dimnames = dn)
var_products_time_ALL_p <- array(NA_real_, dim = dims, dimnames = dn)

for (yr in years) {
  bm <- beta_by_year[[yr]]; pm <- p_by_year[[yr]]
  var_products_time_ALL[, colnames(bm), yr]   <- bm
  var_products_time_ALL_p[, colnames(pm), yr] <- pm
}

## ---- Save ------------------------------------------------------------
save(var_products_time_ALL, var_products_time_ALL_p, file = output_file)
message("Saved ", output_file, " - ", length(all_products), " products x ",
        length(years), " years, model = ", MODEL, ".")
