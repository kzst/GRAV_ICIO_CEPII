## =====================================================================
##  CEPII_PG_OLS.R                                        (script 5 of 10)
##  Log-log OLS gravity regressions on the CEPII BACI data aggregated to
##  2-digit product groups. One regression per (product, year).
##
##  Origin / destination countries are stored in 'sect_i' / 'sect_j'
##  (they are not used as regressors); the product grouping is in 'product'.
##
##  Results: two 3D arrays saved to '3D_grav_products_time_20_OLS.RData':
##    var_products_time_20   : beta values
##    var_products_time_20_p : p-values
##  Dimensions of each array:
##    [[1]] coefficients (12) + "adjreg"   (fit / model p-value)
##    [[2]] products (2-digit)
##    [[3]] years
## =====================================================================

## ---- USER CONFIGURATION (adjust the paths to your environment) -------
MODEL       <- "OLS"
input_file  <- "DYADIC_BACI_2DIGIT_1995_2023.RData"
input_var   <- "DYADIC_BACI_2DIGIT"
output_file <- "3D_grav_products_time_20_OLS.RData"

## ---- Load the shared estimation engine -------------------------------
.here <- "gravity_engine.R"
if (!file.exists(.here)) .here <- file.path("SCRIPTS", "gravity_engine.R")
if (!file.exists(.here)) stop("gravity_engine.R not found; keep it next to this script.")
source(.here)

## ---- Load data -------------------------------------------------------
DAT <- as.data.frame(load_one(input_file, input_var), stringsAsFactors = FALSE)

## ---- Dimensions derived from the data --------------------------------
years    <- sort(unique(as.character(DAT$Year)))
products <- sort(unique(as.character(DAT$product)))

dims <- c(length(OUT_COLS), length(products), length(years))
dn   <- list(OUT_COLS, products, years)
var_products_time_20   <- array(NA_real_, dim = dims, dimnames = dn)
var_products_time_20_p <- array(NA_real_, dim = dims, dimnames = dn)

## ---- One regression per (product, year) -----------------------------
idx_year <- split(seq_len(nrow(DAT)), factor(as.character(DAT$Year), levels = years))
for (y in years) {
  ry <- idx_year[[y]]
  if (length(ry) == 0L) next
  dy <- DAT[ry, , drop = FALSE]
  idx_p <- split(seq_len(nrow(dy)), factor(as.character(dy$product), levels = products))
  for (p in products) {
    rp <- idx_p[[p]]
    if (length(rp) == 0L) next
    r <- estimate_cell(dy[rp, , drop = FALSE], MODEL)
    var_products_time_20[, p, y]   <- cell_beta_vec(r)
    var_products_time_20_p[, p, y] <- cell_pval_vec(r)
  }
  message("Year ", y, " done (", length(products), " products).")
}

## ---- Save ------------------------------------------------------------
save(var_products_time_20, var_products_time_20_p, file = output_file)
message("Saved ", output_file, " - ", length(products), " products x ",
        length(years), " years, model = ", MODEL, ".")
