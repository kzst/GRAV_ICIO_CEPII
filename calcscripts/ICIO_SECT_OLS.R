## =====================================================================
##  ICIO_SECT_OLS.R                                       (script 3 of 10)
##  Log-log OLS gravity regressions on the SECTORAL OECD ICIO dyadic data.
##  One regression per (sector, year).
##
##  Result: a single list 'grav_sect_time' with $beta and $p_value, each a
##  3D array with dimensions:
##    [[1]] coefficients (12) + "adjreg"   (fit / model p-value)
##    [[2]] sectors
##    [[3]] years
##  saved to '3D_grav_sect_time_1_Y_OLS.RData'.
## =====================================================================

## ---- USER CONFIGURATION (adjust the paths to your environment) -------
MODEL       <- "OLS"
input_file  <- "DYADIC_ICIO_SECT_1995_2020.RData"
input_var   <- "DYADIC_ICIO_SECT"
output_file <- "3D_grav_sect_time_1_Y_OLS.RData"

## ---- Load the shared estimation engine -------------------------------
.here <- "gravity_engine.R"
if (!file.exists(.here)) .here <- file.path("SCRIPTS", "gravity_engine.R")
if (!file.exists(.here)) stop("gravity_engine.R not found; keep it next to this script.")
source(.here)

## ---- Load data -------------------------------------------------------
DAT <- as.data.frame(load_one(input_file, input_var), stringsAsFactors = FALSE)

## ---- Dimensions derived from the data --------------------------------
years   <- sort(unique(as.character(DAT$Year)))
sectors <- sort(unique(as.character(DAT$sector)))

dims  <- c(length(OUT_COLS), length(sectors), length(years))
dn    <- list(OUT_COLS, sectors, years)
beta    <- array(NA_real_, dim = dims, dimnames = dn)
p_value <- array(NA_real_, dim = dims, dimnames = dn)

## ---- One regression per (sector, year) ------------------------------
idx_year <- split(seq_len(nrow(DAT)), factor(as.character(DAT$Year), levels = years))
for (y in years) {
  ry <- idx_year[[y]]
  if (length(ry) == 0L) next
  dy <- DAT[ry, , drop = FALSE]
  idx_s <- split(seq_len(nrow(dy)), factor(as.character(dy$sector), levels = sectors))
  for (s in sectors) {
    rs <- idx_s[[s]]
    if (length(rs) == 0L) next
    r <- estimate_cell(dy[rs, , drop = FALSE], MODEL)
    beta[, s, y]    <- cell_beta_vec(r)
    p_value[, s, y] <- cell_pval_vec(r)
  }
  message("Year ", y, " done (", length(sectors), " sectors).")
}

## ---- Save ------------------------------------------------------------
grav_sect_time <- list(beta = beta, p_value = p_value)
save(grav_sect_time, file = output_file)
message("Saved ", output_file, " - ", length(sectors), " sectors x ",
        length(years), " years, model = ", MODEL, ".")
