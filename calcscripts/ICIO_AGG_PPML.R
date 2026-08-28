## =====================================================================
##  ICIO_AGG_PPML.R                                       (script 2 of 10)
##  Yearly PPML (Poisson pseudo-maximum-likelihood) gravity regressions on
##  the AGGREGATED OECD ICIO dyadic data. One regression per year.
##  Dependent variable = Vol in LEVELS (not logged); the continuous
##  regressors (distance, GDP_i, GDP_j) are log-transformed.
##
##  Result: a single list 'Y_2D_GRAV_TIME' with
##    $beta    : data.frame [years x (12 coefficients + adjreg)]
##    $p_value : data.frame [years x (12 p-values  + model p-value)]
##  saved to '2D_GRAV_TIME_1_Y_PPML.RData'.
## =====================================================================

## ---- USER CONFIGURATION (adjust the paths to your environment) -------
MODEL       <- "PPML"
input_file  <- "DYADIC_ICIO_AGG_1995_2020.RData"   # input .RData
input_var   <- "DYADIC_ICIO_AGG"                   # object name inside it
output_file <- "2D_GRAV_TIME_1_Y_PPML.RData"       # output .RData

## ---- Load the shared estimation engine -------------------------------
.here <- "gravity_engine.R"
if (!file.exists(.here)) .here <- file.path("calcscripts", "gravity_engine.R")
if (!file.exists(.here)) stop("gravity_engine.R not found; keep it next to this script.")
source(.here)

## ---- Load data -------------------------------------------------------
DAT <- as.data.frame(load_one(input_file, input_var), stringsAsFactors = FALSE)

## ---- Dimensions derived from the data (robust to missing years) ------
years <- sort(unique(as.character(DAT$Year)))

beta    <- data.frame(matrix(NA_real_, length(years), length(OUT_COLS)))
p_value <- data.frame(matrix(NA_real_, length(years), length(OUT_COLS)))
rownames(beta) <- years; colnames(beta) <- OUT_COLS
rownames(p_value) <- years; colnames(p_value) <- OUT_COLS

## ---- One regression per year ----------------------------------------
idx <- split(seq_len(nrow(DAT)), factor(as.character(DAT$Year), levels = years))
for (y in years) {
  rows <- idx[[y]]
  if (length(rows) == 0L) next
  r <- estimate_cell(DAT[rows, , drop = FALSE], MODEL)
  beta[y, ]    <- cell_beta_vec(r)
  p_value[y, ] <- cell_pval_vec(r)
  message("Year ", y, " done (n = ", length(rows), ").")
}

## ---- Save ------------------------------------------------------------
Y_2D_GRAV_TIME <- list(beta = beta, p_value = p_value)
save(Y_2D_GRAV_TIME, file = output_file)
message("Saved ", output_file, " - ", length(years), " years, model = ", MODEL, ".")
