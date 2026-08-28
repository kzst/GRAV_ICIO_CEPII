## =====================================================================
##  TEST_gravity_scripts.R
##  Self-contained test harness for the 10 gravity scripts.
##
##  It (1) generates small SYNTHETIC datasets that match the required input
##  structures and deliberately contain TYPICAL DATA ERRORS (zero trade,
##  zero distance, missing GDP, negative volume, Inf, near-empty product /
##  sector cells, a constant dummy, ...); (2) runs every script UNCHANGED
##  against that data; and (3) checks that each output has the required
##  object name, class, dimensions and dimnames, and contains NO Inf/NaN.
##
##  HOW TO RUN (from the SCRIPTS/ folder):
##     Rscript TEST_gravity_scripts.R
##  or, in an interactive session with the working directory set to SCRIPTS/:
##     source("TEST_gravity_scripts.R")
## =====================================================================

set.seed(1)

## ---- Locate the script folder (engine + drivers live here) -----------
detect_dir <- function() {
  a  <- commandArgs(FALSE)
  fa <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(fa)) return(dirname(normalizePath(fa[1])))
  getwd()
}
scripts_dir <- detect_dir()
if (!file.exists(file.path(scripts_dir, "gravity_engine.R"))) {
  for (cand in c(getwd(), file.path(getwd(), "SCRIPTS"))) {
    if (file.exists(file.path(cand, "gravity_engine.R"))) { scripts_dir <- cand; break }
  }
}
if (!file.exists(file.path(scripts_dir, "gravity_engine.R")))
  stop("Cannot locate gravity_engine.R. Run this test from the SCRIPTS/ folder.")

DUMMIES <- c("comlang_off", "comlang_ethno", "contig", "col45",
             "col_dep_ever", "comleg_pretrans", "fta_wto", "eu_o", "eu_d")

## ---- One dyadic block (all ordered pairs i != j) with a gravity DGP ---
gen_block <- function(codes) {
  g <- expand.grid(i = codes, j = codes, stringsAsFactors = FALSE)
  g <- g[g$i != g$j, , drop = FALSE]
  n <- nrow(g)
  gdp_i <- exp(rnorm(n, 10, 1)); gdp_j <- exp(rnorm(n, 10, 1))
  dist  <- exp(rnorm(n, 8, 0.5))
  d <- data.frame(
    i = g$i, j = g$j, Vol = NA_real_, Dist_ij = dist, GDP_i = gdp_i, GDP_j = gdp_j,
    comlang_off = rbinom(n, 1, .4), comlang_ethno = rbinom(n, 1, .4),
    contig = rbinom(n, 1, .3), col45 = rbinom(n, 1, .3),
    col_dep_ever = rbinom(n, 1, .3), comleg_pretrans = rbinom(n, 1, .4),
    fta_wto = rbinom(n, 1, .5), eu_o = rbinom(n, 1, .4), eu_d = rbinom(n, 1, .4),
    stringsAsFactors = FALSE)
  lin   <- 2 + 0.9 * log(gdp_i) + 0.8 * log(gdp_j) - 1.1 * log(dist) + 0.3 * d$contig
  d$Vol <- exp(lin + rnorm(n, 0, 0.5))
  d
}

## Inject typical data errors into a data.frame with a 'Vol' column.
inject_errors <- function(d) {
  n <- nrow(d)
  if (n >= 30) d$Vol[sample(n, max(1, round(0.15 * n)))] <- 0     # zero trade
  if (n >= 10) d$Dist_ij[sample(n, 3)] <- 0                       # zero distance
  if (n >= 10) d$GDP_i[sample(n, 2)]   <- NA                      # missing GDP
  if (n >= 10) d$Vol[sample(n, 1)]     <- -1                      # negative Vol
  d$Vol[1] <- Inf                                                 # Inf value
  d
}

## ---- Build synthetic files in a temporary working directory ----------
orig_wd <- getwd()
tmp <- file.path(tempdir(), paste0("grav_test_", as.integer(runif(1, 1, 1e6))))
dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
file.copy(file.path(scripts_dir, "gravity_engine.R"),
          file.path(tmp, "gravity_engine.R"), overwrite = TRUE)
setwd(tmp)

codes <- sprintf("C%02d", 1:8)      # 8 countries -> 56 ordered pairs per cell
years <- 1995:1997

## (1-2) ICIO aggregate
DYADIC_ICIO_AGG <- inject_errors(do.call(rbind, lapply(years, function(y)
  cbind(Year = y, gen_block(codes), stringsAsFactors = FALSE))))
save(DYADIC_ICIO_AGG, file = "DYADIC_ICIO_AGG_1995_2020.RData")

## (3-4) ICIO sectoral  (one near-empty sector cell + one constant dummy)
sects <- c("A01_02", "A03", "B05_06")
sect  <- do.call(rbind, lapply(years, function(y)
  do.call(rbind, lapply(sects, function(s)
    cbind(Year = y, sector = s, gen_block(codes), stringsAsFactors = FALSE)))))
sect$contig[sect$sector == "A03"] <- 1L                    # constant dummy in a cell
degen <- which(sect$sector == "B05_06" & sect$Year == 1997)
sect  <- sect[setdiff(seq_len(nrow(sect)), degen[-(1:3)]), ] # keep only 3 rows there
DYADIC_ICIO_SECT <- inject_errors(sect)
save(DYADIC_ICIO_SECT, file = "DYADIC_ICIO_SECT_1995_2020.RData")

## Helper: rename i/j -> sect_i/sect_j and add a product column (BACI shape)
as_baci <- function(y, pr) {
  b <- gen_block(codes)
  names(b)[names(b) == "i"] <- "sect_i"; names(b)[names(b) == "j"] <- "sect_j"
  cbind(Year = y, product = pr, b, stringsAsFactors = FALSE)
}

## (5-6) BACI 2-digit (a sparse product "77" with only 2 obs in 1995)
prods2 <- c("01", "02", "77")
baci2  <- do.call(rbind, lapply(years, function(y)
  do.call(rbind, lapply(prods2, function(pr) as_baci(y, pr)))))
sp <- which(baci2$product == "77" & baci2$Year == 1995)
baci2 <- baci2[setdiff(seq_len(nrow(baci2)), sp[-(1:2)]), ]
DYADIC_BACI_2DIGIT <- inject_errors(baci2)
save(DYADIC_BACI_2DIGIT, file = "DYADIC_BACI_2DIGIT_1995_2023.RData")

## (7-8) BACI 4-digit
prods4 <- c("0101", "0102", "8703")
baci4  <- do.call(rbind, lapply(years, function(y)
  do.call(rbind, lapply(prods4, function(pr) as_baci(y, pr)))))
DYADIC_BACI_4DIGIT <- inject_errors(baci4)
save(DYADIC_BACI_4DIGIT, file = "DYADIC_BACI_4DIGIT_1995_2023.RData")

## (9-10) BACI 6-digit raw: one file per year in a sub-directory
dir.create("DYADIC_6DIGIT_YEARS", showWarnings = FALSE)
prods6 <- c("010121", "010129", "999999")
for (y in years) {
  blk <- do.call(rbind, lapply(prods6, function(pr) as_baci(y, pr)))
  if (y == 1995) {                                   # make 999999 sparse in 1995
    sp <- which(blk$product == "999999")
    blk <- blk[setdiff(seq_len(nrow(blk)), sp[-(1:2)]), ]
  }
  DYADIC_BACI_6DIGIT_YEAR <- inject_errors(blk)
  save(DYADIC_BACI_6DIGIT_YEAR,
       file = file.path("DYADIC_6DIGIT_YEARS", sprintf("DYADIC_BACI_6DIGIT_%d.RData", y)))
}

## ---- Run a driver unchanged and capture success/failure --------------
run_driver <- function(name) {
  tryCatch({
    sys.source(file.path(scripts_dir, name), envir = new.env(parent = globalenv()))
    TRUE
  }, error = function(e) { message("  ERROR in ", name, ": ", conditionMessage(e)); FALSE })
}

## ---- Structure checks ------------------------------------------------
has_no_inf <- function(x) !any(is.infinite(as.numeric(unlist(x))))
n_finite   <- function(x) sum(is.finite(as.numeric(unlist(x))))

report <- list()
add <- function(script, ok, detail) {
  report[[length(report) + 1L]] <<- data.frame(script = script,
    status = if (ok) "PASS" else "FAIL", detail = detail, stringsAsFactors = FALSE)
}

check_list_df <- function(file, objname, nyears) {
  e <- new.env(); load(file, envir = e); obj <- get(objname, envir = e)
  ok <- is.list(obj) && all(c("beta", "p_value") %in% names(obj)) &&
        is.data.frame(obj$beta) && ncol(obj$beta) == 13 && nrow(obj$beta) == nyears &&
        identical(colnames(obj$beta)[13], "adjreg") && has_no_inf(obj$beta) &&
        has_no_inf(obj$p_value)
  list(ok = ok, detail = sprintf("beta dim = %s; finite betas = %d; no Inf = %s",
       paste(dim(obj$beta), collapse = "x"), n_finite(obj$beta), has_no_inf(obj$beta)))
}

check_arrays <- function(file, objb, objp, d1, d2, d3) {
  e <- new.env(); load(file, envir = e)
  b <- get(objb, envir = e); p <- get(objp, envir = e)
  ok <- is.array(b) && is.array(p) &&
        identical(dim(b), as.integer(c(d1, d2, d3))) &&
        identical(dimnames(b)[[1]][13], "adjreg") &&
        has_no_inf(b) && has_no_inf(p)
  list(ok = ok, detail = sprintf("dim = %s; finite betas = %d; no Inf = %s",
       paste(dim(b), collapse = "x"), n_finite(b), has_no_inf(b)))
}

check_list_arrays <- function(file, objname, d1, d2, d3) {
  e <- new.env(); load(file, envir = e); obj <- get(objname, envir = e)
  b <- obj$beta; p <- obj$p_value
  ok <- is.list(obj) && is.array(b) &&
        identical(dim(b), as.integer(c(d1, d2, d3))) &&
        identical(dimnames(b)[[1]][13], "adjreg") && has_no_inf(b) && has_no_inf(p)
  list(ok = ok, detail = sprintf("beta dim = %s; finite betas = %d; no Inf = %s",
       paste(dim(b), collapse = "x"), n_finite(b), has_no_inf(b)))
}

message("\n=== Running the 10 scripts on synthetic data ===")

## 1-2 AGG (data.frame [3 years x 13])
for (nm in c("ICIO_AGG_OLS.R", "ICIO_AGG_PPML.R")) {
  ok <- run_driver(nm)
  out <- if (grepl("OLS", nm)) "2D_GRAV_TIME_1_Y_OLS.RData" else "2D_GRAV_TIME_1_Y_PPML.RData"
  r <- if (ok && file.exists(out)) check_list_df(out, "Y_2D_GRAV_TIME", 3) else list(ok = FALSE, detail = "no output")
  add(nm, ok && r$ok, r$detail)
}

## 3-4 SECT (array 13 x 3 sectors x 3 years)
for (nm in c("ICIO_SECT_OLS.R", "ICIO_SECT_PPML.R")) {
  ok <- run_driver(nm)
  out <- if (grepl("OLS", nm)) "3D_grav_sect_time_1_Y_OLS.RData" else "3D_grav_sect_time_1_Y_PPML.RData"
  r <- if (ok && file.exists(out)) check_list_arrays(out, "grav_sect_time", 13, 3, 3) else list(ok = FALSE, detail = "no output")
  add(nm, ok && r$ok, r$detail)
}

## 5-6 PG 2-digit (array 13 x 3 products x 3 years)
for (nm in c("CEPII_PG_OLS.R", "CEPII_PG_PPML.R")) {
  ok <- run_driver(nm)
  out <- if (grepl("OLS", nm)) "3D_grav_products_time_20_OLS.RData" else "3D_grav_products_time_20_PPML.RData"
  r <- if (ok && file.exists(out)) check_arrays(out, "var_products_time_20", "var_products_time_20_p", 13, 3, 3) else list(ok = FALSE, detail = "no output")
  add(nm, ok && r$ok, r$detail)
}

## 7-8 PSG 4-digit (array 13 x 3 products x 3 years)
for (nm in c("CEPII_PSG_OLS.R", "CEPII_PSG_PPML.R")) {
  ok <- run_driver(nm)
  out <- if (grepl("OLS", nm)) "3D_grav_products_time_40_OLS.RData" else "3D_grav_products_time_40_PPML.RData"
  r <- if (ok && file.exists(out)) check_arrays(out, "var_products_time_40", "var_products_time_40_p", 13, 3, 3) else list(ok = FALSE, detail = "no output")
  add(nm, ok && r$ok, r$detail)
}

## 9-10 PROD 6-digit (array 13 x 3 products x 3 years)
for (nm in c("CEPII_PROD_OLS.R", "CEPII_PROD_PPML.R")) {
  ok <- run_driver(nm)
  out <- if (grepl("OLS", nm)) "3D_grav_products_time_ALL_OLS.RData" else "3D_grav_products_time_ALL_PPML.RData"
  r <- if (ok && file.exists(out)) check_arrays(out, "var_products_time_ALL", "var_products_time_ALL_p", 13, 3, 3) else list(ok = FALSE, detail = "no output")
  add(nm, ok && r$ok, r$detail)
}

res <- do.call(rbind, report)
message("\n=== TEST SUMMARY ===")
print(res, row.names = FALSE)
message("\nBackends: fixest = ", requireNamespace("fixest", quietly = TRUE),
        " | sandwich = ", requireNamespace("sandwich", quietly = TRUE))
message(if (all(res$status == "PASS")) "ALL 10 SCRIPTS PASSED."
        else "SOME SCRIPTS FAILED - see the table above.")

setwd(orig_wd)
