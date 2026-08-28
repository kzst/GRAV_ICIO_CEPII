## =====================================================================
##  OLS_PPML_COMPARISON.R                                (script 11 of 11)
##  Compares the OLS and PPML estimates of the 12 gravity-model slope
##  coefficients (distance, the two GDP elasticities, and the nine CEPII
##  bilateral controls) produced by the 10 paired calc scripts, for every
##  database, product/sector, and year. This script does NOT estimate
##  anything itself -- it only reads the .RData outputs already produced
##  by the 10 OLS/PPML scripts in this folder and builds the comparison
##  table(s) needed for Appendix \ref{app:ols_ppml_comparison} of
##  GRAV_ICIO_CEPII_PAPER.Rmd.md ("Comparison of OLS and PPML Estimates
##  for the Bilateral Control Variables").
##
##  RUN THE 10 CALC SCRIPTS FIRST. This script looks for some or all of
##  the following 10 files (produced by ICIO_AGG_OLS.R / _PPML.R,
##  ICIO_SECT_OLS.R / _PPML.R, CEPII_PG_OLS.R / _PPML.R, CEPII_PSG_OLS.R /
##  _PPML.R, CEPII_PROD_OLS.R / _PPML.R):
##    2D_GRAV_TIME_1_Y_OLS.RData            2D_GRAV_TIME_1_Y_PPML.RData
##    3D_grav_sect_time_1_Y_OLS.RData       3D_grav_sect_time_1_Y_PPML.RData
##    3D_grav_products_time_20_OLS.RData    3D_grav_products_time_20_PPML.RData
##    3D_grav_products_time_40_OLS.RData    3D_grav_products_time_40_PPML.RData
##    3D_grav_products_time_ALL_OLS.RData   3D_grav_products_time_ALL_PPML.RData
##  Missing pairs are SKIPPED with a message (not an error), so this can
##  already be run once even a single database pair is ready; rerun it
##  later once more pairs finish -- the tables simply become more complete.
##  These 10 files are expected in the 'data/' subdirectory (a sibling of
##  this calcscripts/ folder) -- see 'search_dirs' below if yours differ.
##
##  OUTPUT (written to 'out_dir', default "output/", i.e. R1/output/):
##    TABLE_S_OLS_PPML_COMPARISON.xlsx    - the same comparison, as an
##                                           Excel workbook (ReadMe +
##                                           Pooled + By_Database sheets,
##                                           plus the full row-level panel
##                                           when small enough for one
##                                           Excel sheet)
##    ols_ppml_comparison_long.csv        - tidy panel, one row per
##                                           (database, product/sector,
##                                           year, coefficient); this is
##                                           the raw material everything
##                                           else below is computed from
##    ols_ppml_comparison_by_database.csv - metrics x (database, coef.)
##    ols_ppml_comparison_pooled.csv      - metrics x coefficient, pooled
##                                           across all loaded databases
##    ols_ppml_comparison_table.tex       - two ready-to-paste booktabs
##                                           tables (core elasticities by
##                                           database; control-variable
##                                           sign stability, pooled)
##
##  METRICS reported per (database x coefficient) and pooled:
##    n                   - number of comparable, non-missing cells
##    correlation         - Pearson correlation of the OLS vs. PPML point
##                           estimates across cells
##    mean_abs_diff / median_abs_diff / rmse - level of the OLS-PPML gap
##    sign_agreement_pct  - % of cells where OLS and PPML agree in sign
##    both_sig_pct        - % of cells significant (p < alpha) in BOTH
##    sig_disagree_pct    - % of cells significant in exactly one of the two
##
##  HOW TO RUN (from the calcscripts/ folder, once outputs exist):
##     Rscript OLS_PPML_COMPARISON.R
##  or, interactively with the working directory set to calcscripts/ (or
##  wherever the 10 output files were saved -- see 'search_dirs' below):
##     source("OLS_PPML_COMPARISON.R")
##  In THIS project, run it from the paper's Rmd (working directory R1/,
##  data expected in R1/data/) as:
##     source("calcscripts/OLS_PPML_COMPARISON.r")
##
##  Dependencies: NONE beyond base R (no fixest/sandwich needed here --
##  this script only reshapes and summarizes already-estimated results).
## =====================================================================

## ---- USER CONFIGURATION (adjust to your environment) -------------------
## Folders to search for the 10 OLS/PPML .RData outputs, tried in order;
## the first folder that contains a given file wins. Primary location is
## 'data/', matching the actual project layout (R1/data/ holds the 10
## .RData files produced by the calc scripts, under exactly the names
## those scripts save them as; R1/calcscripts/ holds this script itself).
## The two lines below also try to resolve 'data/' as a sibling of this
## script's OWN folder, so the search still works if this script is
## source()'d from a working directory other than R1/ (e.g. run directly
## from inside calcscripts/); this is a pure best-effort addition and
## silently yields NA (skipped) if the running context does not expose
## it, in which case only the cwd-relative candidates below are used.
.get_self_path <- function() {
  for (fr in sys.frames()) {
    ofile <- fr$ofile
    if (!is.null(ofile)) return(normalizePath(ofile, mustWork = FALSE))
  }
  NA_character_
}
.self_path    <- tryCatch(.get_self_path(), error = function(e) NA_character_)
.sibling_data <- if (!is.na(.self_path)) file.path(dirname(.self_path), "..", "data") else NA_character_

search_dirs <- c("data", .sibling_data, ".", "output", file.path("..", "output"), "..")
search_dirs <- unique(search_dirs[!is.na(search_dirs)])

## Folder where this script's own outputs are written (created if
## needed); "output" resolves to R1/output/ when this script is sourced
## from the Rmd sitting in R1/ (source("calcscripts/OLS_PPML_COMPARISON.r")),
## matching where every other calc script in this project already writes
## its own tables and figures.
out_dir <- "output"

## Set to FALSE to skip writing the (potentially large, millions-of-rows)
## row-level long CSV for the 6-digit CEPII_PROD family. The summary
## tables (by-database and pooled) are written either way.
write_long_csv <- TRUE

## Significance threshold used for the *_sig_pct metrics.
alpha <- 0.05

## =====================================================================
##  Nothing below this line needs to be edited.
## =====================================================================

## ---- Coefficient metadata (MUST match gravity_engine.R exactly) -------
COEF_NAMES <- c("dist", "GDPi", "GDPj", "comlang_off", "comlang_ethno",
                "contig", "col45", "col_dep_ever", "comleg_pretrans",
                "fta_wto", "eu_o", "eu_d")
CORE_COEFS <- c("dist", "GDPi", "GDPj")
CTRL_COEFS <- setdiff(COEF_NAMES, CORE_COEFS)
COEF_LABEL <- c(
  dist            = "Distance",
  GDPi            = "Origin GDP",
  GDPj            = "Destination GDP",
  comlang_off     = "Common official language",
  comlang_ethno   = "Common ethnic language",
  contig          = "Contiguity",
  col45           = "Colonial relationship post-1945",
  col_dep_ever    = "Ever-colonial relationship",
  comleg_pretrans = "Common legal origin (pre-transition)",
  fta_wto         = "Joint FTA/WTO membership",
  eu_o            = "EU membership (origin)",
  eu_d            = "EU membership (destination)")

## ---- Database family registry ------------------------------------------
## pattern "df"     : list(beta = data.frame[years x 13], p_value = ...)
##                    (ICIO_AGG; no product/sector dimension)
## pattern "list3d" : list(beta = array[13, group, years], p_value = ...)
##                    (ICIO_SECT)
## pattern "obj3d"  : two bare top-level arrays  <name>  and  <name>_p
##                    (CEPII_PG / CEPII_PSG / CEPII_PROD)
families <- list(
  list(db = "ICIO_AGG",  pattern = "df",
       ols_file  = "2D_GRAV_TIME_1_Y_OLS.RData",  ols_obj  = "Y_2D_GRAV_TIME",
       ppml_file = "2D_GRAV_TIME_1_Y_PPML.RData", ppml_obj = "Y_2D_GRAV_TIME",
       group_label = NA_character_),
  list(db = "ICIO_SECT", pattern = "list3d",
       ols_file  = "3D_grav_sect_time_1_Y_OLS.RData",  ols_obj  = "grav_sect_time",
       ppml_file = "3D_grav_sect_time_1_Y_PPML.RData", ppml_obj = "grav_sect_time",
       group_label = "sector"),
  list(db = "CEPII_PG",  pattern = "obj3d",
       ols_file  = "3D_grav_products_time_20_OLS.RData",  ols_obj  = "var_products_time_20",
       ppml_file = "3D_grav_products_time_20_PPML.RData", ppml_obj = "var_products_time_20",
       group_label = "product"),
  list(db = "CEPII_PSG", pattern = "obj3d",
       ols_file  = "3D_grav_products_time_40_OLS.RData",  ols_obj  = "var_products_time_40",
       ppml_file = "3D_grav_products_time_40_PPML.RData", ppml_obj = "var_products_time_40",
       group_label = "product"),
  list(db = "CEPII_PROD", pattern = "obj3d",
       ols_file  = "3D_grav_products_time_ALL_OLS.RData",  ols_obj  = "var_products_time_ALL",
       ppml_file = "3D_grav_products_time_ALL_PPML.RData", ppml_obj = "var_products_time_ALL",
       group_label = "product")
)
db_order <- vapply(families, function(f) f$db, character(1))

## ---- Small helpers -------------------------------------------------------
find_file <- function(fname, dirs) {
  for (d in dirs) {
    p <- file.path(d, fname)
    if (file.exists(p)) return(p)
  }
  NA_character_
}

## Load a named object from an .RData file into a fresh environment.
load_named <- function(path, objname) {
  e <- new.env()
  load(path, envir = e)
  if (!objname %in% ls(e))
    stop("Object '", objname, "' not found in ", path,
         " (found: ", paste(ls(e), collapse = ", "), ")")
  get(objname, envir = e)
}

## The p-value companion object for the "obj3d" pattern always has this name.
p_obj_name <- function(objname) paste0(objname, "_p")

## Replace any non-finite numeric value with NA (never let Inf/NaN leak
## into the CSV / LaTeX outputs), mirroring gravity_engine.R's own
## na_if_nonfinite() philosophy.
sanitize_nonfinite <- function(df, skip = c("database", "coefficient", "coef_label", "group_type")) {
  for (col in setdiff(names(df), skip)) {
    x <- df[[col]]
    if (is.numeric(x)) { x[!is.finite(x)] <- NA_real_; df[[col]] <- x }
  }
  df
}

## ---- Tidy-panel builders: one row per (group, year, coefficient) -------
## "df" pattern (ICIO_AGG): beta / p_value are data.frames [years x 13].
tidy_df_pattern <- function(fam, ols, ppml) {
  years <- intersect(rownames(ols$beta), rownames(ppml$beta))
  if (length(years) == 0L) {
    warning(fam$db, ": OLS and PPML share no common years; skipped.")
    return(NULL)
  }
  do.call(rbind, lapply(COEF_NAMES, function(cf) {
    data.frame(database = fam$db, group_type = NA_character_, group = NA_character_,
               year = years, coefficient = cf,
               beta_ols  = ols$beta[years, cf],    beta_ppml = ppml$beta[years, cf],
               p_ols     = ols$p_value[years, cf], p_ppml    = ppml$p_value[years, cf],
               stringsAsFactors = FALSE)
  }))
}

## "list3d" / "obj3d" patterns: beta / p_value (or the bare arrays) are 3D
## arrays [13, group, years] with matching dimnames. Both arrays are
## indexed by NAME using the SAME (groups, years) vectors, so the two
## results are always aligned even if they store dimensions in a
## different internal order; as.vector() on the resulting [group x year]
## matrix flattens it column-major, matching expand.grid(group=groups,
## year=years) row order (group varies fastest) -- this indexing +
## flattening + alignment logic was cross-checked against a hand-built
## synthetic case (with a known, deliberately planted sign mismatch)
## before this script was finalized.
tidy_array_pattern <- function(fam, ols_beta, ols_p, ppml_beta, ppml_p) {
  groups <- intersect(dimnames(ols_beta)[[2]], dimnames(ppml_beta)[[2]])
  years  <- intersect(dimnames(ols_beta)[[3]], dimnames(ppml_beta)[[3]])
  if (length(groups) == 0L || length(years) == 0L) {
    warning(fam$db, ": OLS and PPML share no common product/sector-year cells; skipped.")
    return(NULL)
  }
  grid <- expand.grid(group = groups, year = years, stringsAsFactors = FALSE)
  do.call(rbind, lapply(COEF_NAMES, function(cf) {
    data.frame(database = fam$db, group_type = fam$group_label,
               group = grid$group, year = grid$year, coefficient = cf,
               beta_ols  = as.vector(ols_beta[cf, groups, years]),
               beta_ppml = as.vector(ppml_beta[cf, groups, years]),
               p_ols     = as.vector(ols_p[cf, groups, years]),
               p_ppml    = as.vector(ppml_p[cf, groups, years]),
               stringsAsFactors = FALSE)
  }))
}

## ---- Load every available family and stack into one tidy panel --------
message("=== OLS_PPML_COMPARISON.R: loading available outputs ===")
all_rows <- list()
missing  <- character(0)

for (fam in families) {
  ols_path  <- find_file(fam$ols_file,  search_dirs)
  ppml_path <- find_file(fam$ppml_file, search_dirs)
  if (is.na(ols_path) || is.na(ppml_path)) {
    missing <- c(missing, fam$db)
    message("  [SKIP] ", fam$db, ": missing ",
            if (is.na(ols_path)) fam$ols_file else fam$ppml_file,
            " in ", paste(search_dirs, collapse = ", "))
    next
  }
  message("  [OK]   ", fam$db, ": ", fam$ols_file, " + ", fam$ppml_file)

  rows <- tryCatch({
    if (fam$pattern == "df") {
      ols  <- load_named(ols_path,  fam$ols_obj)
      ppml <- load_named(ppml_path, fam$ppml_obj)
      tidy_df_pattern(fam, ols, ppml)
    } else if (fam$pattern == "list3d") {
      ols  <- load_named(ols_path,  fam$ols_obj)
      ppml <- load_named(ppml_path, fam$ppml_obj)
      tidy_array_pattern(fam, ols$beta, ols$p_value, ppml$beta, ppml$p_value)
    } else { # "obj3d"
      ols_beta  <- load_named(ols_path,  fam$ols_obj)
      ols_p     <- load_named(ols_path,  p_obj_name(fam$ols_obj))
      ppml_beta <- load_named(ppml_path, fam$ppml_obj)
      ppml_p    <- load_named(ppml_path, p_obj_name(fam$ppml_obj))
      tidy_array_pattern(fam, ols_beta, ols_p, ppml_beta, ppml_p)
    }
  }, error = function(e) {
    message("  [ERROR] ", fam$db, ": ", conditionMessage(e)); NULL
  })
  if (!is.null(rows)) all_rows[[fam$db]] <- rows
}

if (length(all_rows) == 0L)
  stop("None of the 10 OLS/PPML output files were found in: ",
       paste(search_dirs, collapse = ", "),
       ". Run the 10 calc scripts first, then either copy their outputs ",
       "next to this script or edit 'search_dirs' above to point at them.")

long <- do.call(rbind, all_rows)
rownames(long) <- NULL
long$coef_label <- COEF_LABEL[long$coefficient]

if (length(missing))
  message("\nNOTE: ", length(missing), " of ", length(families),
          " databases have no output yet and were skipped: ",
          paste(missing, collapse = ", "),
          ". The tables below cover the remaining ",
          length(families) - length(missing), " database(s) only -- ",
          "rerun this script after the missing calc scripts finish.")

## ---- Derived per-cell comparison fields --------------------------------
long$diff       <- long$beta_ppml - long$beta_ols
sign_ols        <- sign(long$beta_ols)
sign_ppml       <- sign(long$beta_ppml)
long$sign_agree <- ifelse(sign_ols == 0 | sign_ppml == 0, NA, sign_ols == sign_ppml)
long$both_sig     <- long$p_ols < alpha & long$p_ppml < alpha
long$sig_disagree <- (long$p_ols < alpha) != (long$p_ppml < alpha)

## ---- Summary metrics for one subset of rows ----------------------------
summarize_pair <- function(d) {
  ok <- is.finite(d$beta_ols) & is.finite(d$beta_ppml)
  n  <- sum(ok)
  if (n < 2L) {
    return(data.frame(n = n, correlation = NA_real_, mean_abs_diff = NA_real_,
                       median_abs_diff = NA_real_, rmse = NA_real_,
                       sign_agreement_pct = NA_real_, both_sig_pct = NA_real_,
                       sig_disagree_pct = NA_real_))
  }
  d  <- d[ok, , drop = FALSE]
  sa <- d$sign_agree[!is.na(d$sign_agree)]
  cc <- suppressWarnings(stats::cor(d$beta_ols, d$beta_ppml))
  data.frame(
    n                  = n,
    correlation        = if (is.finite(cc)) cc else NA_real_,
    mean_abs_diff       = mean(abs(d$diff)),
    median_abs_diff     = stats::median(abs(d$diff)),
    rmse                = sqrt(mean(d$diff^2)),
    sign_agreement_pct  = if (length(sa)) 100 * mean(sa) else NA_real_,
    both_sig_pct        = 100 * mean(d$both_sig,     na.rm = TRUE),
    sig_disagree_pct    = 100 * mean(d$sig_disagree, na.rm = TRUE)
  )
}

## By database x coefficient.
groups_db_coef <- split(long, list(long$database, long$coefficient), drop = TRUE)
by_db_coef <- do.call(rbind, lapply(groups_db_coef, function(d) {
  cbind(database = d$database[1], coefficient = d$coefficient[1],
        coef_label = d$coef_label[1], summarize_pair(d))
}))
rownames(by_db_coef) <- NULL
by_db_coef <- by_db_coef[order(match(by_db_coef$coefficient, COEF_NAMES),
                                match(by_db_coef$database, db_order)), ]

## Pooled across all loaded databases, by coefficient only.
groups_coef <- split(long, long$coefficient, drop = TRUE)
pooled_coef <- do.call(rbind, lapply(groups_coef, function(d) {
  cbind(coefficient = d$coefficient[1], coef_label = d$coef_label[1], summarize_pair(d))
}))
rownames(pooled_coef) <- NULL
pooled_coef <- pooled_coef[match(COEF_NAMES, pooled_coef$coefficient), ]

by_db_coef  <- sanitize_nonfinite(by_db_coef)
pooled_coef <- sanitize_nonfinite(pooled_coef)

## ---- Sanity checks (fail loudly on structural bugs, not silently) -----
in_pct_range <- function(x) x >= 0 & x <= 100
stopifnot(
  all(COEF_NAMES %in% unique(long$coefficient)),
  all(is.na(by_db_coef$sign_agreement_pct)  | in_pct_range(by_db_coef$sign_agreement_pct)),
  all(is.na(pooled_coef$sign_agreement_pct) | in_pct_range(pooled_coef$sign_agreement_pct))
)
message("\nSanity checks passed.")

## ---- Write CSV outputs --------------------------------------------------
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

if (write_long_csv)
  write.csv(long, file.path(out_dir, "ols_ppml_comparison_long.csv"), row.names = FALSE)
write.csv(by_db_coef,  file.path(out_dir, "ols_ppml_comparison_by_database.csv"), row.names = FALSE)
write.csv(pooled_coef, file.path(out_dir, "ols_ppml_comparison_pooled.csv"),      row.names = FALSE)

message("\n=== Wrote ===")
if (write_long_csv)
  message("  ", file.path(out_dir, "ols_ppml_comparison_long.csv"),      " (", nrow(long), " rows)")
message("  ", file.path(out_dir, "ols_ppml_comparison_by_database.csv"), " (", nrow(by_db_coef), " rows)")
message("  ", file.path(out_dir, "ols_ppml_comparison_pooled.csv"),      " (", nrow(pooled_coef), " rows)")

## ---- Excel workbook (preferred human-readable comparison table) -------
xlsx_path <- file.path(out_dir, "TABLE_S_OLS_PPML_COMPARISON.xlsx")
EXCEL_ROW_LIMIT <- 1000000L  # safety margin under Excel's 1,048,576 hard limit
xlsx_ok <- FALSE

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  message("\n[WARNING] Package 'openxlsx' is not installed -- ", xlsx_path,
          " was NOT written.\n          Install it with install.packages(\"openxlsx\") ",
          "and re-run this script;\n          the CSV files above already contain the same numbers.")
} else {
  xlsx_ok <- tryCatch({
    wb <- openxlsx::createWorkbook()

    hdr_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D9E1F2")
    pct_style <- openxlsx::createStyle(numFmt = '0.0"%"')
    num_style <- openxlsx::createStyle(numFmt = "0.000")

    write_sheet <- function(df, sheet_name) {
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet_name, df, headerStyle = hdr_style)
      pct_cols <- grep("_pct$", names(df))
      num_cols <- match(c("correlation", "mean_abs_diff", "median_abs_diff", "rmse"), names(df))
      num_cols <- num_cols[!is.na(num_cols)]
      if (nrow(df) > 0) {
        if (length(pct_cols))
          openxlsx::addStyle(wb, sheet_name, pct_style, rows = 2:(nrow(df) + 1L),
                              cols = pct_cols, gridExpand = TRUE)
        if (length(num_cols))
          openxlsx::addStyle(wb, sheet_name, num_style, rows = 2:(nrow(df) + 1L),
                              cols = num_cols, gridExpand = TRUE)
      }
      tryCatch(openxlsx::freezePane(wb, sheet_name, firstRow = TRUE), error = function(e) NULL)
      tryCatch(openxlsx::setColWidths(wb, sheet_name, cols = seq_len(ncol(df)), widths = "auto"),
                error = function(e) NULL)
    }

    status_df <- data.frame(
      database = db_order,
      status   = ifelse(db_order %in% missing, "Missing (calc script not yet run)", "Loaded"),
      stringsAsFactors = FALSE)

    readme <- data.frame(
      Field = c("Generated", "",
                "n", "correlation", "mean_abs_diff / median_abs_diff / rmse",
                "sign_agreement_pct", "both_sig_pct", "sig_disagree_pct",
                "", "Database status"),
      Description = c(
        as.character(Sys.time()), "",
        "Number of comparable, non-missing (OLS, PPML) cells.",
        "Pearson correlation of the OLS vs. PPML point estimates across cells.",
        "Level of the OLS-PPML gap (mean / median absolute difference; root-mean-square error).",
        "Share of cells (%) where OLS and PPML point estimates agree in sign.",
        "Share of cells (%) significant at the 5% level under BOTH estimators.",
        "Share of cells (%) significant under exactly one of the two estimators.",
        "", "see table below"),
      stringsAsFactors = FALSE)

    openxlsx::addWorksheet(wb, "ReadMe")
    openxlsx::writeData(wb, "ReadMe", readme, headerStyle = hdr_style)
    openxlsx::writeData(wb, "ReadMe", status_df, startRow = nrow(readme) + 3, headerStyle = hdr_style)
    tryCatch(openxlsx::setColWidths(wb, "ReadMe", cols = 1:2, widths = c(22, 70)), error = function(e) NULL)

    write_sheet(pooled_coef, "Pooled")
    write_sheet(by_db_coef,  "By_Database")

    if (nrow(long) <= EXCEL_ROW_LIMIT) {
      write_sheet(long, "Long_Panel")
      long_note <- paste0(nrow(long), " rows included")
    } else {
      long_note <- paste0(nrow(long), " rows -- too large for one Excel sheet; ",
                           "see ols_ppml_comparison_long.csv instead")
    }

    openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
    message("  ", xlsx_path, " (ReadMe + Pooled + By_Database",
            if (nrow(long) <= EXCEL_ROW_LIMIT) " + Long_Panel" else "", " sheets; ", long_note, ")")
    TRUE
  }, error = function(e) {
    message("\n[WARNING] Could not write ", xlsx_path, " (", conditionMessage(e), ").\n",
            "          The CSV files above already contain the same numbers.")
    FALSE
  })
}

## ---- LaTeX tables (booktabs / tabularx, matching the paper's style) ----
fmt_num <- function(x, digits = 3) ifelse(is.na(x), "--", formatC(x, format = "f", digits = digits))
fmt_pct <- function(x, digits = 1) ifelse(is.na(x), "--", paste0(formatC(x, format = "f", digits = digits), "\\%"))
esc_us  <- function(x) gsub("_", "\\_", x, fixed = TRUE)

tex <- c(
  "% Auto-generated by OLS_PPML_COMPARISON.R -- paste into the",
  "% 'Comparison of OLS and PPML Estimates for the Bilateral Control",
  "% Variables' appendix (label app:ols_ppml_comparison) and remove the",
  "% current editorial placeholder note there.",
  "\\begin{table}[htbp]",
  "\\centering",
  paste0("\\caption{Comparison of \\gls{ols} and \\gls{ppml} estimates of the core ",
         "distance and \\gls{gdp} elasticities, by database. $N$ is the number of ",
         "comparable (product/sector, year) cells; correlation and the mean absolute ",
         "difference are computed on those cells; sign agreement is the share of ",
         "cells where the \\gls{ols} and \\gls{ppml} point estimates have the same ",
         "sign; ``both sig.'' is the share of cells significant at the 5\\% level ",
         "under both estimators.}"),
  "\\label{tab:ols_ppml_core}",
  "\\begin{tabularx}{\\textwidth}{l l r r r r r}",
  "\\toprule",
  "Coefficient & Database & $N$ & Correlation & Mean $|\\Delta|$ & Sign agreem. & Both sig. \\\\",
  "\\midrule"
)

for (cf in CORE_COEFS) {
  sub <- by_db_coef[by_db_coef$coefficient == cf, ]
  sub <- sub[match(db_order, sub$database), ]
  sub <- sub[!is.na(sub$database), ]
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    tex <- c(tex, sprintf("%s & %s & %d & %s & %s & %s & %s \\\\",
      if (i == 1) COEF_LABEL[[cf]] else "", esc_us(r$database), r$n,
      fmt_num(r$correlation), fmt_num(r$mean_abs_diff),
      fmt_pct(r$sign_agreement_pct), fmt_pct(r$both_sig_pct)))
  }
  pr <- pooled_coef[pooled_coef$coefficient == cf, ]
  if (nrow(pr) == 1L)
    tex <- c(tex, sprintf("%s & \\emph{Pooled} & %d & %s & %s & %s & %s \\\\",
      "", pr$n, fmt_num(pr$correlation), fmt_num(pr$mean_abs_diff),
      fmt_pct(pr$sign_agreement_pct), fmt_pct(pr$both_sig_pct)))
  if (cf != CORE_COEFS[length(CORE_COEFS)]) tex <- c(tex, "\\addlinespace")
}
tex <- c(tex, "\\bottomrule", "\\end{tabularx}", "\\end{table}", "")

tex <- c(tex,
  "\\begin{table}[htbp]",
  "\\centering",
  paste0("\\caption{Sign stability of the nine CEPII bilateral control coefficients ",
         "($\\delta_{m,t}^{k}$) between the \\gls{ols} and \\gls{ppml} specifications, ",
         "pooled across all loaded databases, products/sectors, and years.}"),
  "\\label{tab:ols_ppml_controls}",
  "\\begin{tabularx}{\\textwidth}{X r r r}",
  "\\toprule",
  "Control variable & $N$ & Sign agreem. & Both sig. \\\\",
  "\\midrule"
)
for (cf in CTRL_COEFS) {
  r <- pooled_coef[pooled_coef$coefficient == cf, ]
  if (nrow(r) == 1L)
    tex <- c(tex, sprintf("%s & %d & %s & %s \\\\",
      COEF_LABEL[[cf]], r$n, fmt_pct(r$sign_agreement_pct), fmt_pct(r$both_sig_pct)))
}
tex <- c(tex, "\\bottomrule", "\\end{tabularx}", "\\end{table}")

writeLines(tex, file.path(out_dir, "ols_ppml_comparison_table.tex"))
message("  ", file.path(out_dir, "ols_ppml_comparison_table.tex"), " (2 booktabs tables)")

message("\n=== Pooled summary (coefficient level) ===")
print(pooled_coef[, c("coef_label", "n", "correlation", "mean_abs_diff",
                       "sign_agreement_pct", "both_sig_pct")], row.names = FALSE)

message("\nNext step: open ", file.path(out_dir, "ols_ppml_comparison_table.tex"),
        ", paste the two tables into the appendix section titled",
        "\n'Comparison of OLS and PPML Estimates for the Bilateral Control Variables'",
        " (label app:ols_ppml_comparison) in GRAV_ICIO_CEPII_PAPER.Rmd.md,",
        "\nand remove the current editorial placeholder note there.")

## ---- Bundle results into ONE object (kept out of the global namespace
## except for this single list) for convenient reuse if this script is
## source()'d from inside a larger session, e.g. the paper's Rmd itself.
OLS_PPML_COMPARISON <- list(
  long              = long,
  by_database       = by_db_coef,
  pooled            = pooled_coef,
  missing_databases = missing,
  xlsx_written      = xlsx_ok,
  xlsx_path         = xlsx_path,
  out_dir           = out_dir)
message("\nAll results also available in the list 'OLS_PPML_COMPARISON' ",
        "(elements: long, by_database, pooled, missing_databases, ",
        "xlsx_written, xlsx_path, out_dir).")
