# plot_granger_rf_shadow.R
#
# Three-way visualization companion to plot_granger_rf() (see plot_granger_rf.R),
# purpose-built for the output of granger_rf_shadow() (see granger_rf_shadow.R,
# which addresses The_World_Economy_Rev.md, Point 6).
#
# WHY A SEPARATE FUNCTION (not a modification of plot_granger_rf.R):
# Same rationale as granger_rf_shadow.R itself -- your existing, already-
# validated plot_granger_rf() stays untouched. This new function is meant to
# be dropped into the SAME slot in the "rfnetcaus" chunk once you switch the
# reported figure to the shadow-corrected fit.
#
# WHAT'S NEW COMPARED TO plot_granger_rf():
# plot_granger_rf() only knows a binary state (Causal / No causality, i.e.
# importance > 0 or not). This function draws THREE states, all derived from
# the SAME shadow-augmented random forest fit used inside granger_rf_shadow()
# (so the counts below are directly comparable to res$n_flagged_raw /
# res$n_flagged_shadow):
#
#   - "Causal (shadow-robust)"    granger_precedes_shadow == TRUE     -> green
#   - "Borderline (fails shadow)" granger_precedes_raw & !...shadow   -> amber
#   - "No causality"              granger_precedes_raw == FALSE      -> red
#
# The "amber" bars are exactly the variables that the OLD decision rule
# (granger_rf()'s "importance > drop_threshold", drop_threshold = 0 by
# default) would have reported as causal, but which do NOT survive the
# shadow-feature significance test. Plotting them in their ORIGINAL best-lag
# column (not collapsed to column 0) makes the correction visually legible:
# a reader can see exactly which bars "used to be green" and now are not.
#
# USAGE (drop-in replacement inside the "rfnetcaus" chunk):
#   fit_dist_rf_shadow <- granger_rf_shadow(Y_dist, X)
#   plot_granger_rf_shadow(fit_dist_rf_shadow, y_label = "Dist", cex_base = 0.7,
#                           sizes = list(var_rect_w = 0.95, var_rect_h = 0.70,
#                                        y_rect_w = 0.50, y_rect_h = 0.90,
#                                        hist_left = -1.20, hist_max_width = 0.95,
#                                        label_offset = 0.06))
#
# IMPORTANT: this script has not been executed against your actual data. The
# geometry/plotting logic mirrors plot_granger_rf.R line-for-line wherever the
# two-state and three-state versions coincide; please sanity-check visually
# on one known example before using it in the final figure.
# ---------------------------------------------------------------------

plot_granger_rf_shadow <- function(
    res,                        # output of granger_rf_shadow()
    y_label   = "y",
    sort_order = c("desc", "asc"),
    connected  = FALSE,
    palette = list(
      causal      = "#A8D5BA",  # pastel green  - shadow-robust
      borderline  = "#F5D08A",  # pastel amber  - passes raw rule, fails shadow test
      nocausal    = "#F4A6A6",  # pastel red    - no importance at all
      y_fill      = "#BFD7EA",  # pastel blue
      rect_border = "#555555",
      curve_color = "#777777",
      vline       = "grey85",
      sep_hline   = "grey80",
      text        = "#212121"
    ),
    sizes = list(
      var_rect_w     = 0.35,
      var_rect_h     = 0.70,
      y_rect_w       = 0.50,
      y_rect_h       = 0.90,
      hist_left      = -1.20,
      hist_max_width = 0.95,
      label_offset   = 0.06
    ),
    title   = NULL,
    margins = c(3, 4, 3, 3),
    cex_base = 1
) {
  sort_order <- match.arg(sort_order)

  # Helper: draw a smooth quadratic Bezier curve between two points (identical
  # to plot_granger_rf.R)
  draw_curve <- function(x0, y0, x1, y1, curvature = 0.25, col = palette$curve_color, lwd = 1) {
    dx <- x1 - x0; dy <- y1 - y0
    d  <- sqrt(dx*dx + dy*dy)
    if (d == 0) return(invisible())
    ux <- -dy / d; uy <- dx / d
    xm <- (x0 + x1) / 2; ym <- (y0 + y1) / 2
    offset <- curvature * d
    xc <- xm + offset * ux; yc <- ym + offset * uy
    t  <- seq(0, 1, length.out = 40)
    xb <- (1 - t)^2 * x0 + 2 * (1 - t) * t * xc + t^2 * x1
    yb <- (1 - t)^2 * y0 + 2 * (1 - t) * t * yc + t^2 * y1
    lines(xb, yb, col = col, lwd = lwd)
  }

  dec <- res$decision
  vi  <- res$variable_importance
  if (is.null(dec) || is.null(vi)) {
    stop("Input 'res' must be a granger_rf_shadow() result containing 'decision' and 'variable_importance'.")
  }

  df <- merge(
    dec[, c("variable", "best_lag", "expected_lag", "granger_precedes_shadow", "granger_precedes_raw")],
    vi[,  c("variable", "importance")],
    by = "variable", all.x = TRUE
  )
  df$importance[is.na(df$importance)] <- 0
  # res$variable_importance$importance is already the VI_total share, summed
  # to 1 across all variables with VI_total > 0 (see granger_rf_shadow.R), so
  # no further re-normalization is needed here.
  df$importance_plot <- df$importance

  # Three-way status, all derived from the SAME shadow-augmented RF fit
  df$status <- ifelse(df$granger_precedes_shadow, "Causal (shadow-robust)",
                ifelse(df$granger_precedes_raw,   "Borderline (fails shadow)",
                                                    "No causality"))
  df$rank_key <- ifelse(df$granger_precedes_shadow, 2L,
                  ifelse(df$granger_precedes_raw,   1L, 0L))

  # Use best_lag (not expected_lag) for column placement whenever the
  # variable had ANY positive raw importance (causal OR borderline). This is
  # what makes "amber" bars appear in the lag column the OLD rule would have
  # assigned them, instead of being collapsed to column 0 like genuine
  # non-causal variables.
  df$plot_lag <- ifelse(df$granger_precedes_raw, df$best_lag, 0L)

  if (sort_order == "desc") {
    ord <- order(-df$rank_key, -df$importance_plot, df$variable)
  } else {
    ord <- order(-df$rank_key,  df$importance_plot, df$variable)
  }
  df <- df[ord, , drop = FALSE]

  n_var <- nrow(df)
  if (n_var == 0) stop("No variables to plot.")
  df$ypos <- seq_len(n_var)

  max_nonzero_lag <- suppressWarnings(max(df$plot_lag[df$plot_lag > 0], na.rm = TRUE))
  if (!is.finite(max_nonzero_lag)) max_nonzero_lag <- 1L
  Lmax <- as.integer(max_nonzero_lag)

  df$x_col <- ifelse(
    df$plot_lag <= 0,
    0,
    Lmax - pmax(1L, as.integer(df$plot_lag)) + 1L
  )

  var_w <- sizes$var_rect_w
  var_h <- sizes$var_rect_h
  df$xmin <- df$x_col - var_w/2
  df$xmax <- df$x_col + var_w/2
  df$ymin <- df$ypos  - var_h/2
  df$ymax <- df$ypos  + var_h/2

  hist_left  <- sizes$hist_left
  hist_width <- sizes$hist_max_width
  df$bar_len <- hist_width * df$importance_plot
  df$bxmin   <- hist_left
  df$bxmax   <- hist_left + df$bar_len
  df$bymin   <- df$ypos - var_h * 0.35
  df$bymax   <- df$ypos + var_h * 0.35

  y_cx <- Lmax + 1.6
  y_cy <- mean(range(df$ypos))
  y_w  <- sizes$y_rect_w
  y_h  <- sizes$y_rect_h
  y_rect <- c(xmin = y_cx - y_w/2, xmax = y_cx + y_w/2,
              ymin = y_cy - y_h/2, ymax = y_cy + y_h/2)

  df_edges <- df[df$status == "Causal (shadow-robust)" & df$plot_lag > 0, , drop = FALSE]

  vlines_at <- seq(-0.5, Lmax + 0.5, by = 1)

  lag_labels <- data.frame(
    lag = 0:Lmax,
    x   = c(0, Lmax - (1:Lmax) + 1L),
    y   = n_var + 0.9,
    lab = as.character(0:Lmax)
  )

  # Two separators: Causal | Borderline | No-causality
  n_causal <- sum(df$status == "Causal (shadow-robust)")
  n_border <- sum(df$status == "Borderline (fails shadow)")
  sep_y1 <- if (n_causal > 0 && n_causal < n_var) (n_causal + 0.5) else NA_real_
  sep_y2 <- if ((n_causal + n_border) > 0 && (n_causal + n_border) < n_var) (n_causal + n_border + 0.5) else NA_real_

  oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar))
  par(mar = margins)

  x_left  <- sizes$hist_left - 0.1
  x_right <- y_cx + y_w/2 + 0.6
  y_bottom <- 0.5
  y_top    <- n_var + 1.2

  plot.new()
  plot.window(xlim = c(x_left, x_right), ylim = c(y_bottom, y_top), xaxs = "i", yaxs = "i")
  if (!is.null(title)) mtext(title, side = 3, line = 1, font = 2, cex = 1.0 * cex_base, col = palette$text)

  for (vx in vlines_at) {
    segments(vx, y_bottom, vx, y_top, col = palette$vline, lwd = 1)
  }

  fill_for <- function(s) switch(s,
    "Causal (shadow-robust)"    = palette$causal,
    "Borderline (fails shadow)" = palette$borderline,
    palette$nocausal)

  for (i in seq_len(n_var)) {
    fill_col <- fill_for(df$status[i])
    rect(df$bxmin[i], df$bymin[i], df$bxmax[i], df$bymax[i], col = fill_col, border = NA)
    lbl <- sprintf("%.1f%%", df$importance_plot[i] * 100)
    text(df$bxmax[i] + sizes$label_offset, df$ypos[i], labels = lbl,
         cex = 0.85 * cex_base, adj = c(0, 0.5), col = palette$text)
  }

  for (i in seq_len(n_var)) {
    fill_col <- fill_for(df$status[i])
    rect(df$xmin[i], df$ymin[i], df$xmax[i], df$ymax[i],
         col = fill_col, border = palette$rect_border)
    text(df$x_col[i], df$ypos[i], labels = df$variable[i],
         cex = 0.95 * cex_base, col = palette$text)
  }

  rect(y_rect["xmin"], y_rect["ymin"], y_rect["xmax"], y_rect["ymax"],
       col = palette$y_fill, border = palette$rect_border)
  text(y_cx, y_cy, labels = y_label, cex = 1.0 * cex_base, col = palette$text)

  if (connected && nrow(df_edges) > 0) {
    for (i in seq_len(nrow(df_edges))) {
      draw_curve(df_edges$xmax[i], df_edges$ypos[i], y_rect["xmin"], y_cy,
                 curvature = 0.25, col = palette$curve_color, lwd = 1)
    }
  }

  text(lag_labels$x, lag_labels$y, labels = lag_labels$lab,
       cex = 0.9 * cex_base, col = "#444444")

  if (is.finite(sep_y1)) segments(x_left, sep_y1, x_right, sep_y1, col = palette$sep_hline, lwd = 2/3)
  if (is.finite(sep_y2)) segments(x_left, sep_y2, x_right, sep_y2, col = palette$sep_hline, lwd = 2/3, lty = 2)

  invisible(list(
    df = df,
    xlim = c(x_left, x_right),
    ylim = c(y_bottom, y_top),
    y_rect = y_rect
  ))
}

# ---------------------------------------------------------------------
# Suggested legend text for the figure caption/note, e.g.:
#   "Green: Granger-precedes under the shadow-feature test; amber: positive
#    importance but does not exceed its own shadow (permuted) copies, i.e.
#    would have been called causal under the pre-correction rule; red: no
#    importance. Percentages are each variable's share of total positive
#    variable importance."
# ---------------------------------------------------------------------
