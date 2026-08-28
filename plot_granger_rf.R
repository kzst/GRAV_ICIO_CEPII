plot_granger_rf <- function(
    res,
    y_label   = "y",
    sort_order = c("desc", "asc"),
    connected  = FALSE,
    palette = list(
      causal      = "#A8D5BA",  # pastel green
      nocausal    = "#F4A6A6",  # pastel red
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
  
  # Helper: draw a smooth quadratic Bezier curve between two points
  draw_curve <- function(x0, y0, x1, y1, curvature = 0.25, col = palette$curve_color, lwd = 1) {
    dx <- x1 - x0; dy <- y1 - y0
    d  <- sqrt(dx*dx + dy*dy)
    if (d == 0) return(invisible())
    # Perpendicular unit vector
    ux <- -dy / d; uy <- dx / d
    # Control point at mid + offset along perpendicular
    xm <- (x0 + x1) / 2; ym <- (y0 + y1) / 2
    offset <- curvature * d
    xc <- xm + offset * ux; yc <- ym + offset * uy
    # Quadratic Bezier sampling
    t  <- seq(0, 1, length.out = 40)
    xb <- (1 - t)^2 * x0 + 2 * (1 - t) * t * xc + t^2 * x1
    yb <- (1 - t)^2 * y0 + 2 * (1 - t) * t * yc + t^2 * y1
    lines(xb, yb, col = col, lwd = lwd)
  }
  
  # Extract and merge data
  df_exp <- res$expected_lags
  df_vi  <- res$variable_importance
  if (is.null(df_exp) || is.null(df_vi)) {
    stop("Input 'res' must contain 'expected_lags' and 'variable_importance'.")
  }
  
  df <- merge(
    df_exp[, c("variable", "expected_lag", "drop")],
    df_vi[,  c("variable", "importance")],
    by = "variable", all.x = TRUE
  )
  df$importance[is.na(df$importance)] <- 0
  
  # Raw importance (may contain negatives from permutation importance)
  df$importance_raw <- df$importance
  
  # Renormalize strictly positive importances to sum to 1 (100%)
  pos_total <- sum(df$importance_raw[df$importance_raw > 0], na.rm = TRUE)
  if (pos_total > 0) {
    df$importance_plot <- ifelse(df$importance_raw > 0, df$importance_raw / pos_total, 0)
  } else {
    df$importance_plot <- 0
  }
  
  # Causality flag for plotting: red if drop==TRUE or importance_raw<=0
  df$is_causal_plot <- (!df$drop) & (df$importance_raw > 0)
  df$status <- ifelse(df$is_causal_plot, "Causal", "No causality")
  
  # Order: causal group on top, non-causal bottom; within group sort by importance asc/desc
  if (sort_order == "desc") {
    ord <- order(-as.integer(df$is_causal_plot), -df$importance_plot, df$variable)
  } else {
    ord <- order(-as.integer(df$is_causal_plot),  df$importance_plot, df$variable)
  }
  df <- df[ord, , drop = FALSE]
  
  n_var <- nrow(df)
  if (n_var == 0) stop("No variables to plot.")
  df$ypos <- seq_len(n_var)
  
  # Lag columns: 0 leftmost (non-causal), then 1..Lmax towards y on the right
  max_nonzero_lag <- suppressWarnings(max(df$expected_lag[df$expected_lag > 0], na.rm = TRUE))
  if (!is.finite(max_nonzero_lag)) max_nonzero_lag <- 1L
  Lmax <- as.integer(max_nonzero_lag)
  
  df$x_col <- ifelse(
    df$expected_lag <= 0,
    0,
    Lmax - pmax(1L, as.integer(df$expected_lag)) + 1L
  )
  
  # Rectangles for variables
  var_w <- sizes$var_rect_w
  var_h <- sizes$var_rect_h
  df$xmin <- df$x_col - var_w/2
  df$xmax <- df$x_col + var_w/2
  df$ymin <- df$ypos  - var_h/2
  df$ymax <- df$ypos  + var_h/2
  
  # Histogram bars (left), lengths from importance_plot (sum of positives = 1)
  hist_left  <- sizes$hist_left
  hist_width <- sizes$hist_max_width
  df$bar_len <- hist_width * df$importance_plot
  df$bxmin   <- hist_left
  df$bxmax   <- hist_left + df$bar_len
  df$bymin   <- df$ypos - var_h * 0.35
  df$bymax   <- df$ypos - 0 + var_h * 0.35
  
  # y-node rectangle (right side)
  y_cx <- Lmax + 1.6
  y_cy <- mean(range(df$ypos))
  y_w  <- sizes$y_rect_w
  y_h  <- sizes$y_rect_h
  y_rect <- c(xmin = y_cx - y_w/2, xmax = y_cx + y_w/2,
              ymin = y_cy - y_h/2, ymax = y_cy + y_h/2)
  
  # Edges for causal variables only if connected = TRUE
  df_edges <- df[df$is_causal_plot & df$expected_lag > 0, , drop = FALSE]
  
  # Vertical lines for columns -0.5 ... Lmax + 0.5
  vlines_at <- seq(-0.5, Lmax + 0.5, by = 1)
  
  # Lag labels above columns (0..Lmax)
  lag_labels <- data.frame(
    lag = 0:Lmax,
    x   = c(0, Lmax - (1:Lmax) + 1L),
    y   = n_var + 0.9,
    lab = as.character(0:Lmax)
  )
  
  # Horizontal separator between groups
  n_causal <- sum(df$is_causal_plot)
  sep_y <- if (n_causal > 0 && n_causal < n_var) (n_causal + 0.5) else NA_real_
  
  # Prepare plotting region
  oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar))
  par(mar = margins)
  
  x_left  <- sizes$hist_left - 0.1
  x_right <- y_cx + y_w/2 + 0.6
  y_bottom <- 0.5
  y_top    <- n_var + 1.2
  
  plot.new()
  plot.window(xlim = c(x_left, x_right), ylim = c(y_bottom, y_top), xaxs = "i", yaxs = "i")
  # Title if any
  if (!is.null(title)) mtext(title, side = 3, line = 1, font = 2, cex = 1.0 * cex_base, col = palette$text)
  
  # Vertical separators (lag columns)
  for (vx in vlines_at) {
    segments(vx, y_bottom, vx, y_top, col = palette$vline, lwd = 1)
  }
  
  # Left histogram bars and labels
  for (i in seq_len(n_var)) {
    fill_col <- if (df$is_causal_plot[i]) palette$causal else palette$nocausal
    rect(df$bxmin[i], df$bymin[i], df$bxmax[i], df$bymax[i],
         col = fill_col, border = NA)
    # percentage label
    lbl <- sprintf("%.1f%%", df$importance_plot[i] * 100)
    text(df$bxmax[i] + sizes$label_offset, df$ypos[i], labels = lbl,
         cex = 0.85 * cex_base, adj = c(0, 0.5), col = palette$text)
  }
  
  # Variable rectangles and labels
  for (i in seq_len(n_var)) {
    fill_col <- if (df$is_causal_plot[i]) palette$causal else palette$nocausal
    rect(df$xmin[i], df$ymin[i], df$xmax[i], df$ymax[i],
         col = fill_col, border = palette$rect_border)
    text(df$x_col[i], df$ypos[i], labels = df$variable[i],
         cex = 0.95 * cex_base, col = palette$text)
  }
  
  # y-node rectangle and label
  rect(y_rect["xmin"], y_rect["ymin"], y_rect["xmax"], y_rect["ymax"],
       col = palette$y_fill, border = palette$rect_border)
  text(y_cx, y_cy, labels = y_label, cex = 1.0 * cex_base, col = palette$text)
  
  # Curved connections (optional)
  if (connected && nrow(df_edges) > 0) {
    for (i in seq_len(nrow(df_edges))) {
      x0 <- df_edges$xmax[i]
      y0 <- df_edges$ypos[i]
      x1 <- y_rect["xmin"]
      y1 <- y_cy
      draw_curve(x0, y0, x1, y1, curvature = 0.25, col = palette$curve_color, lwd = 1)
    }
  }
  
  # Lag labels above columns
  text(lag_labels$x, lag_labels$y, labels = lag_labels$lab,
       cex = 0.9 * cex_base, col = "#444444")
  
  # Horizontal separator between causal and non-causal groups
  if (is.finite(sep_y)) {
    segments(x_left, sep_y, x_right, sep_y, col = palette$sep_hline, lwd = 2/3)
  }
  
  invisible(list(
    df = df,
    xlim = c(x_left, x_right),
    ylim = c(y_bottom, y_top),
    y_rect = y_rect
  ))
}