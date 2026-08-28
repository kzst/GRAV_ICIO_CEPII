plot_granger_vars <- function(
    res,
    view = c("main", "graph", "both"),
    main_title = NULL,
    graph_title = NULL,
    y_label = "y",
    sort_order = c("desc", "asc"),
    connected = FALSE,
    dep_graph_layout = c("tree", "lag", "fr", "kk", "lgl", "circle"),
    dep_graph_root = NULL,
    palette = list(
      causal      = "#A8D5BA",
      nocausal    = "#F4A6A6",
      y_fill      = "#BFD7EA",
      rect_border = "#555555",
      curve_color = "#777777",
      vline       = "grey85",
      sep_hline   = "grey80",
      text        = "#212121",
      bic_color   = "darkred"
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
    margins = c(3, 4, 3, 3),
    cex_base = 1
) {
  view <- match.arg(view)
  sort_order <- match.arg(sort_order)
  # dep_graph_layout lehet előre definiált választás, string, vagy függvény
  dep_layout_choice <- dep_graph_layout[1]
  
  # ---- Extract inputs ----
  expected <- res$expected_lags
  vi <- res$variable_importance
  if (is.null(expected) || nrow(expected) == 0) stop("res$expected_lags is missing or empty.")
  if (is.null(vi)) vi <- data.frame(variable = expected$variable, importance = 0)
  
  # Merge and sanitize
  df <- merge(
    expected[, c("variable", "expected_lag", "drop", "p_value")],
    vi[, c("variable", "importance")],
    by = "variable", all.x = TRUE
  )
  df$importance[!is.finite(df$importance)] <- 0
  df$importance <- as.numeric(df$importance)
  df$is_causal_plot <- !df$drop
  
  # Ordering
  if (sort_order == "desc") {
    ord <- order(-as.integer(df$is_causal_plot), -df$importance, df$variable)
  } else {
    ord <- order(-as.integer(df$is_causal_plot),  df$importance, df$variable)
  }
  df <- df[ord, , drop = FALSE]
  n_var <- nrow(df)
  df$ypos <- seq_len(n_var)
  
  # Main panel geometry
  Lmax <- max(c(1, df$expected_lag), na.rm = TRUE)
  df$x_col <- ifelse(df$expected_lag <= 0, 0L, Lmax - pmax(1L, as.integer(df$expected_lag)) + 1L)
  var_w <- sizes$var_rect_w; var_h <- sizes$var_rect_h
  df$xmin <- df$x_col - var_w/2; df$xmax <- df$x_col + var_w/2
  df$ymin <- df$ypos - var_h/2; df$ymax <- df$ypos + var_h/2
  hist_left <- sizes$hist_left; hist_width <- sizes$hist_max_width
  df$bar_len <- hist_width * df$importance
  df$bxmin <- hist_left; df$bxmax <- hist_left + df$bar_len
  df$bymin <- df$ypos - var_h * 0.35; df$bymax <- df$ypos + var_h * 0.35
  y_cx <- Lmax + 1.6; y_cy <- mean(range(df$ypos))
  y_w  <- sizes$y_rect_w; y_h  <- sizes$y_rect_h
  y_rect <- c(xmin = y_cx - y_w/2, xmax = y_cx + y_w/2, ymin = y_cy - y_h/2,  ymax = y_cy + y_h/2)
  vlines_at <- seq(-0.5, Lmax + 0.5, by = 1)
  n_causal <- sum(df$is_causal_plot)
  sep_y <- if (n_causal > 0 && n_causal < n_var) (n_causal + 0.5) else NA_real_
  
  draw_curve <- function(x0, y0, x1, y1, curvature = 0.25, col = palette$curve_color, lwd = 1) {
    dx <- x1 - x0; dy <- y1 - y0; d <- sqrt(dx*dx + dy*dy); if (d == 0) return(invisible())
    ux <- -dy / d; uy <- dx / d
    xm <- (x0 + x1)/2; ym <- (y0 + y1)/2
    offset <- curvature * d
    xc <- xm + offset * ux; yc <- ym + offset * uy
    t  <- seq(0, 1, length.out = 40)
    xb <- (1 - t)^2 * x0 + 2 * (1 - t) * t * xc + t^2 * x1
    yb <- (1 - t)^2 * y0 + 2 * (1 - t) * t * yc + t^2 * y1
    lines(xb, yb, col = col, lwd = lwd)
  }
  
  oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar))
  if (view == "both") par(mfrow = c(1, 2))
  
  # ---- Main panel ----
  if (view %in% c("main", "both")) {
    par(mar = margins)
    x_left  <- sizes$hist_left - 0.1
    x_right <- y_cx + y_w/2 + 0.6
    y_bottom <- 0.5
    y_top    <- n_var + 1.2
    
    plot.new()
    plot.window(xlim = c(x_left, x_right), ylim = c(y_bottom, y_top), xaxs = "i", yaxs = "i")
    if (!is.null(main_title)) mtext(main_title, side = 3, line = 1, font = 2, cex = 1.0 * cex_base, col = palette$text)
    
    for (vx in vlines_at) segments(vx, y_bottom, vx, y_top, col = palette$vline, lwd = 1)
    
    # Left bars
    for (i in seq_len(n_var)) {
      fill_col <- if (df$is_causal_plot[i]) palette$causal else palette$nocausal
      rect(df$bxmin[i], df$bymin[i], df$bxmax[i], df$bymax[i], col = fill_col, border = NA)
      lbl <- if (df$importance[i] > 0) sprintf("%.1f%%", df$importance[i] * 100) else "0.0%"
      text(df$bxmax[i] + sizes$label_offset, df$ypos[i], labels = lbl,
           cex = 0.85 * cex_base, adj = c(0, 0.5), col = palette$text)
    }
    
    # Variable rectangles
    for (i in seq_len(n_var)) {
      fill_col <- if (df$is_causal_plot[i]) palette$causal else palette$nocausal
      rect(df$xmin[i], df$ymin[i], df$xmax[i], df$ymax[i], col = fill_col, border = palette$rect_border)
      text(df$x_col[i], df$ypos[i], labels = df$variable[i], cex = 0.95 * cex_base, col = palette$text)
    }
    
    # y node
    rect(y_rect["xmin"], y_rect["ymin"], y_rect["xmax"], y_rect["ymax"], col = palette$y_fill, border = palette$rect_border)
    text(y_cx, y_cy, labels = y_label, cex = 1.0 * cex_base, col = palette$text)
    
    if (isTRUE(connected)) {
      for (i in which(df$is_causal_plot)) {
        x0 <- df$xmax[i]; y0 <- df$ypos[i]
        x1 <- y_rect["xmin"]; y1 <- y_cy
        draw_curve(x0, y0, x1, y1, curvature = 0.22, col = palette$curve_color, lwd = 1)
      }
    }
    
    # Lag labels and separator
    lag_labels <- data.frame(lag = 0:Lmax, x = c(0, Lmax - (1:Lmax) + 1L), y = n_var + 0.9, lab = as.character(0:Lmax))
    text(lag_labels$x, lag_labels$y, labels = lag_labels$lab, cex = 0.9 * cex_base, col = "#444444")
    if (is.finite(sep_y)) segments(x_left, sep_y, x_right, sep_y, col = palette$sep_hline, lwd = 1)
  }
  
  # ---- Explanatory graph panel ----
  if (view %in% c("graph", "both")) {
    g <- res$explanatory_graph
    
    BIN_GRANGER<-as.matrix(as_adjacency_matrix(g))
    rownames(BIN_GRANGER)<-colnames(BIN_GRANGER)<-rownames(as.matrix(V(g)))
    resG <- iBBiG(binaryMatrix=BIN_GRANGER,nModules = 1,alpha=0.3,pop_size = 100,mutation = 0.3,stagnation = 50,selection_pressure = 1.2,max_sp = 15,success_ratio = 0.8)
    ADJ_GRANGER<-BIN_GRANGER[resG@RowxNumber,resG@NumberxCol]
    
    
    
    nd_all <- df[, c("variable", "is_causal_plot", "expected_lag", "x_col", "ypos", "importance")]
    nd_all$importance[!is.finite(nd_all$importance)] <- 0
    nd_all$importance <- as.numeric(nd_all$importance)
    
    ed <- NULL
    if (!is.null(g)) ed <- igraph::as_data_frame(g, what = "edges")
    nodes_use <- character(0)
    if (!is.null(ed) && nrow(ed) > 0) nodes_use <- unique(c(as.character(ed$from), as.character(ed$to)))
    nd <- nd_all[nd_all$variable %in% nodes_use, , drop = FALSE]
    has_nodes <- nrow(nd) > 0
    
    par(mar = margins)
    plot.new()
    if (!is.null(graph_title)) mtext(graph_title, side = 3, line = 1, font = 2, cex = 1.0 * cex_base, col = palette$text)
    
    if (!has_nodes) {
      plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
      text(0.5, 0.5, labels = "No edges among explanatory variables", cex = 0.95 * cex_base, col = "#666666")
      return(invisible(list(df_main = df, y_rect = y_rect, Lmax = Lmax, view = view)))
    }
    
    # Build induced subgraph in node order
    node_order <- nd$variable
    g_sub <- igraph::induced_subgraph(g, vids = node_order)
    
    if (dep_layout_choice == "lag") {
      # Lag-based columns layout (as before)
      if (sort_order == "desc") {
        nd <- nd[order(nd$x_col, -as.integer(nd$is_causal_plot), -nd$importance, nd$variable), ]
      } else {
        nd <- nd[order(nd$x_col, -as.integer(nd$is_causal_plot),  nd$importance, nd$variable), ]
      }
      nd$slot_in_col <- ave(seq_len(nrow(nd)), nd$x_col, FUN = seq_along)
      max_slots <- tapply(nd$slot_in_col, nd$x_col, max); max_slots[is.na(max_slots)] <- 0
      col_heights <- pmax(1, max_slots)
      nd$ypos_g <- ave(nd$slot_in_col, nd$x_col, FUN = function(s) {
        s <- as.numeric(s); m <- max(s); if (!is.finite(m) || m == 0) return(rep(1, length(s)))
        (m + 1) - s + 0.5
      })
      max_col <- max(nd$x_col, na.rm = TRUE); if (!is.finite(max_col)) max_col <- 1L
      max_y_stack <- max(col_heights, na.rm = TRUE); if (!is.finite(max_y_stack)) max_y_stack <- 1
      
      x_left_g  <- -0.6
      x_right_g <- max_col + 0.6
      y_bottom_g <- 0.5
      y_top_g    <- max_y_stack + 1.2
      plot.window(xlim = c(x_left_g, x_right_g), ylim = c(y_bottom_g, y_top_g), xaxs = "i", yaxs = "i")
      
      vlines_at_g <- seq(-0.5, max_col + 0.5, by = 1)
      for (vx in vlines_at_g) segments(vx, y_bottom_g, vx, y_top_g, col = palette$vline, lwd = 1)
      lag_labels_g <- data.frame(lag = 0:Lmax, x = c(0, Lmax - (1:Lmax) + 1L), y = y_top_g - 0.2, lab = as.character(0:Lmax))
      text(lag_labels_g$x, lag_labels_g$y, labels = lag_labels_g$lab, cex = 0.9 * cex_base, col = "#444444")
      
      var_w <- sizes$var_rect_w; var_h <- sizes$var_rect_h
      nd$xmin_g <- nd$x_col - var_w/2; nd$xmax_g <- nd$x_col + var_w/2
      nd$ymin_g <- nd$ypos_g - var_h/2; nd$ymax_g <- nd$ypos_g + var_h/2
      
      ed_sub <- ed[ed$from %in% nd$variable & ed$to %in% nd$variable, , drop = FALSE]
      if (nrow(ed_sub) > 0) {
        pos_map <- setNames(seq_len(nrow(nd)), nd$variable)
        for (k in seq_len(nrow(ed_sub))) {
          from <- as.character(ed_sub$from[k]); to <- as.character(ed_sub$to[k])
          i <- pos_map[[from]]; j <- pos_map[[to]]
          x0 <- nd$xmax_g[i]; y0 <- (nd$ymin_g[i] + nd$ymax_g[i]) / 2
          x1 <- nd$xmin_g[j]; y1 <- (nd$ymin_g[j] + nd$ymax_g[j]) / 2
          if ((from %in% rownames(ADJ_GRANGER))&&(to %in% colnames(ADJ_GRANGER))){
            w<-5
            lwd_edge <- 0.8 + 0.4 * as.numeric(w)
            arrows(x0, y0, x1, y1, length = 0.07, lwd = lwd_edge, col = palette$bic_color, angle = 20)
          }else{
            w<-1
            lwd_edge <- 0.8 + 0.4 * as.numeric(w)
            arrows(x0, y0, x1, y1, length = 0.07, lwd = lwd_edge, col = palette$curve_color, angle = 20)
          }
          w <- ed_sub$weight[k]; if (!is.finite(w)) w <- 1
          #lwd_edge <- 0.8 + 0.4 * as.numeric(w)
          #arrows(x0, y0, x1, y1, length = 0.07, lwd = lwd_edge, col = palette$curve_color, angle = 20)
          xm <- (x0 + x1)/2; ym <- (y0 + y1)/2
          text(xm, ym + 0.05, labels = as.character(w), cex = 0.75 * cex_base, col = "#444444")
        }
      }
      
      for (i in seq_len(nrow(nd))) {
        fill_col <- if (nd$is_causal_plot[i]) palette$causal else palette$nocausal
        rect(nd$xmin_g[i], nd$ymin_g[i], nd$xmax_g[i], nd$ymax_g[i], col = fill_col, border = palette$rect_border)
        text(nd$x_col[i], nd$ypos_g[i], labels = nd$variable[i], cex = 0.9 * cex_base, col = palette$text)
      }
      
    } else {
      # --- igraph-based generic layouts (including tree and arbitrary user-provided) ---
      # Resolve layout function
      layout_fun <- NULL
      if (is.function(dep_layout_choice)) {
        layout_fun <- dep_layout_choice
      } else {
        # Predefined shortcuts
        if (dep_layout_choice %in% c("tree","fr","kk","lgl","circle")) {
          if (dep_layout_choice == "tree")       layout_fun <- function(g) igraph::layout_as_tree(g, root = {
            roots <- dep_graph_root
            if (is.null(roots)) {
              indeg <- igraph::degree(g, mode = "in")
              roots <- names(indeg)[indeg == 0]
              if (length(roots) == 0) roots <- igraph::V(g)$name[1]
            }
            igraph::V(g)[igraph::V(g)$name %in% roots]
          })
          if (dep_layout_choice == "fr")         layout_fun <- igraph::layout_with_fr
          if (dep_layout_choice == "kk")         layout_fun <- igraph::layout_with_kk
          if (dep_layout_choice == "lgl")        layout_fun <- igraph::layout_with_lgl
          if (dep_layout_choice == "circle")     layout_fun <- igraph::layout_in_circle
        } else if (is.character(dep_layout_choice)) {
          # Try exact name in igraph
          if (exists(dep_layout_choice, where = asNamespace("igraph"), inherits = FALSE)) {
            layout_fun <- get(dep_layout_choice, envir = asNamespace("igraph"))
          } else {
            # Try with "layout_" prefix variations
            trial_names <- c(paste0("layout_", dep_layout_choice),
                             paste0("layout_with_", dep_layout_choice),
                             paste0("layout_as_", dep_layout_choice),
                             paste0("layout_in_", dep_layout_choice),
                             "layout_nicely")
            for (nm in trial_names) {
              if (exists(nm, where = asNamespace("igraph"), inherits = FALSE)) {
                layout_fun <- get(nm, envir = asNamespace("igraph"))
                break
              }
            }
          }
          if (is.null(layout_fun)) layout_fun <- igraph::layout_nicely
        } else {
          layout_fun <- igraph::layout_nicely
        }
      }
      
      coords <- layout_fun(g_sub)
      
      # Normalize coords -> [0,1], then pad by node box half-width/height so nothing clips
      xr <- range(coords[,1], na.rm = TRUE); yr <- range(coords[,2], na.rm = TRUE)
      if (diff(xr) == 0) xr <- xr + c(-0.5, 0.5)
      if (diff(yr) == 0) yr <- yr + c(-0.5, 0.5)
      x_norm <- (coords[,1] - xr[1]) / diff(xr)
      y_norm <- (coords[,2] - yr[1]) / diff(yr)
      
      # Node positions
      node_names <- igraph::V(g_sub)$name
      pos_map <- setNames(seq_along(node_names), node_names)
      nd_pos_idx <- pos_map[nd$variable]
      nd$xg <- x_norm[nd_pos_idx]
      nd$yg <- y_norm[nd_pos_idx]
      
      # Box size in normalized units and padding
      w_norm <- 0.10; h_norm <- 0.08
      pad_x <- w_norm/2 + 0.04
      pad_y <- h_norm/2 + 0.04
      xlim <- c(0 - pad_x, 1 + pad_x)
      ylim <- c(0 - pad_y, 1 + pad_y)
      
      plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i")
      
      # Edges
      ed_sub <- ed[ed$from %in% nd$variable & ed$to %in% nd$variable, , drop = FALSE]
      if (nrow(ed_sub) > 0) {
        for (k in seq_len(nrow(ed_sub))) {
          from <- as.character(ed_sub$from[k]); to <- as.character(ed_sub$to[k])
          i <- pos_map[[from]]; j <- pos_map[[to]]
          x0 <- x_norm[i]; y0 <- y_norm[i]
          x1 <- x_norm[j]; y1 <- y_norm[j]
          if ((from %in% rownames(ADJ_GRANGER))&&(to %in% colnames(ADJ_GRANGER))){
            w<-5
            lwd_edge <- 0.8 + 0.4 * as.numeric(w)
            arrows(x0, y0, x1, y1, length = 0.06, lwd = lwd_edge, col = palette$bic_color, angle = 20)
          }else{
            w<-1
            lwd_edge <- 0.8 + 0.4 * as.numeric(w)
            arrows(x0, y0, x1, y1, length = 0.06, lwd = lwd_edge, col = palette$curve_color, angle = 20)
          }
          w <- ed_sub$weight[k]; if (!is.finite(w)) w <- 1
          #lwd_edge <- 0.8 + 0.4 * as.numeric(w)
          #arrows(x0, y0, x1, y1, length = 0.06, lwd = lwd_edge, col = palette$curve_color, angle = 20)
          text((x0 + x1)/2, (y0 + y1)/2 + 0.02, labels = as.character(w), cex = 0.75 * cex_base, col = "#444444")
        }
      }
      
      # Nodes
      nd$xmin_g <- nd$xg - w_norm/2; nd$xmax_g <- nd$xg + w_norm/2
      nd$ymin_g <- nd$yg - h_norm/2; nd$ymax_g <- nd$yg + h_norm/2
      for (i in seq_len(nrow(nd))) {
        fill_col <- if (nd$is_causal_plot[i]) palette$causal else palette$nocausal
        rect(nd$xmin_g[i], nd$ymin_g[i], nd$xmax_g[i], nd$ymax_g[i], col = fill_col, border = palette$rect_border)
        text(nd$xg[i], nd$yg[i], labels = nd$variable[i], cex = 0.9 * cex_base, col = palette$text)
      }
    }
    
    
  }
  
  invisible(list(df_main = df, y_rect = y_rect, Lmax = Lmax, view = view, dep_graph_layout = dep_layout_choice))
}