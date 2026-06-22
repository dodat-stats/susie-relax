#' Plot SuSiE Alpha (PIP) with Credible Sets
#'
#' Visualizes posterior inclusion probabilities (PIPs) and overlays
#' identified credible sets with distinct colors. Optionally highlights
#' true causal effects.
#'
#' @param alpha L x J matrix of posterior inclusion weights, a susieR fit, or
#'   a SuSiE-relax-lite fit returned by `fit_susie_relax_lite()`.
#' @param R Optional J x J covariance matrix. If provided, credible sets are
#'   filtered by a purity threshold.
#' @param b Optional numeric vector of length J indicating true effects
#'   (nonzero entries are highlighted in red).
#' @param main Character. Plot title.
#' @param coverage Numeric. Coverage threshold for credible sets (default 0.95).
#' @param min_abs_corr Numeric. Minimum absolute correlation for purity
#'   filtering (default 0.5).
#' @param ... Additional arguments passed to [plot()].
#'
#' @return Invisibly returns a list of credible sets.
#'
#' @export
susie_plot_alpha <- function(alpha, R = NULL, b = NULL,
                             main = "SuSiE PIP Visualization",
                             coverage = 0.95, min_abs_corr = 0.5, ...) {

  model <- alpha
  if (is.matrix(model)) {
    alpha <- model
    p <- 1 - apply(1 - alpha, 2, prod)
  } else if (is.list(model) && !is.null(model$alpha)) {
    alpha <- as.matrix(model$alpha)
    if (!is.null(model$pip)) {
      p <- as.numeric(model$pip)
    } else {
      p <- 1 - apply(1 - alpha, 2, prod)
    }
  } else {
    stop("alpha must be an L x p matrix or a fitted object with an alpha matrix")
  }

  # Identify credible sets (CS) and calculate purity.
  cs_list <- list()
  purity_list <- numeric()

  for (l in 1:nrow(alpha)) {
    alpha_l <- alpha[l, ]
    ord <- order(alpha_l, decreasing = TRUE)
    cum_alpha <- cumsum(alpha_l[ord])
    idx <- which(cum_alpha >= coverage)[1]

    if (max(alpha_l) > 1e-4) {
      cs_vars <- ord[1:idx]

      if (!is.null(R)) {
        if (length(cs_vars) == 1) {
          purity <- 1.0
        } else {
          R_sub <- abs(R[cs_vars, cs_vars])
          purity <- min(R_sub[upper.tri(R_sub)])
        }
        if (purity >= min_abs_corr) {
          cs_list[[paste0("L", l)]] <- cs_vars
          purity_list[paste0("L", l)] <- purity
        }
      } else {
        cs_list[[paste0("L", l)]] <- cs_vars
        purity_list[paste0("L", l)] <- NA_real_
      }
    }
  }

  # Setup base plot parameters.
  pos <- 1:length(p)
  if (is.null(b)) {
    b <- rep(0, length(p))
  }
  if (length(b) != length(p)) {
    stop("b must have the same length as the PIP vector")
  }

  args <- list(...)
  if (!("xlab" %in% names(args))) args$xlab <- "variable"
  if (!("ylab" %in% names(args))) args$ylab <- "PIP"
  if (!("pch" %in% names(args))) args$pch <- 16
  if (!("ylim" %in% names(args))) args$ylim <- c(0, max(1, p, na.rm = TRUE))
  args$x <- pos
  args$y <- p
  args$main <- main

  # Draw base PIP plot.
  do.call(graphics::plot, args)

  color <- c(
    "dodgerblue2", "green4", "#6A3D9A", "#FF7F00", "gold1",
    "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F",
    "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1",
    "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )

  legend_text <- list(col = vector(), cs_label = vector())

  # Overlay credible set colors.
  for (i in rev(seq_along(cs_list))) {
    cs_vars <- cs_list[[i]]
    cs_name <- names(cs_list)[i]

    if (length(cs_vars) > 0) {
      col_to_use <- head(color, 1)
      graphics::points(pos[cs_vars], p[cs_vars], col = col_to_use, cex = 1.5, lwd = 2.5)

      legend_text$col <- append(col_to_use, legend_text$col)
      label <- if (is.na(purity_list[cs_name])) {
        sprintf("%s: C=%d", cs_name, length(cs_vars))
      } else {
        sprintf("%s: C=%d/R=%.3f", cs_name, length(cs_vars), purity_list[cs_name])
      }
      legend_text$cs_label <- append(label, legend_text$cs_label)

      color <- c(color[-1], color[1])
    }
  }

  # Legend.
  if (length(legend_text$col) > 0) {
    graphics::legend("topright",
                     legend = legend_text$cs_label,
                     bty = "n",
                     col = legend_text$col,
                     cex = 0.65,
                     pch = 15)
  }

  # Highlight true effects.
  causal <- which(b != 0)
  if (length(causal) > 0) {
    graphics::points(pos[causal], p[causal], col = 2, pch = 16)
    graphics::rug(pos[causal], col = 2, lwd = 1.5)
  }

  return(invisible(list(pip = p, cs = cs_list, purity = purity_list)))
}

#' Plot SuSiE-relax-lite PIPs with a susieR-like interface.
#'
#' @param model Fit returned by `fit_susie_relax_lite()`.
#' @param y Currently supports `"PIP"` and `"log10PIP"`.
#' @param R Optional LD/correlation matrix for credible-set purity.
#' @param b Optional true effect vector; nonzero entries are highlighted.
#' @param ... Additional arguments passed to [plot()].
#'
#' @return Invisibly returns the plotted PIP vector and credible sets.
susie_plot_relax_lite <- function(model, y = "PIP", R = NULL, b = NULL,
                                  main = NULL, coverage = 0.95,
                                  min_abs_corr = 0.5, ...) {
  if (!is.list(model) || is.null(model$alpha)) {
    stop("model must be a fit returned by fit_susie_relax_lite()")
  }

  if (is.null(main)) {
    main <- if (y == "log10PIP") {
      "SuSiE-relax-lite log10(PIP)"
    } else {
      "SuSiE-relax-lite PIP"
    }
  }

  if (y == "PIP") {
    susie_plot_alpha(
      model, R = R, b = b, main = main,
      coverage = coverage, min_abs_corr = min_abs_corr, ...
    )
  } else if (y == "log10PIP") {
    pip <- as.numeric(model$pip)
    graphics::plot(
      seq_along(pip), log10(pip),
      xlab = "variable",
      ylab = "log10(PIP)",
      pch = 16,
      main = main,
      ...
    )
    if (!is.null(b)) {
      causal <- which(b != 0)
      graphics::points(causal, log10(pip[causal]), col = 2, pch = 16)
      graphics::rug(causal, col = 2, lwd = 1.5)
    }
    invisible(list(pip = pip))
  } else {
    stop("Only y = 'PIP' and y = 'log10PIP' are supported")
  }
}
