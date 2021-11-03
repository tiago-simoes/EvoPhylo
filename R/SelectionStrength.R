#t-tests
f1 <- function(RatesByClade, posterior.clockrate.all.) {

  post.df <- length(posterior.clockrate.all.) - 1
  post.mean <- mean(posterior.clockrate.all.)

  RatesByClade$abs_rate <- RatesByClade$rate * post.mean

  post.se <- sd(posterior.clockrate.all.)/sqrt(length(posterior.clockrate.all.))
  post.ts <- abs(post.mean - RatesByClade$abs_rate)/post.se
  pvals <- 2*pt(post.ts, df = post.df, lower.tail = FALSE)

  out <- data.frame(RatesByClade$clade,
                    RatesByClade$nodes,
                    RatesByClade$clock,
                    RatesByClade$rate,
                    RatesByClade$abs_rate,
                    RatesByClade$abs_rate,
                    pvals)
  names(out) <- c("clade", "nodes", "clock", "Relative rate", "Absolute rate (mean)", "null", "p.value")
  out
}

#Plot tree with colored thresholds
f2 <- function(tree, posterior.clockrate.all., clock = 1, summary = "mean", threshold = c("1 SD", "2 SD"),
               low = "blue", mid = "gray90", high = "red", xlim = NULL,
               geo_skip = c("Quaternary", "Holocene", "Late Pleistocene")) {

  #Process threshold
  if (!is.character(threshold)) stop("'threshold' must be a character vector.", call. = FALSE)
  thresh_conf <- endsWith(threshold, "%")
  thresh_sd <- endsWith(tolower(threshold), "sd")

  if (any(!thresh_conf & !thresh_sd)) stop("All entries in 'threshold' must end in '%' for confidence intervals or 'SD' for standard deviations.", call. = FALSE)

  thresh_vals <- character(length(threshold))

  if (any(thresh_conf)) thresh_vals[thresh_conf] <- substring(threshold[thresh_conf], 1, nchar(threshold[thresh_conf]) - 1)
  if (any(thresh_sd)) thresh_vals[thresh_sd] <- substring(threshold[thresh_sd], 1, nchar(threshold[thresh_sd]) - 2)

  if (anyNA(suppressWarnings(as.numeric(thresh_vals)))) {
    stop("All entries in 'threshold' must be confidence levels or the number of standard deviations to use as the thresholds.", call. = FALSE)
  }
  thresh_vals <- as.numeric(thresh_vals)

  if (any(thresh_conf)) {
    n <- length(posterior.clockrate.all.)
    tcrits <- qt(.5*(1 + thresh_vals[thresh_conf]/100), n - 1)
  }

  breaks <- numeric(2*length(threshold))
  if (any(thresh_sd)) breaks[c(thresh_sd, thresh_sd)] <- mean(posterior.clockrate.all.) + c(-thresh_vals[thresh_sd], thresh_vals[thresh_sd]) * sd(posterior.clockrate.all.)
  if (any(thresh_conf)) breaks[c(thresh_conf, thresh_conf)] <- mean(posterior.clockrate.all.) + c(-tcrits, tcrits)*sd(posterior.clockrate.all.)/sqrt(n)
  breaks <- breaks/mean(posterior.clockrate.all.)

  labels <- character(2*length(threshold))
  if (any(thresh_sd)) labels[c(thresh_sd, thresh_sd)] <- c(paste0("-", round(thresh_vals[thresh_sd],2), " SD"),
                                                           paste0("+", round(thresh_vals[thresh_sd],2), " SD"))
  if (any(thresh_conf)) labels[c(c(thresh_conf, thresh_conf))] <- c(paste0("Lower ", round(thresh_vals[thresh_conf],2), "%CI"),
                                                                    paste0("Upper", round(thresh_vals[thresh_conf],2), "%CI"))

  break_order <- order(breaks)
  breaks <- breaks[break_order]
  labels <- labels[break_order]

  #Drop extant "dummy" tip
  tree <- treeio::drop.tip(tree, "Dummyextant")

  p <- unglue::unglue_data(names(tree@data), "rate{model}rlens{clock}_{summary}")
  rownames(p) <- names(tree@data)

  p <- p[rowSums(is.na(p)) < ncol(p),,drop=FALSE]

  summary <- match.arg(summary, c("mean", "median"))

  #Process lens
  if (all(p[["clock"]] == "")) {
    p[["clock"]] <- 1L
    if (!isTRUE(clock == 1)) {
      warning("Only one clock is available; ignoring 'clock'.", call. = FALSE)
    }
    clock <- 1L
  }
  else if (!is.numeric(clock) || !clock %in% as.numeric(p[["clock"]])) {
    stop(paste0("Only the following values are allowed for 'clock': ", paste(unique(p[["clock"]]), collapse = ", ")), call. = FALSE)
  }

  rate_var <- rownames(p)[as.numeric(p[["clock"]]) == clock & p[["summary"]] == summary]

  tree@data$`prob+-sd` <- as.factor(tree@data$`prob+-sd`)

  char_col <- vapply(tree@data, is.character, logical(1L))
  tree@data[char_col] <- lapply(tree@data[char_col], as.numeric)

  offset <- min(tree@data$age_median)

  if (is.null(xlim)) {
    x1 <- -round(max(tree@data$age_median) + 15, -1)
    x2 <- -round(min(tree@data$age_median) - 15, -1)
  }
  else {
    x1 <- round(xlim[1], -1)
    x2 <- round(xlim[2], -1)
  }


  #Create integer version of variable split up by breaks
  tree@data$clockfac <- as.numeric(cut(tree@data[[rate_var]], breaks = c(-Inf, breaks, Inf)))

  selection_plot <- ggtree::ggtree(tree, layout = "rectangular", ladderize = TRUE, right = TRUE,
                                                position = position_nudge(x = -offset),
                                                size = 2,
                                                aes(color = clockfac)) +
    ggtree::geom_tiplab(size = 3, linesize = 0.01, fontface = "italic", color="black", offset = -offset + .5) +
    scale_color_steps2("Background Rate\nThreshold",
                       low = low, high = high, mid = mid,
                       midpoint = mean(c(1, length(breaks) + 1)),
                       breaks = seq_along(breaks) + .5,
                       labels = labels,
                       limits = c(1.5 - 1e-8, length(breaks) + .5 + 1e-8)) +
    deeptime::coord_geo(xlim = c(x1, x2), ylim = c(0, treeio::Ntip(tree) + 2), expand = FALSE,
                        dat = list("epochs", "periods"), abbrv = list(TRUE, FALSE),
                        skip = geo_skip,
                        pos = list("bottom", "bottom"), alpha = 1, height = unit(1, "line"),
                        rot = 0, size = list(2.5, 3), neg = TRUE) +
    scale_x_continuous(n.breaks = 7, labels = abs) +
    ggtree::theme_tree2() +
    labs(title = "Selection Strength") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.position = c(.9, .2),
          legend.title = element_text(size = 10, face = "bold"))

  selection_plot <- ggtree::revts(selection_plot)

  return(selection_plot)

}
