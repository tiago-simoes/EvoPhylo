#Reshape rate_table with a column for each rate to long
clock_reshape <- function(rate_table) {
  rate_table_long <- reshape(rate_table, direction = "long",
                             varying = which(startsWith(names(rate_table), "rates")),
                             v.names = "rate",
                             timevar = "clock",
                             idvar = "nodes",
                             sep = "_")
  rate_table_long[["clock"]] <- factor(rate_table_long[["clock"]])
  rownames(rate_table_long) <- NULL
  attr(rate_table_long, "reshapeLong") <- NULL

  rate_table_long
}


#t-tests
#rate_table = rate_table_means
#posterior.clockrate = samples$clockrate?
get_pwt_rates <- function(rate_table, posterior) {

  if (missing(rate_table) || !is.data.frame(rate_table)) {
    stop("'rate_table' must be a data frame.", call. = FALSE)
  }
  if (!any(startsWith(names(rate_table), "rates"))) {
    stop("'rate_table' must contain \"rates\" columns containing clockrate summaries.", call. = FALSE)
  }
  if (!hasName(rate_table, "clade")) {
    stop("A 'clade' column must be present in 'rate_table'.", call. = FALSE)
  }
  if (missing(posterior) || !is.data.frame(posterior)) {
    stop("'posterior' must be a data frame.", call. = FALSE)
  }
  if (!hasName(posterior, "clockrate.all.")) {
    stop("A 'clockrate.all.' column must be present in 'posterior'.", call. = FALSE)
  }

  posterior.clockrate <- posterior$clockrate.all.
  post.df <- length(posterior.clockrate) - 1
  post.mean <- mean(posterior.clockrate)

  rate_table_long <- clock_reshape(rate_table)
  rate_table_long$abs_rate <- rate_table_long$rate * post.mean

  post.se <- sd(posterior.clockrate)/sqrt(length(posterior.clockrate))
  post.ts <- abs(post.mean - rate_table_long$abs_rate)/post.se
  pvals <- 2*pt(post.ts, df = post.df, lower.tail = FALSE)

  out <- data.frame(rate_table_long$clade,
                    rate_table_long$nodes,
                    rate_table_long$clock,
                    rate_table_long$rate,
                    rate_table_long$abs_rate,
                    rate_table_long$abs_rate,
                    pvals)
  names(out) <- c("clade", "nodes", "clock", "relative rate", "absolute rate (mean)", "null", "p.value")
  out
}

#Plot tree with colored thresholds
plot_treerates_sgn <- function(tree, posterior, clock = 1, summary = "mean", threshold = c("1 SD", "2 SD"),
                               drop.dummyextant = TRUE,
                               low = "blue", mid = "gray90", high = "red",
                               branch_size = 2, tip_size = 2,
                               xlim = NULL, nbreaks = 10, geo_size=list(2, 3),
                               geo_skip = c("Quaternary", "Holocene", "Late Pleistocene")) {

  #Process threshold
  if (length(threshold) > 0) {
    if (missing(posterior) || !is.data.frame(posterior)) {
      stop("'posterior' must be a data frame.", call. = FALSE)
    }
    if (!hasName(posterior, "clockrate.all.")) {
      stop("A 'clockrate.all.' column must be present in 'posterior'.", call. = FALSE)
    }

    if (!is.character(threshold)) stop("'threshold' must be a character vector.", call. = FALSE)
    thresh_conf <- endsWith(threshold, "%")
    thresh_sd <- endsWith(tolower(threshold), "sd")

    if (any(!thresh_conf & !thresh_sd)) {
      stop("All entries in 'threshold' must end in '%' for confidence intervals or 'SD' for standard deviations.", call. = FALSE)
    }

    thresh_vals <- character(length(threshold))

    if (any(thresh_conf)) thresh_vals[thresh_conf] <- substring(threshold[thresh_conf], 1, nchar(threshold[thresh_conf]) - 1)
    if (any(thresh_sd)) thresh_vals[thresh_sd] <- substring(threshold[thresh_sd], 1, nchar(threshold[thresh_sd]) - 2)

    if (anyNA(suppressWarnings(as.numeric(thresh_vals)))) {
      stop("All entries in 'threshold' must be confidence levels or the number of standard deviations to use as the thresholds.", call. = FALSE)
    }
    thresh_vals <- as.numeric(thresh_vals)

    #Use relative clockrate
    posterior.rel.clockrate <- posterior$clockrate.all./mean(posterior$clockrate.all.)
    mean.posterior.rel.clockrate <- 1

    breaks <- numeric(2*length(threshold))
    labels <- character(2*length(threshold))

    if (any(thresh_sd)) {
      breaks[c(thresh_sd, thresh_sd)] <- mean.posterior.rel.clockrate + c(-thresh_vals[thresh_sd], thresh_vals[thresh_sd]) * sd(posterior.rel.clockrate)
      labels[c(thresh_sd, thresh_sd)] <- c(sprintf("-%s SD", round(thresh_vals[thresh_sd], 2)),
                                           sprintf("+%s SD", round(thresh_vals[thresh_sd], 2)))
    }
    if (any(thresh_conf)) {
      n <- length(posterior.rel.clockrate)
      tcrits <- qt(.5*(1 + thresh_vals[thresh_conf]/100), n - 1)
      breaks[c(thresh_conf, thresh_conf)] <- mean.posterior.rel.clockrate + c(-tcrits, tcrits) * sd(posterior.rel.clockrate)/sqrt(n)
      labels[c(thresh_conf, thresh_conf)] <- c(sprintf("Lower %s%%CI", round(thresh_vals[thresh_conf],2)),
                                               sprintf("Lower %s%%CI", round(thresh_vals[thresh_conf],2)))
    }

    break_order <- order(breaks)
    breaks <- breaks[break_order]
    labels <- labels[break_order]
  }
  else {
    breaks <- numeric(0)
    labels <- character(0)
  }

  #Drop extant "dummy" tip
  if (drop.dummyextant) {
    tree <- treeio::drop.tip(tree, "Dummyextant")
  }

  p <- unglue::unglue_data(names(tree@data), "rate<model>rlens<clock>_<summary>",
                           open = "<", close = ">")
  rownames(p) <- names(tree@data)

  p <- p[rowSums(is.na(p)) < ncol(p),,drop=FALSE]

  p$clock <- gsub("\\{|\\}", "", p$clock)

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
                                   size = branch_size,
                                   mapping = aes(color = .data$clockfac)) +
    ggtree::geom_tiplab(size = tip_size, linesize = 0.01, fontface = "italic",
                        color = "black", offset = -offset + .5) +
    scale_color_steps2("Background Rate\nThreshold",
                       low = low, high = high,
                       mid = mid,
                       midpoint = mean(c(1, length(breaks) + 1)),
                       breaks = seq_along(breaks) + .5,
                       labels = labels,
                       limits = c(1.5 - 1e-8, length(breaks) + .5 + 1e-8)) +
    deeptime::coord_geo(xlim = c(x1, x2), ylim = c(-1, treeio::Ntip(tree) + 2), expand = FALSE,
                        dat = list("epochs", "periods"), abbrv = list(TRUE, FALSE),
                        skip = geo_skip,
                        pos = list("bottom", "bottom"), alpha = 1, height = unit(1, "line"),
                        rot = 0, size = geo_size, neg = TRUE) +
    scale_x_continuous(n.breaks = nbreaks, labels = abs) +
    ggtree::theme_tree2() +
    labs(title = "Selection Strength") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.position = c(.05, .25),
          legend.title = element_text(size = 8, face = "bold"),
          legend.key.size = unit(0.5,'cm'),
          legend.text = element_text(size=7))

  selection_plot <- ggtree::revts(selection_plot)

  return(selection_plot)
}
