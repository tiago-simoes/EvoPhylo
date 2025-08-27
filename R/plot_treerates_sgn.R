#' Plot Bayesian evolutionary tree with rate thresholds for selection mode
#'
#' Plots the summary Bayesian evolutionary tree with branches, according to
#' user-defined thresholds (in units of standard deviations) used to infer the
#' strength and mode of selection.
#'
#' @details
#' Plots the phylogentic tree contained in `tree` using
#' [ggtree::ggtree()]. Branches undergoing
#' accelerating evolutionary rates (e.g., `"1 SD"`, `"3 SD"`, or
#' `"5 SD"` relative to the background rate) for each morphological clock
#' partition suggest directional (or positive) selection for that morphological
#' partition in that branch of the tree. Branches undergoing decelerating
#' evolutionary rates (e.g., `"1 SD"`, `"3 SD"`, or `"5 SD"`
#' relative to the background rate) for each morphological clock partition
#' suggest stabilizing selection for that morphological partition in that
#' branch of the tree. For details on rationale, see Simões & Pierce (2021).
#'
#' Please double check that the distribution of background rates (mean rates
#' for the tree) sampled from the posterior follow the assumptions of a normal
#' distribution (e.g., check for normality of distribution in Tracer).
#' Otherwise, displayed results may not have a valid interpretation.
#'
#' @inheritParams plot_back_rates
#' @param tree A `tidytree` object; the output of a call to
#' [treeio::read.beast()]. Summary trees from Mr.
#' Bayes will include branch specific rates for all clock partitions, and the
#' partition to be plotted will be specified using the "clock" argument. On the
#' other hand, BEAST2 will output one separate summary tree file for each clock
#' partition. For the latter, the tree file for the partition of interest
#' should be provided for plotting.
#' @param summary Only when using Mr. Bayes trees. The rate summary stats
#' chosen to calculate selection mode. Only rates "mean" and "median" are
#' allowed. Default is "mean".
#' @param drop.dummyextant `logical`; Only when using Mr. Bayes trees.
#' Whether to drop the "Dummyextant" tip (if present) from the tree before
#' plotting the tree. Default is `TRUE`.
#' @param threshold A vector of threshold values. Default is to display thresholds of ±1 relative standard deviation (SD) of the relative posterior clock rates. Should be specified as a number of standard deviations (e.g., `"1 SD"`) or the confidence level for a confidence internal around the mean relative posterior clockrate (e.g., `"95\%"`). Multiple values are allowed to produce a plot with multiple thresholds. Set to `NULL` to omit thresholds.
#' @param low,mid,high Colors passed to [ggplot2::scale_color_steps2()] to control the colors of the branches based on which thresholds are exceeded. When no thresholds are supplied, use `mid` to control the color of the tree.
#' @param branch_size The thickness of the lines that form the tree.
#' @param tip_size The font size for the tips of the tree.
#' @param xlim The x-axis limits. Should be two negative numbers (though the
#' axis labels will be in absolute value, i.e., Ma).
#' @param nbreaks The number of interval breaks in the geological timescale.
#' @param geo_size The font size for the labels in the geological scale. The
#' first value in `list()` is the font size for geological epochs and the
#' second value is for geological periods. Passed directly to the `size`
#' argument of [deeptime::coord_geo()].
#' @param geo_skip A vector of interval names indicating which intervals should
#' not be labeled. Passed directly to the `skip` argument of
#' [deeptime::coord_geo()].
#'
#' @returns
#' A `ggtree` object, which inherits from `ggplot`.
#'
#' @seealso
#' `vignette("rates-selection")` for the use of this function as
#' part of an analysis pipeline.
#'
#' [ggtree::ggtree()], [deeptime::coord_geo()]
#'
#' @references
#' Simões, T. R. and S. E. Pierce (2021). Sustained High Rates of
#' Morphological Evolution During the Rise of Tetrapods. *Nature Ecology & Evolution* 5: 1403–1414.
#'
#' @examples
#' ## MrBayes example
#' # Load example tree and posterior
#' data("tree3p")
#' data("posterior3p")
#'
#' plot_treerates_sgn(
#'   type = "MrBayes",
#'   tree3p, posterior3p,          #MrBayes tree file with data for all partitions
#'   trans = "none",
#'   summary = "mean",             #MrBayes specific argument
#'   drop.dummyextant = TRUE,      #MrBayes specific argument
#'   clock = 1,                           #Show rates for clock partition 1
#'   threshold = c("1 SD", "3 SD"),       #sets background rate threshold for selection mode
#'   branch_size = 1.5, tip_size = 3,                          #sets size for tree elements
#'   xlim = c(-450, -260), nbreaks = 8, geo_size = list(3, 3)) #sets limits and breaks for geoscale
#'
#' \dontrun{
#' ## BEAST2 example
#' tree_clock1 <- system.file("extdata", "Penguins_MCC_morpho_part1", package = "EvoPhylo")
#' tree_clock1 <- treeio::read.beast(tree_clock1)
#' posterior <- system.file("extdata", "Penguins_log.log", package = "EvoPhylo")
#' posterior <- read.table(posterior, header = TRUE)
#'
#' plot_treerates_sgn(
#'   type = "BEAST2",
#'   tree_clock1, posterior,                 #BEAST2 tree file with data for partition 1
#'   trans = "log10",
#'   clock = 1,                              #Show rates for clock partition 1
#'   threshold = c("1 SD", "3 SD"),          #sets background rate threshold for selection mode
#'   branch_size = 1.5, tip_size = 3,                        #sets size for tree elements
#'   xlim = c(-70, 30), nbreaks = 8, geo_size = list(3, 3))  #sets limits and breaks for geoscale
#' }

#' @export
plot_treerates_sgn <- function(type = c("MrBayes", "BEAST2"),
                               tree, posterior,
                               trans = c("none", "log", "log10"),
                               summary = "mean", drop.dummyextant = TRUE,
                               clock = 1, threshold = c("1 SD", "2 SD"),
                               low = "blue", mid = "gray90", high = "red",
                               branch_size = 2, tip_size = 2,
                               xlim = NULL, nbreaks = 10, geo_size = list(2, 3),
                               geo_skip = c("Quaternary", "Holocene", "Late Pleistocene")) {


  if (!is.character(type) || length(type) != 1L || !type %in% c("MrBayes", "BEAST2")) {
    stop("Bad type call")
  }

  #Drop extant "dummy" tip
  if (type == "MrBayes" && drop.dummyextant) {
    tree <- treeio::drop.tip(tree, "Dummyextant")
  }

  #Process threshold
  if (length(threshold) > 0) {
    if (type == "BEAST2" && (missing(posterior) || !is.data.frame(posterior))) {
      stop("'posterior' must be a data frame.", call. = FALSE)
    }

    if (type == "MrBayes") {
      if (missing(posterior) || !is.data.frame(posterior)) {
        stop("'posterior' must be a data frame.", call. = FALSE)
      }
      if (hasName(posterior, "clockrate.all.")) {
        names(posterior)[which(names(posterior) == "clockrate.all.")] <- "clockrate"
      }
      if (!hasName(posterior, "clockrate")) {
        stop("A 'clockrate' column must be present in 'posterior'.", call. = FALSE)
      }
    }

    if (!is.character(threshold)) {
      stop("'threshold' must be a character vector.", call. = FALSE)
    }

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

    if (type == "BEAST2") {
      #get BEAST2 relative background clock rate (for the desired clock partition) and data transform
      posterior.clockrate<-get_clockrate_posterior(posterior)                     #get rate table from posterior sample using 'get_clockrate_posterior' helper
      posterior.clockrate.long<-posterior_clockrate_reshape(posterior.clockrate)  #convert posterior mean rates table to long using 'posterior_clockrate_reshape' helper
      posterior.final<-posterior.clockrate.long[posterior.clockrate.long$clock == clock,]                       #keep values for desired clock

      if (trans == "none") {
        #get relative background clock rate
        posterior.rel.clockrate <- posterior.final$rates.post / mean(posterior.final$rates.post)
      }
      else if (trans == "log10") {
        posterior.final$rates.post.log <- log10(posterior.final$rates.post)
        #get relative background clock rate
        posterior.rel.clockrate <- posterior.final$rates.post.log/mean(posterior.final$rates.post.log)
      }
      else {
        # ln transform data
        posterior.final$rates.post.log<-log(posterior.final$rates.post)
        #get relative background clock rate
        posterior.rel.clockrate <- posterior.final$rates.post.log/mean(posterior.final$rates.post.log)
      }
    }
    else {
      #get Mr. Bayes relative background clock rate (shared among all partitions)
      if (trans == "none") {
        posterior.rel.clockrate <- posterior$clockrate / mean(posterior$clockrate)
      }
      else if (trans == "log10"){
        posterior$clockrate.log <- log10(posterior$clockrate)
        #get relative background clock rate
        posterior.rel.clockrate <- posterior$clockrate.log / mean(posterior$clockrate.log)
      }
      else {
        # ln transform data
        posterior$clockrate.log <- log(posterior$clockrate)
        #get relative background clock rate
        posterior.rel.clockrate <- posterior$clockrate.log / mean(posterior$clockrate.log)
      }
    }

    mean.posterior.rel.clockrate <- 1

    breaks <- numeric(2 * length(threshold))
    labels <- character(2 * length(threshold))

    # Get threshold sd/conf values for each clock partition
    if (any(thresh_sd)) {
      breaks[c(thresh_sd, thresh_sd)] <- mean.posterior.rel.clockrate + c(-thresh_vals[thresh_sd], thresh_vals[thresh_sd]) * sd(posterior.rel.clockrate)
      labels[c(thresh_sd, thresh_sd)] <- c(sprintf("-%s SD", round(thresh_vals[thresh_sd], 2)),
                                           sprintf("+%s SD", round(thresh_vals[thresh_sd], 2)))
    }

    if (any(thresh_conf)) {
      n <- length(posterior.rel.clockrate)
      tcrits <- qt(.5 * (1 + thresh_vals[thresh_conf] / 100), n - 1)
      breaks[c(thresh_conf, thresh_conf)] <- mean.posterior.rel.clockrate + c(-tcrits, tcrits) * sd(posterior.rel.clockrate)/sqrt(n)
      labels[c(thresh_conf, thresh_conf)] <- c(sprintf("Lower %s%%CI", round(thresh_vals[thresh_conf],2)),
                                               sprintf("Lower %s%%CI", round(thresh_vals[thresh_conf],2)))
    }

    break_order <- order(breaks)
    breaks <- breaks[break_order]
    labels <- labels[break_order]
  }
  else {
    breaks <- numeric()
    labels <- character()
  }

  if (type == "MrBayes") {
    #Getting multiple clock rates
    p <- unglue::unglue_data(names(tree@data), "rate<model>Brlens<clock>_<summary>",
                             open = "<", close = ">")
    rownames(p) <- names(tree@data)

    p <- p[rowSums(is.na(p)) < ncol(p),,drop=FALSE]

    p$clock <- gsub("\\{|\\}", "", p$clock)

    summary <- match.arg(summary, c("mean", "median"))

    #Process lens
    if (all(p[["clock"]] == "") || all(p[["clock"]] == "all")) {
      if (!isTRUE(clock == 1)) {
        warning("Only one clock is available; ignoring 'clock' value.", call. = FALSE)
      }

      if (all(p[["clock"]] == "all")) {
        warning("All data partitions were analyzed as a single clock ('all'); ignoring value in argument 'clock'.", call. = FALSE)
      }

      p[["clock"]] <- 1L
      clock <- 1L
    }
    else if (!is.numeric(clock) || !clock %in% as.numeric(p[["clock"]])) {
      stop(paste0("All 'clock' values must be numeric, but the following values were found in the tree file: ", toString(unique(p[["clock"]]))), call. = FALSE)
    }

    rate_var <- rownames(p)[as.numeric(p[["clock"]]) == clock & p[["summary"]] == summary]

    tree@data$`prob+-sd` <- as.factor(tree@data$`prob+-sd`)

    char_col <- vapply(tree@data, is.character, logical(1L))
    tree@data[char_col] <- lapply(tree@data[char_col], as.numeric)

    offset <- min(tree@data$age_median)
  }
  else {
    offset <- min(tree@data$height_median)
  }

  if (type == "MrBayes" && is.null(xlim)) {
    x1 <- -round(max(tree@data$age_median) + 15, -1)
    x2 <- -round(min(tree@data$age_median) - 15, -1)
  }
  else if (type == "BEAST2" && is.null(xlim)) {
    x1 <- -round(max(tree@data$height_median) + 15, -1)
    x2 <- -round(min(tree@data$height_median) - 15, -1)
  }
  else {
    x1 <- round(xlim[1], -1)
    x2 <- round(xlim[2], -1)
  }


  #Create integer version of rate variable split up by breaks
  #MrBayes already uses relative branch rates (normalized)
  if (type == "BEAST2") {
    #get relative branch rates (normalize) and split up by breaks
    tree@data$rel.rate <- tree@data$rate / mean(posterior.final$rates.post)
    tree@data$clockfac <- as.numeric(cut(tree@data$rel.rate, breaks = c(-Inf, breaks, Inf)))
  }
  else {
    tree@data$clockfac <- as.numeric(cut(tree@data[[rate_var]], breaks = c(-Inf, breaks, Inf)))
  }

  #Make tree plot
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
    labs(title = sprintf("Selection for partition %s", clock), call. = FALSE) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.position = c(.05, .25),
          legend.title = element_text(size = 8, face = "bold"),
          legend.key.size = unit(0.5,'cm'),
          legend.text = element_text(size=7))

  ggtree::revts(selection_plot)
}
