#All functions documented with examples

#Rate Table - BEAST2 version
#' Extract evolutionary rates from Bayesian clock trees produced by BEAST2
#'
#' BEAST2 stores the rates for each clock in a separate file. All trees need to be loaded using \code{treeio::read.beast}.
#'
#' @param ... \code{treedata} objects containing the summary trees with associated data on the rates for each separate clock.
#' @param summary summary metric used for the rates. Currently supported: \code{"mean"} or \code{"median"}, default \code{"median"}.
#' @param drop_dummy if not \code{NULL}, will drop the dummy extant tip with the given label from the BEAST2 summary trees prior to extracting the clock rates (when present). Default is \code{NULL}.
#'
#' @return A data frame with a column containing the node identifier (\code{node}) and one column containing the clock rates for each tree provided, in the same order as the trees.
#'
#' @seealso [get_clockrate_table_MrBayes()] for the equivalent function for MrBayes output files.
#' @seealso [clockrate_summary()] for summarizing and examining properties of the resulting rate table. Note that clade membership for each node must be customized (manually added) before these functions can be used, since this is tree and dataset dependent.
#'
#' @export
#'
#' @examples
#' #Import all clock summary trees produced by BEAST2 from your local directory
#' \dontrun{
#' tree_clock1 <- treeio::read.beast("tree_file_clock1.tre")
#' tree_clock2 <- treeio::read.beast("tree_file_clock2.tre")
#' }
#'
#' #Or use the example BEAST2 multiple clock trees that accompany EvoPhylo.
#' data(tree_clock1)
#' data(tree_clock2)
#'
#' # obtain the rate table from BEAST2 trees
#' rate_table <- get_clockrate_table_BEAST2(tree_clock1, tree_clock2, summary = "mean")
#'
#' @md
get_clockrate_table_BEAST2 <- function(..., summary = "median", drop_dummy = NULL) {
  trees <- list(...)
  if(length(trees) == 0) stop("No trees provided")

  summary <- match.arg(summary, c("mean", "median"))
  name <- if(summary == "mean") "rate" else "rate_median"

  if (!is.null(drop_dummy)) {
    trees <- lapply(trees, function(tr) treeio::drop.tip(tr, drop_dummy))
  }

  rate_table <- data.frame(nodes = as.integer(trees[[1]]@data$node))
  for(tr in trees) {
    data <- tr@data[match(rate_table$nodes, as.integer(tr@data$node)), ] #get rates in same order as 1st column
    rate_table <- cbind(rate_table, data[[name]])
  }
  colnames(rate_table)[2:ncol(rate_table)] <- paste0("rates", 1:(ncol(rate_table) - 1))
  row.names(rate_table) <- NULL

  rate_table
}

#Rate Table - MrBayes version
get_clockrate_table_MrBayes <- function(tree, summary = "median", drop_dummy = NULL) {

  if (!is.null(drop_dummy)) {
    tree <- treeio::drop.tip(tree, drop_dummy)
  }

  nodes <- as.integer(tree@data$node)

  p <- unglue::unglue_data(names(tree@data), "rate<model>rlens<clock>_<summary>",
                           open = "<", close = ">")
  rownames(p) <- names(tree@data)

  p <- p[rowSums(is.na(p)) < ncol(p),,drop=FALSE]

  p$clock <- gsub("\\{|\\}", "", p$clock)

  summary <- match.arg(summary, c("mean", "median"))

  rates <- rownames(p)[p$summary == summary]

  rate_table <- setNames(data.frame(nodes, tree@data[rates]),
                        c("nodes", paste0("rates", p[rates, "clock"])))

  for (i in seq_len(ncol(rate_table))[-1]) rate_table[[i]] <- as.numeric(rate_table[[i]])

  return(rate_table)
}

#Summary stats for clades
clockrate_summary <- function(rate_table, file = NULL, digits = 3) {

  if (!is.data.frame(rate_table)) {
    stop("'rate_table' must be a data frame.", call. = FALSE)
  }
  if (!any(startsWith(names(rate_table), "rates"))) {
    stop("'rate_table' must contain \"rates\" columns containing clockrate summaries.", call. = FALSE)
  }
  if (!hasName(rate_table, "clade")) {
    stop("A 'clade' column must be present in the data.", call. = FALSE)
  }

  clocks <- sort(gsub("rates", "", names(rate_table)[startsWith(names(rate_table), "rates")], fixed = TRUE))

  clades <- sort(unique(rate_table$clade))
  if ("Other" %in% clades) clades <- c(setdiff(clades, "Other"), "Other")

  out <- do.call("rbind", lapply(clades, function(cl) {
    in_cl <- which(rate_table$clade == cl)
    cbind(data.frame(clade = cl,
                     clock = clocks),
          do.call("rbind", lapply(clocks, function(r) {
            oneSummary(rate_table[[paste0("rates", r)]][in_cl], digits = digits)
          }))
    )
  }))

  out$clade <- factor(out$clade, levels = clades)

  out <- out[with(out, order(clock, clade)),]

  rownames(out) <- NULL

  if (all(clocks == "")) out$clock <- NULL

  if (length(file) > 0) {
    write.csv(out, file = file)
    invisible(out)
  }
  else {
    return(out)
  }
}

#Density plot of rates by clade
clockrate_dens_plot <- function(rate_table, clock = NULL, stack = FALSE, nrow = 1, scales = "fixed") {

  if (!is.data.frame(rate_table)) {
    stop("'rate_table' must be a data frame.", call. = FALSE)
  }
  if (!any(startsWith(names(rate_table), "rates"))) {
    stop("'rate_table' must contain \"rates\" columns containing clockrate summaries.", call. = FALSE)
  }
  if (!hasName(rate_table, "clade")) {
    stop("A 'clade' column must be present in the data.", call. = FALSE)
  }

  clock_cols <- which(startsWith(names(rate_table), "rates"))

  if (is.null(clock)) {
    rate_table <- rate_table[c("clade", "nodes", names(rate_table)[clock_cols])]
  }
  else {
    if (!is.numeric(clock) || anyNA(clock) ||
        !all(clock %in% as.numeric(gsub("rates", "", names(rate_table)[clock_cols])))) {
      stop(paste0("'clock' must be a vector of clock indices. In this data, ",
                  ngettext(length(clock_cols), "1 clock is ",
                           paste(length(clock_cols), "clocks are")),
                  " available."),
           call. = FALSE)
    }
    rate_table <- rate_table[c("clade", "nodes", paste0("rates", clock))]
  }

  clades <- sort(unique(rate_table$clade))
  if ("Other" %in% clades) clades <- c(setdiff(clades, "Other"), "Other")
  rate_table$clade <- factor(rate_table$clade, levels = clades)

  rt <- clock_reshape(rate_table)
  levels(rt$clock) <- paste("Clock", levels(rt$clock))

  rateplot <- ggplot(data = rt, mapping = aes(x = .data$rate, fill = .data$clade,
                                              color = .data$clade)) +
    geom_hline(yintercept = 0) +
    geom_density(position = if (stack) "stack" else "identity",
                 alpha = if (stack) 1 else .3) +
    scale_x_continuous() +
    labs(x = "Rate", y = "Density", fill = "Clade", color = "Clade") +
    theme_bw() +
    theme(legend.position = "top")

  # if (nlevels(rt$clock) > 1) {
    rateplot <- rateplot  + facet_wrap(~.data$clock, nrow = nrow, scales = scales)
  # }
  rateplot
}

#Plot linear model and Pearson correlation of one rate against another
clockrate_reg_plot <- function(rate_table, clock_x, clock_y, method = "lm", show_lm = TRUE, ...) {

  if (!is.data.frame(rate_table)) {
    stop("'rate_table' must be a data frame.", call. = FALSE)
  }
  if (!any(startsWith(names(rate_table), "rates"))) {
    stop("'rate_table' must contain \"rates\" columns containing clockrate summaries.", call. = FALSE)
  }

  clock_cols <- which(startsWith(names(rate_table), "rates"))

  if (length(clock_cols) <= 1) {
    stop("At least two clock rates must be present in the input.", call. = FALSE)
  }

  clocks <- sort(gsub("rates", "", names(rate_table)[clock_cols], fixed = TRUE))

  if (missing(clock_x) && missing(clock_y)) {
    clock_x <- clocks[1]
    clock_y <- clocks[2]
  }
  else if (missing(clock_x)) {
    if (length(clock_y) != 1 || !paste0("rates", clock_y) %in% names(rate_table)[clock_cols]) {
      stop("'clock_y' must be the index of a clock rate in the input.", call. = FALSE)
    }
    clock_x <- setdiff(clocks, as.character(clock_y))[1]
  }
  else if (missing(clock_y)) {
    if (length(clock_x) != 1 || !paste0("rates", clock_x) %in% names(rate_table)[clock_cols]) {
      stop("'clock_x' must be the index of a clock rate in the input.", call. = FALSE)
    }
    clock_y <- setdiff(clocks, as.character(clock_x))[1]
  }
  else {
    if (length(clock_x) != 1 || length(clock_y) != 1 ||
        !all(paste0("rates", c(clock_x, clock_y)) %in% names(rate_table)[clock_cols])) {
      stop("'clock_x' and 'clock_y' must be the indices of clock rates in the input.", call. = FALSE)
    }
  }

  names(rate_table)[names(rate_table) == paste0("rates", clock_x)] <- "clock_x"
  names(rate_table)[names(rate_table) == paste0("rates", clock_y)] <- "clock_y"

  regplot <- ggplot(rate_table, aes(x = clock_x, y = clock_y)) +
    geom_point() +
    geom_smooth(method = method, formula = y ~ x, ...) +
    scale_x_continuous() +
    scale_y_continuous() +
    labs(x = paste("Clock", clock_x), y = paste("Clock", clock_y)) +
    theme_bw()

  if (show_lm) {
    r <- cor(rate_table$clock_x, rate_table$clock_y)

    #Extract underlying ggplot data to place correlation in correct place in plot
    ggbd <- ggplot_build(regplot)$data

    ggbd1 <- ggbd[[1]] #geom_point data
    ggbd2 <- ggbd[[2]] #geom_smooth data

    min_x <- min(min(ggbd1$x), min(ggbd2$x))
    max_x <- max(max(ggbd1$x), max(ggbd2$x))

    min_y <- min(min(ggbd1$y), min(ggbd2$y),
                 if (hasName(ggbd2, "ymin")) min(ggbd2$ymin)) #FALSE when se = FALSE
    max_y <- max(max(ggbd1$y), max(ggbd2$y),
                 if (hasName(ggbd2, "ymax")) max(ggbd2$ymax)) #FALSE when se = FALSE

    regplot <- regplot +
      annotate("text", label = c(paste0("italic(R)^2 == ", round(r^2, 2)),
                                  paste0("italic(r) == ", round(r, 2))), parse = TRUE,
               x = .3*min_x + .7*max_x,
               y = c(.85*min_y + .15*max_y,
                     .8*min_y + .2*max_y))
  }

  regplot
}


