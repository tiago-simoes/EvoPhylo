#Rate Table
get_clockrate_table <- function(tree, summary = "mean", file = NULL) {

  tree <- treeio::drop.tip(tree, "Dummyextant")

  nodes <- as.integer(tree@data$node)

  p <- unglue::unglue_data(names(tree@data), "rate{model}rlens{clock}_{summary}")
  rownames(p) <- names(tree@data)

  p <- p[rowSums(is.na(p)) < ncol(p),,drop=FALSE]

  summary <- match.arg(summary, c("mean", "median"))

  rates <- rownames(p)[p$summary == summary]

  RateTable <- setNames(data.frame(nodes, tree@data[rates]),
                        c("nodes", paste0("rates", p[rates, "clock"])))

  for (i in names(RateTable)[-1]) RateTable[[i]] <- as.numeric(RateTable[[i]])

  if (length(file) > 0) {
    write.csv(RateTable, file = file)
    invisible(RateTable)
  }
  else {
    return(RateTable)
  }
}

#Summary stats for clades
clockrate_summary <- function(rate_table, file = NULL, digits = 3) {

  if (!"clade" %in% names(rate_table)) {
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
  rownames(out) <- NULL

  out <- out[with(out, order(clock, clade)),]

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

  rate_table_long <- reshape(rate_table, direction = "long",
                             idvar = c("clade", "nodes"),
                             v.names = "rates",
                             varying = startsWith(names(rate_table), "rates"),
                             timevar = "clock")
  rate_table_long$clock <- factor(rate_table_long$clock)
  levels(rate_table_long$clock) <- paste("Clock", levels(rate_table_long$clock))


  rateplot <- ggplot(rate_table_long, aes(x = rates, fill = clade, color = clade)) +
    geom_density(position = if (stack) "stack" else "identity",
                 alpha = if (stack) 1 else .3) +
    scale_x_continuous() +
    labs(x = "Rate", y = "Density", fill = "Clade", color = "Clade") +
    theme_bw() +
    theme(legend.position = "top")

  if (nlevels(rate_table_long$clock) > 1) {
    rateplot <- rateplot  + facet_wrap(vars(clock), nrow = nrow, scales = scales)
  }
  rateplot
}

#Regression plot of one rate against another
clockrate_reg_plot <- function(rate_table, clock_x, clock_y, method = "lm", corr = TRUE) {

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
    geom_smooth(method = method, se = TRUE, formula = y ~ x) +
    scale_x_continuous() +
    scale_y_continuous() +
    labs(x = paste("Clock", clock_x), y = paste("Clock", clock_y)) +
    theme_bw()

  if (corr) {
    r <- cor(rate_table$clock_x, rate_table$clock_y)

    ggbd1 <- ggplot_build(regplot)$data[[1]]
    ggbd2 <- ggplot_build(regplot)$data[[2]]

    min_x <- min(min(ggbd1$x), min(ggbd2$x))
    max_x <- max(max(ggbd1$x), max(ggbd2$x))

    min_y <- min(min(ggbd1$y), min(ggbd2$y), min(ggbd2$ymin))
    max_y <- max(max(ggbd1$y), max(ggbd2$y), max(ggbd2$ymax))

    regplot <- regplot +
      annotate("label", label = paste0("italic(r) == ", round(r, 2)), parse = TRUE,
               x = .3*min_x + .7*max_x,
               y = .9*min_y + .1*max_y)
  }

  regplot
}

