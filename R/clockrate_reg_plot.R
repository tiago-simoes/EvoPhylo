#' Plot regression lines between sets of rates
#'
#' Displays a scatterplot and fits regression line of one set of clock rates
#' against another, optionally displaying their Pearson correlation coefficient
#' (r) and R-squared values (R^2).
#'
#' @details
#' `clockrate_reg_plot()` can only be used when multiple clocks are
#' present in the clock rate table. Unlike [clockrate_summary()] and
#' [clockrate_dens_plot()], no "clade" column is required.
#'
#' @param rate_table A table of clock rates, such as from the output of
#' [get_clockrate_table_MrBayes()].
#' @param clock_x,clock_y The clock rates that should go on the x- and y-axes,
#' respectively.
#' @param method The method (function) used fit the regression of one clock on
#' the other. Check the `method` argument in [ggplot2::geom_smooth()] for all options. Default
#' is `"lm"` for a linear regression model. `"glm"` and
#' `"loess"` are alternative options.
#' @param show_lm Whether to display the Pearson correlation coefficient (r)
#' and R-squared values (R^2) between two sets of clock rates.
#' @param \dots Other arguments passed to [ggplot2::geom_smooth()].
#'
#' @returns
#' A `ggplot` object, which can be modified using \pkg{ggplot2}
#' functions.
#'
#' @seealso
#' `vignette("rates-selection")` for the use of this function as
#' part of an analysis pipeline.
#'
#' [ggplot2::geom_point()], [ggplot2::geom_smooth()]
#'
#' @examples
#' data("RateTable_Means_3p_Clades")
#'
#' #Plot correlations between clocks 1 and 3
#' clockrate_reg_plot(RateTable_Means_3p_Clades,
#'                    clock_x = 1, clock_y = 3)
#'
#' #Use arguments supplied to geom_smooth():
#' clockrate_reg_plot(RateTable_Means_3p_Clades,
#'                    clock_x = 1, clock_y = 3,
#'                    color = "red", se = FALSE)
#'

#' @export
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

  regplot <- ggplot(rate_table, aes(x = .data$clock_x, y = .data$clock_y)) +
    geom_point() +
    geom_smooth(method = method, formula = y ~ x, ...) +
    scale_x_continuous() +
    scale_y_continuous() +
    labs(x = paste("Clock", clock_x),
         y = paste("Clock", clock_y)) +
    theme_bw()

  if (!show_lm) {
    return(regplot)
  }

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

  regplot +
    annotate("text", label = c(paste0("italic(R)^2 == ", round(r^2, 2)),
                               paste0("italic(r) == ", round(r, 2))), parse = TRUE,
             x = .3*min_x + .7*max_x,
             y = c(.85*min_y + .15*max_y,
                   .8*min_y + .2*max_y))
}
