#' Plot clock rate distributions
#'
#' Plots the distribution density of clock rates by clock and clade. The input
#' must have a "clade" column.
#'
#' @details
#' The user must manually add clades to the rate table produced by
#' [get_clockrate_table_MrBayes()] before it can be used with this
#' function. This can be done manually with in R, such as by using a graphical
#' user interface for editing data like the \pkg{DataEditR} package, or by
#' writing the rate table to a spreadsheet and reading it back in after adding
#' the clades. The example below uses a table that has had the clades added.
#'
#' @param rate_table A data frame of clock rates, such as from the output of
#' [get_clockrate_table_MrBayes()] with an extra "clade" column.
#' @param clock Which clock rates will be plotted. If unspecified, all clocks
#' are plotted.
#' @param stack Whether to display stacked density plots (`TRUE`) or
#' overlapping density plots (`FALSE`).
#' @param nrow When plotting rates for more than one clock, how many rows
#' should be filled by the plots. This is passed to [ggplot2::facet_wrap()].
#' @param scales When plotting rates for more than one clock, whether the axis
#' scales should be "fixed" (default) across clocks or allowed to vary ("free",
#' "free_x", or "free_y"). This is passed to [ggplot2::facet_wrap()].
#'
#' @returns
#' A `ggplot` object, which can be modified using \pkg{ggplot2}
#' functions.
#'
#' @seealso
#' `vignette("rates-selection")` for the use of this function as
#' part of an analysis pipeline.
#'
#' [get_clockrate_table_MrBayes()], [ggplot2::geom_density()]
#'
#' @examples
#' data("RateTable_Means_3p_Clades")
#'
#' # Overlapping plots
#' clockrate_dens_plot(RateTable_Means_3p_Clades, stack = FALSE,
#'                     nrow = 1, scales = "fixed")
#'
#' # Stacked density for all three clocks, changing the color
#' # palette to viridis using ggplot2 functions
#' clockrate_dens_plot(RateTable_Means_3p_Clades,
#'                     clock = 1:3, nrow = 1, stack = TRUE,
#'                     scales = "fixed") +
#'   ggplot2::scale_color_viridis_d() +
#'   ggplot2::scale_fill_viridis_d()
#'

#' @export
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

  rateplot <- ggplot(data = rt, mapping = aes(x = .data$rates, fill = .data$clade,
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
