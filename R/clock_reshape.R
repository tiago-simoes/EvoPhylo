#' Convert clock rate tables from wide to long format
#'
#' Converts clock rate tables, such as those produced by
#' [clockrate_summary()] and imported back after including clade
#' names, from wide to long format.
#'
#' @details
#' This function will convert clock rate tables from wide to long format, with
#' a new column "clock" containing the clock partition from where each rate
#' estimate was obtained as a factor. The long format is necessary for
#' downstream analyses of selection strength (mode), as similarly done by
#' [FBD_reshape()] for posterior parameter log files.
#'
#' @param rate_table A data frame of clock rates, such as from the output of
#' [get_clockrate_table_MrBayes()] with an extra "clade" column.
#'
#' @returns
#' A data frame containing a single "value" column (for all rate
#' values) and one column for the "clock" variable (indicating to which clock
#' partition each rate values refers to).
#'
#' @seealso
#' `vignette("rates-selection")` for the use of this function as
#' part of an analysis pipeline.
#'
#' [get_clockrate_table_MrBayes()], [summary()],
#' [clockrate_summary()], [FBD_reshape()]
#'
#' @examples
#' ## The example dataset rate_table_clades_means3
#' ## has clades and 3 clock rate columns:
#' data("rate_table_clades_means3")
#'
#' ## Reshape a clock rate table with clade names to long format
#' \dontrun{
#' rates_by_clade <- clock_reshape(rate_table_clades_means3)
#' }

#' @export
clock_reshape <- function(rate_table) {
  rate_table_long <- reshape(rate_table, direction = "long",
                             varying = which(startsWith(names(rate_table), "rates")),
                             v.names = "rates",
                             timevar = "clock",
                             idvar = "nodes",
                             sep = "_")
  rate_table_long[["clock"]] <- factor(rate_table_long[["clock"]])
  rownames(rate_table_long) <- NULL
  attr(rate_table_long, "reshapeLong") <- NULL

  rate_table_long
}
