#' Compute rate summary statistics across clades and clocks
#'
#' Computes summary statistics for each clade and/or each clock partition. The
#' input must have a "clade" column.
#'
#' @details
#' The user must manually add clades to the rate table produced by
#' [get_clockrate_table_MrBayes()] before it can be used with this
#' function. This can be doen manually within R, such as by using a graphical
#' user interface for editing data like the \pkg{DataEditR} package, or by
#' writing the rate table to a spreadsheet and reading it back in after adding
#' the clades. The example below uses a table that has had the clades added.
#'
#' @param rate_table A data frame of clock rates, such as from the output of
#' [get_clockrate_table_MrBayes()] with an extra `"clade"`
#' column.
#' @param file An optional file path where the resulting table will be stored
#' using [write.csv()].
#' @param digits The number of digits to round the summary results to. Default
#' is 3. See [round()].
#'
#' @returns
#' A data frame containing a row for each clade and each clock with
#' summary statistics (n, mean, standard deviation, minimum, 1st quartile,
#' median, third quartile, maximum).
#'
#' @seealso
#' `vignette("rates-selection")` for the use of this function as
#' part of an analysis pipeline.
#'
#' [get_clockrate_table_MrBayes()], [summary()]
#'
#' @examples
#' data("RateTable_Means_3p_Clades")
#'
#' clockrate_summary(RateTable_Means_3p_Clades)
#'

#' @export
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

  if (all(clocks == "")) {
    out$clock <- NULL
  }

  if (length(file) == 0L) {
    return(out)
  }

  write.csv(out, file = file)
  invisible(out)
}
