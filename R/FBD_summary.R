#' Summarize FBD posterior parameter estimates
#'
#' Produces numerical summaries of each fossilized birth–death process (FBD)
#' posterior parameter by time bin.
#'
#' @inheritParams FBD_dens_plot
#' @param file An optional file path where the resulting table will be stored
#' using [write.csv()].
#' @param digits The number of digitis to round the summary results to. Default
#' is 3. See [round()].
#'
#' @returns
#' A data frame with a row for each parameter and time bin, and columns
#' for different summary statistics. These include the number of data points
#' (`n`) and the mean, standard deviation (`sd`), minimum value
#' (`min`), first quartile (`Q1`), median, third quartile
#' (`Q3`), and maximum value (`max`). When `file` is not
#' `NULL`, a .csv file containing this data frame will be saved to the
#' filepath specified in `file` and the output will be returned invisibly.
#'
#' @seealso
#' `vignette("fbd-params")` for the use of this function as part
#' of an analysis pipeline.
#'
#' [combine_log()] for producing a single data set of parameter
#' posterior samples from individual parameter log files.
#'
#' [FBD_reshape()] for converting posterior parameter table from wide
#' to long format.
#'
#' [FBD_dens_plot()], [FBD_normality_plot()],
#' [FBD_tests1()], and [FBD_tests2()] for other functions
#' used to summarize and display the distributions of the parameters.
#'
#' @examples
#' data("posterior3p")
#'
#' posterior3p_long <- FBD_reshape(posterior3p)
#'
#' FBD_summary(posterior3p_long)

#' @export
FBD_summary <- function(posterior, file = NULL, digits = 3) {

  if (!is.data.frame(posterior) || !"Time_bin" %in% names(posterior) || is.null(attr(posterior, "variables"))) {
    stop("'posterior' must be a data frame of MCMC posterior samples of FBD parameters.", call. = FALSE)
  }

  parameters <- attr(posterior, "variables")

  time.bins <- sort(unique(posterior$Time_bin))

  out <- expand.grid(parameter = parameters, Time_bin = time.bins,
                     stringsAsFactors = FALSE)

  summary.list <- lapply(seq_len(nrow(out)), function(i) {
    in_t <- which(posterior[["Time_bin"]] == out[["Time_bin"]][i])

    oneSummary(posterior[[out[["parameter"]][i]]][in_t], digits = digits)
  })

  out <- cbind(out, do.call("rbind", summary.list))

  out$parameter <- factor(out$parameter, levels = parameters)

  out <- out[with(out, order(parameter, Time_bin)),]

  rownames(out) <- NULL

  if (length(file) == 0L) {
    return(out)
  }

  write.csv(out, file = file)
  invisible(out)
}

#New summary function
oneSummary <- function(x, digits = 3) {
  qq <- unname(quantile(x, c(0, .25, .5, .75, 1)))
  d <- data.frame(
    n = length(x),
    mean = round(mean(x), digits),
    sd = round(sd(x), digits),
    min = round(qq[1], digits),
    Q1 = round(qq[2], digits),
    median = round(qq[3], digits),
    Q3 = round(qq[4], digits),
    max = round(qq[5], digits)
  )
  rownames(d) <- NULL

  d
}
