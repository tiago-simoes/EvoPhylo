#' Convert an FBD posterior parameter table from wide to long format
#'
#' Converts FBD posterior parameter table, such as those imported using
#' [combine_log()], from wide to long format.
#'
#' @details
#' The posterior parameters log files produced by Bayesian evolutionary
#' analyses using skyline birth-death tree models, including the skyline FBD
#' model, result into two or more estimates for each FBD parameter, one for
#' each time bin. This function will convert a table of parameters with skyline
#' FBD parameters from wide to long format, with one row per generation per
#' time bin and a new column "Time_bin" containing the respective time bins as
#' a factor. The long format is necessary for downstream analyses using
#' [FBD_summary()], [FBD_dens_plot()],
#' [FBD_normality_plot()], [FBD_tests1()], or
#' [FBD_tests2()], as similarly done by [clock_reshape()]
#' for clock rate tables.
#'
#' The format of the log files can either be specified using the
#' `variables` and `log.type` or auto-detected by the function. The
#' "posterior" data frame can be obtained by reading in a log file directly
#' (e.g. using the `read.table` function) or by combining several output
#' log files from Mr. Bayes using [combine_log()].
#'
#' @param posterior Single posterior parameter sample dataset with skyline FBD
#' parameters produced with [combine_log()].
#' @param variables Names of FBD rate variables in the log. If `NULL` (default),
#' will attempt to auto-detect the names and log type.
#' @param log.type Name of the software which produced the log (currently
#' supported: MrBayes or BEAST2). Has to be set if `variables` is not
#' `NULL`.
#'
#' @returns
#' A data frame of posterior parameter estimates containing a single
#' "Time_bin" column and one column for each FBD parameter value.
#'
#' @seealso
#' `vignette("fbd-params")` for the use of this function as part
#' of an analysis pipeline.
#'
#' [combine_log()], [reshape()]
#'
#' @examples
#' data("posterior3p")
#'
#' head(posterior3p)
#'
#' ## Reshape FBD table to long format
#' posterior3p_long <- FBD_reshape(posterior3p)
#'
#' head(posterior3p_long)
#'

#' @export
FBD_reshape <- function(posterior, variables = NULL, log.type = c("MrBayes", "BEAST2")) {
  if (!is.data.frame(posterior)) {
    stop("'posterior' must be a data frame.", call. = FALSE)
  }

  if (is.null(variables)) {
    autodetect <- detect_posterior(posterior)
    variables <- autodetect$variables
    log.type <- autodetect$log.type
  }
  else {
    exist <- vapply(variables, function(nm) {
      any(startsWith(names(posterior), nm))
    }, logical(1L))

    if (!all(exist)) {
      stop("Specified variables not found in posterior")
    }

    if (length(log.type) > 1 || !log.type %in% c("MrBayes", "BEAST2")) {
      stop("Log type must be one of 'MrBayes' or 'BEAST2'")
    }
  }

  varying <- lapply(variables, function(v) {
    names(posterior)[startsWith(names(posterior), v)]
  })

  idname <- if (log.type == "MrBayes") "Gen" else "Sample"

  posterior_long <- reshape(posterior, direction = "long",
                            varying = varying,
                            v.names = variables,
                            timevar = "Time_bin",
                            sep = if(log.type == "MrBayes") "_" else ".",
                            idvar = "Gen", ids = posterior[[idname]])

  posterior_long[["Time_bin"]] <- factor(posterior_long[["Time_bin"]])
  rownames(posterior_long) <- NULL
  attr(posterior_long, "reshapeLong") <- NULL

  attr(posterior_long, "log.type") <- log.type
  attr(posterior_long, "variables") <- variables

  posterior_long
}
