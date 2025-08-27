#' Combine and filter (.p) log files from Mr.Bayes
#'
#' Imports parameter (.p) log files from Mr. Bayes and combines them into a
#' single data frame. Samples can be dropped from the start of each log file
#' (i.e., discarded as burn-in) and/or downsampled to reduce the size of the
#' output object.
#'
#' @details
#' `combine_log()` imports log files produced by Mr.Bayes, ignoring the
#' first row of the file (which contains an ID number). The files are appended
#' together, optionally after removing burn-in generations from the beginning
#' and/or by further filtering throughout the rest of each file. When
#' `burnin` is greater than 0, the number or propotion of generations
#' corresponding to the supplied value will be dropped from the beginning of
#' each file as it is read in. For example, setting `burnin = .25` (the
#' default) will drop the first 25\% of generations from each file. When
#' `downsample` is greater than 0, the file will be downsampled until the
#' number or proportion of generations corresponding to the supplied value is
#' reached. For example, if `downsample = 10000` generations (the default)
#' for log files from 4 independent runs (i.e., 4 (.p) files), each log file
#' will be downsampled to 2500 generations, and the final combined data frame
#' will contain 10000 samples, selected in approximately equally spaced
#' intervals from the original data.
#'
#' The output can be supplied to [get_pwt_rates_MrBayes()] and to
#' [FBD_reshape()]. The latter will convert the log data frame from
#' my wide to long format, which is necessary to be used as input for
#' downstream analyses using [FBD_summary()],
#' [FBD_dens_plot()], [FBD_normality_plot()],
#' [FBD_tests1()], or [FBD_tests2()].
#'
#' @param path The path to a folder containing (.p) log files or a character
#' vector of log files to be read.
#' @param burnin Either the number or a proportion of generations to drop from
#' the beginning of each log file.
#' @param downsample Either the number or the proportion of generations the
#' user wants to keep after downsampling for the final (combined) log file.
#' Generations will be dropped in approximately equally-spaced intervals.
#'
#' @returns
#' A data frame with columns corresponding to the columns in the
#' supplied log files and rows containing the sampled parameter values.
#' Examples of the kind of output produced can be accessed using
#' [`data("posterior1p")`][posterior1p] and
#' [`data("posterior3p")`][posterior3p].
#'
#' @seealso
#' `vignette("fbd-params")` for the use of this function as part
#' of an analysis pipeline.
#'
#' [FBD_reshape()], which reshapes a combined parameter log file for
#' use in some other package functions.
#'
#' @examples
#' \dontrun{
#' posterior <- combine_log("path/to/folder", burnin = .25,
#'                          downsample = 10000)
#' }
#'

#' @export
combine_log <- function(path = ".", burnin = .25, downsample = 1e4) {
  #Get FBD parameter estimates from collection of log files (.p)
  files <- NULL
  if (is.character(path)) {
    if (all(utils::file_test(path, op = "-f"))) {
      files <- path
    }
    else if (length(path) == 1L && utils::file_test(path, op = "-d")) {
      files <- paste0(path, "/", list.files(path, pattern = '\\.p$'))
    }
  }

  if (length(files) == 0L) {
    stop("The value supplied to 'path' must be a character vector containing the names of parameter log (.p) files or of a folder containg such files.", call. = FALSE)
  }

  if (length(burnin) != 1 || !is.numeric(burnin) || burnin < 0) {
    stop("'burnin' must be a single number corresponding to the number or percentage of rows to drop from the beginning of each log file.",
         call. = FALSE)
  }

  if (length(downsample) != 1L || !is.numeric(downsample) || downsample < 0) {
    stop("'downsample' must be a single number corresponding to the number or percentage of rows to remain after downsampling from each log file.",
         call. = FALSE)
  }

  L <- lapply(files, function(x) {
    rtest <- read.table(x, skip = 1, header = TRUE, nrows = 3)
    if (!all(c("Gen", "LnL", "LnPr") %in% names(rtest))) {
      return(NULL)
    }

    r <- read.table(x, skip = 1, header = TRUE)

    if (burnin <= 0) {
      return(r)
    }

    if (burnin < 1) {
      b <- seq_len(round(burnin * NROW(r)))
    }
    else {
      b <- seq_len(round(min(burnin, NROW(r))))
    }

    r[-b, , drop = FALSE]
  })

  if (all(lengths(L) == 0L)) {
    stop("No parameter log files were found in the supplied path.", call. = FALSE)
  }

  L[lengths(L) == 0L] <- NULL

  if (!all(vapply(L, function(i) identical(names(i), names(L[[1L]])), logical(1L)))) {
    stop("All parameter log files must have the same column names.", call. = FALSE)
  }

  samples <- do.call("rbind", L)

  if (downsample > 0) {
    if (downsample < 1) {
      d <- as.integer(round(seq(1, NROW(samples), length.out = NROW(samples) * downsample)))
    }
    else {
      d <- as.integer(round(seq(1, NROW(samples), length.out = min(NROW(samples), downsample))))
    }

    samples <- samples[d,,drop = FALSE]
  }

  rownames(samples) <- NULL

  return(samples)
}
