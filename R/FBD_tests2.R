#' Test for differences in FBD parameter values
#'
#' `FBD_tests2()` performs t-tests and Mann-Whitney U-tests to compare the
#' average value of fossilized birth–death process (FBD) parameters between
#' time bins.
#'
#' @details
#' [pairwise.t.test()] and [pairwise.wilcox.test()] are
#' used to calculate, respectively, the t-test and Mann-Whitney U-tests
#' statistics and p-values. Because the power of these tests depends on the
#' number of posterior samples, it can be helpful to examine the distributions
#' of FBD parameter posteriors using [FBD_dens_plot()] instead of
#' relying heavily on the tests.
#'
#' @inheritParams FBD_tests1
#' @param p.adjust.method The method use to adjust the p-values for multiple
#' testing. See [p.adjust()] for details and options. Default if
#' `"fdr"` for the Benjamini-Hochberg false discovery rate correction.
#'
#' @returns
#' A list with an element for each test, each of which contains a list
#' of test results for each parameter. The results are in the form of a data
#' frame containing the sample sizes and unadjusted and adjusted p-values for
#' each comparison.
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
#' used to summarize and display the distributions of the parameter posteriors.
#'
#' [pairwise.t.test()] and [pairwise.wilcox.test()] for the
#' tests used.
#'
#' @examples
#' data("posterior3p")
#'
#' posterior3p_long <- FBD_reshape(posterior3p)
#'
#' FBD_tests2(posterior3p_long)

#' @export
FBD_tests2 <- function(posterior, p.adjust.method = "fdr") {

  if (!is.data.frame(posterior) || !"Time_bin" %in% names(posterior) || is.null(attr(posterior, "variables"))) {
    stop("'posterior' must be a data frame of MCMC posterior samples of FBD parameters.", call. = FALSE)
  }

  parameters <- attr(posterior, "variables")

  posterior$Time_bin <- factor(posterior$Time_bin)

  df_dimnames <- list(NULL, c("parameter", "Time_bin1", "Time_bin2", "n1", "n2", "p-value", "p-value adj"))
  n_tests <- choose(nlevels(posterior[["Time_bin"]]), 2)

  d <- as.data.frame(matrix(nrow = n_tests, ncol = length(df_dimnames[[2]]),
                            dimnames = df_dimnames))

  #T-tests between pairs of Time_bins for each param
  t_test_list <- setNames(lapply(parameters, function(p) {
    t <- pairwise.t.test(posterior[[p]],
                         posterior[["Time_bin"]],
                         p.adjust.method = "none")

    d$parameter <- p

    d[paste0("Time_bin", 1:2)] <- as.data.frame(t(combn(levels(posterior[["Time_bin"]]), 2)))

    for (i in seq_len(nrow(d))) {
      d[["n1"]][i] <- sum(posterior[["Time_bin"]] == d[["Time_bin1"]][i])
      d[["n2"]][i] <- sum(posterior[["Time_bin"]] == d[["Time_bin2"]][i])
      d[["p-value"]][i] <- t$p.value[d[["Time_bin2"]][i], d[["Time_bin1"]][i]]
    }

    d[["p-value adj"]] <- p.adjust(d[["p-value"]], p.adjust.method)
    d
  }), parameters)

  #Mann-Whitney U tests between pairs of Time_bins for each param
  mwu_test_list <- setNames(lapply(parameters, function(p) {

    t <- pairwise.wilcox.test(posterior[[p]],
                              posterior[["Time_bin"]],
                              p.adjust.method = "none")

    d$parameter <- p

    d[paste0("Time_bin", 1:2)] <- as.data.frame(t(combn(levels(posterior[["Time_bin"]]), 2)))

    for (i in seq_len(nrow(d))) {
      d[["n1"]][i] <- sum(posterior[["Time_bin"]] == d[["Time_bin1"]][i])
      d[["n2"]][i] <- sum(posterior[["Time_bin"]] == d[["Time_bin2"]][i])
      d[["p-value"]][i] <- t$p.value[d[["Time_bin2"]][i], d[["Time_bin1"]][i]]
    }


    d[["p-value adj"]] <- p.adjust(d[["p-value"]], p.adjust.method)
    d
  }), parameters)

  list(t_tests = t_test_list,
       mwu_tests = mwu_test_list)
}
