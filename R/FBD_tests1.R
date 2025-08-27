#' Test assumptions of normality and homoscedasticity for FBD posterior
#' parameters
#'
#' Produces tests of normality (within time bin, ignoring time bin, and pooling
#' within-time bin values) and homoscedasticity (homogeneity of variances) for
#' each fossilized birth–death process (FBD) parameter in the posterior
#' parameter log file.
#'
#' @details
#' `FBD_tests1()` performs several tests on the posterior distributions of
#' parameter values within and across time bins. It produces the Shapiro-Wilk
#' test for normality using [shapiro.test()] and the Bartlett and
#' Fligner tests for homogeneity of variance using [bartlett.test()]
#' and [fligner.test()], respectively. Note that these tests are
#' likely to be significant even if the observations are approximately normally
#' distributed or have approximately equal variance; therefore, they should be
#' supplemented with visual inspection using [FBD_normality_plot()].
#'
#' @inheritParams FBD_dens_plot
#' @param downsample Whether to downsample the observations to ensure
#' Shapiro-Wilk normality tests can be run. If `TRUE`, observations will
#' be dropped so that no more than 5000 observations are used for the tests on
#' the full dataset, as required by [shapiro.test()]. They will be
#' dropped in evenly spaced intervals. If `FALSE` and there are more than
#' 5000 observations for any test, that test will not be run.
#'
#' @returns
#' A list containing the results of the three tests with the following
#' elements:
#' \item{shapiro}{A list with an element for each parameter. Each
#' element is a data frame with a row for each time bin and the test statistic
#' and p-value for the Shapiro-Wilk test for normality. In addition, there will
#' be a row for an overall test, combining all observations ignoring time bin,
#' and a test of the residuals, which combines the group-mean-centered
#' observations (equivalent to the residuals in a regression of the parameter
#' on time bin).}
#' \item{bartlett}{A data frame of the Bartlett test for
#' homogeneity of variance across time bins with a row for each parameter and
#' the test statistic and p-value for the test.}
#' \item{fligner}{A data frame of
#' the Fligner test for homogeneity of variance across time bins with a row for
#' each parameter and the test statistic and p-value for the test.}
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
#' [FBD_normality_plot()] for visual assessments.
#'
#' [FBD_tests2()] for tests of differences between parameter means.
#'
#' [shapiro.test()], [bartlett.test()], and
#' [fligner.test()] for the statistical tests used.
#'
#' @examples
#' data("posterior3p")
#'
#' posterior3p_long <- FBD_reshape(posterior3p)
#'
#' FBD_tests1(posterior3p_long)

#' @export
FBD_tests1 <- function(posterior, downsample = TRUE) {
  if (!is.data.frame(posterior) || !"Time_bin" %in% names(posterior) || is.null(attr(posterior, "variables"))) {
    stop("'posterior' must be a data frame of MCMC posterior samples of FBD parameters.", call. = FALSE)
  }

  parameters <- attr(posterior, "variables")

  posterior$Time_bin <- factor(posterior$Time_bin)

  #Bartlett test for homogeneity of variance across Time_bins
  #Run for each param
  d <- as.data.frame(matrix(nrow = 1, ncol = 3,
                            dimnames = list(NULL, c("parameter", "statistic", "p-value"))))

  bartlett_df <- do.call("rbind", lapply(parameters, function(p) {
    d[["parameter"]] <- p
    bt <- bartlett.test(posterior[[p]], posterior$Time_bin)
    d[["statistic"]] <- bt$statistic
    d[["p-value"]] <- bt$p.value
    d
  }))

  #Fligner-Killeen test for homogeneity of variance across Time_bins
  #Run for each param
  fligner_df <- do.call("rbind", lapply(parameters, function(p) {
    d[["parameter"]] <- p
    bt <- fligner.test(posterior[[p]], posterior$Time_bin)
    d[["statistic"]] <- bt$statistic
    d[["p-value"]] <- bt$p.value
    d
  }))

  #Shapiro-Wilk test for normality
  #Run within each Time_bin group, overall, and on residuals for each param
  #Need to downsample for SW test
  max.n <- floor(5000 / nlevels(posterior$Time_bin))

  if (downsample) {
    keep <- unlist(lapply(levels(posterior$Time_bin),
                          function(i) {
                            if (sum(posterior$Time_bin == i) > max.n) {
                              which(posterior$Time_bin == i)[round(seq(1, sum(posterior$Time_bin == i), length.out = max.n))]
                            }
                            else {
                              which(posterior$Time_bin == i)
                            }
                          }))
    posterior <- posterior[keep, , drop = FALSE]
    run.sw.test <- rep(TRUE, nlevels(posterior$Time_bin))
  }
  else {
    t <- table(posterior$Time_bin)

    run.sw.test <- t[levels(posterior$Time_bin)] <= max.n

    if (!any(run.sw.test)) {
      warning("Shapiro-Wilk normality tests require downsampling and will not be run. Set downsample = TRUE to run these tests.", call. = FALSE)
    }
  }

  d <- as.data.frame(matrix(nrow = nlevels(posterior$Time_bin) + 2, ncol = 3,
                            dimnames = list(c(paste("Time bin", levels(posterior$Time_bin)), "Overall", "Residuals"),
                                            c("parameter", "statistic", "p-value"))))
  shapiro_df <- setNames(lapply(parameters, function(p) {
    d[["parameter"]] <- p

    for (i in seq_len(nlevels(posterior$Time_bin))[run.sw.test]) {
      sw <- shapiro.test(posterior[[p]][posterior$Time_bin == levels(posterior$Time_bin)[i]])
      d[["statistic"]][i] <- sw[["statistic"]]
      d[["p-value"]][i] <- sw[["p.value"]]
    }

    if (all(run.sw.test)) {
      #Overall
      sw <- shapiro.test(posterior[[p]])
      d[["statistic"]][nlevels(posterior$Time_bin) + 1] <- sw[["statistic"]]
      d[["p-value"]][nlevels(posterior$Time_bin) + 1] <- sw[["p.value"]]

      #Residuals
      means <- ave(posterior[[p]], posterior$Time_bin)
      sw <- shapiro.test(posterior[[p]] - means)
      d[["statistic"]][nlevels(posterior$Time_bin) + 2] <- sw[["statistic"]]
      d[["p-value"]][nlevels(posterior$Time_bin) + 2] <- sw[["p.value"]]
    }

    d
  }), parameters)

  list(shapiro = shapiro_df,
       bartlett = bartlett_df,
       fligner = fligner_df)
}
