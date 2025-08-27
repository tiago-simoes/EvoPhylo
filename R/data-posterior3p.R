#' Posterior parameter samples (3 clock partions)
#'
#' An example dataset of posterior parameter samples resulting from a
#' clock-based Bayesian inference analysis using the skyline fossilized
#' birth–death process (FBD) tree model with Mr. Bayes after combining all
#' parameter (.p) files into a single data frame with the
#' [combine_log()] function. This particular example was produced by
#' analyzing the data set with three morphological partitions from Simões &
#' Pierce (2021).
#'
#' Datasets like this one can be produced from parameter log (.p) files using
#' [combine_log()]. The number of variables depends on parameter set
#' up, but for clock analyses with Mr. Bayes, will typically include the ones
#' above, possibly also including an `alpha` for each partition, which
#' contains the shape of the gamma distribution governing how much rates vary
#' across characters (when shape of the distribution is unlinked across
#' partitions). When using the traditional FBD model rather than the skyline
#' FBD model used to produce this dataset, there will be only one column for
#' each of `net_speciation`, `relative_extinction` and
#' `relative_fossilization`. When using a single morphological partition,
#' different columns may be present; see [`posterior1p`] for an
#' example with just one partition.
#'
#' @name posterior3p
#' @docType data
#' @usage data("posterior3p")
#' @format A data frame with 4000 observations on several variables estimated
#' for each generation during analysis. The number of variables depends on
#' parameter set up, but for clock analyses with Mr. Bayes, will typically
#' include the following:
#' \describe{
#' \item{`Gen`}{A numeric vector for
#' the generation number}
#' \item{`LnL`}{A numeric vector for the natural
#' log likelihood of the cold chain}
#' \item{`LnPr`}{A numeric vector for
#' the natural log likelihood of the priors}
#' \item{`TH.all.`}{A numeric
#' vector for the total tree height (sum of all branch durations, as
#' chronological units)}
#' \item{`TL.all.`}{A numeric vector for total tree
#' length (sum of all branch lengths, as accumulated substitutions/changes)}
#' \item{`prop_ancfossil.all.`}{A numeric vector indicating the
#' proportion of fossils recovered as ancestors}
#' \item{`sigma.1.`}{A
#' numeric vector for the standard deviation of the lognormal distribution
#' governing how much rates vary across characters for each data
#' partition}
#' \item{`sigma.2.`}{A numeric vector for the
#' standard deviation of the lognormal distribution governing how much rates
#' vary across characters for each data partition}
#' \item{`sigma.3.`}{A numeric vector for the standard
#' deviation of the lognormal distribution governing how much rates vary across
#' characters for each data partition}
#' \item{`m.1.`}{A numeric vector for
#' the rate multiplier parameter for each data partition}
#' \item{`m.2.`}{A numeric vector for the rate multiplier
#' parameter for each data partition}
#' \item{`m.3.`}{A numeric
#' vector for the rate multiplier parameter for each data partition}
#' \item{`net_speciation_1.all.`}{A numeric vector for net speciation
#' estimates for each time bin}
#' \item{`net_speciation_2.all.`}{A numeric
#' vector for net speciation estimates for each time bin}
#' \item{`net_speciation_3.all.`}{A numeric vector for net speciation
#' estimates for each time bin}
#' \item{`net_speciation_4.all.`}{A numeric
#' vector for net speciation estimates for each time bin}
#' \item{`relative_extinction_1.all.`}{A numeric vector for relative
#' extinction estimates for each time bin}
#' \item{`relative_extinction_2.all.`}{A numeric vector for relative
#' extinction estimates for each time bin}
#' \item{`relative_extinction_3.all.`}{A numeric vector for relative
#' extinction estimates for each time bin}
#' \item{`relative_extinction_4.all.`}{A numeric vector for relative
#' extinction estimates for each time bin}
#' \item{`relative_fossilization_1.all.`}{A numeric vector for relative
#' fossilization estimates for each time bin}
#' \item{`relative_fossilization_2.all.`}{A numeric vector for
#' relative fossilization estimates for each time bin}
#' \item{`relative_fossilization_3.all.`}{A numeric vector for
#' relative fossilization estimates for each time bin}
#' \item{`relative_fossilization_4.all.`}{A numeric vector for
#' relative fossilization estimates for each time bin}
#' \item{`tk02var.1.`}{A numeric vector for the variance on the base of
#' the clock rate for each clock partition}
#' \item{`tk02var.2.`}{A numeric vector for the variance on the
#' base of the clock rate for each clock partition}
#' \item{`tk02var.3.`}{A numeric vector for the variance on the
#' base of the clock rate for each clock partition}
#' \item{`clockrate.all.`}{A numeric vector for the base of the clock
#' rate}
#' }
#' @seealso [`posterior1p`] for an example dataset of posterior
#' parameter samples resulting from an analysis with 1 partition rather than 3.
#' @references Simões, T. R. and S. E. Pierce (2021). Sustained High Rates of
#' Morphological Evolution During the Rise of Tetrapods. *Nature Ecology &
#' Evolution* 5: 1403–1414.
#' @keywords datasets
NULL
