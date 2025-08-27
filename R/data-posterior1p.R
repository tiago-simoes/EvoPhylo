#' Posterior parameter samples (single clock)
#'
#' An example dataset of posterior parameter samples resulting from a
#' clock-based Bayesian inference analysis using the skyline fossilized
#' birth–death process (FBD) tree model with Mr. Bayes after combining all
#' parameter (.p) files into a single data frame with the
#' [combine_log()] function. This particular example was produced by
#' analyzing the data set with a single morphological partition from Simões &
#' Pierce (2021).
#'
#' Datasets like this one can be produced from parameter log (.p) files using
#' [combine_log()]. The number of variables depends on parameter set
#' up, but for clock analyses with Mr. Bayes, will typically include the ones
#' above, possibly also including `alpha`, which contains the shape of the
#' gamma distribution governing how much rates vary across characters. When
#' using the traditional FBD model rather than the skyline FBD model used to
#' produce this dataset, there will be only one column for each of
#' `net_speciation`, `relative_extinction` and
#' `relative_fossilization`. When using more than one morphological
#' partition, different columns may be present; see [`posterior3p`]
#' for an example with 3 partitions.
#'
#' @name posterior1p
#' @docType data
#' @usage data("posterior1p")
#' @format A data frame with 4000 observations on several variables estimated
#' for each generation during analysis:
#' \describe{
#' \item{Gen}{A numeric vector
#' for the generation number}
#' \item{`LnL`}{A numeric vector for the
#' natural log likelihood of the cold chain}
#' \item{`LnPr`}{A numeric
#' vector for the natural log likelihood of the priors}
#' \item{`TH`}{A
#' numeric vector for the total tree height (sum of all branch durations, as
#' chronological units)}
#' \item{`TL`}{A numeric vector for total tree
#' length (sum of all branch lengths, as accumulated substitutions/changes)}
#' \item{`prop_ancfossil`}{A numeric vector indicating the proportion of
#' fossils recovered as ancestors} \item{`sigma`}{A numeric vector for
#' the standard deviation of the lognormal distribution governing how much
#' rates vary across characters.} \item{`net_speciation_1`}{A numeric
#' vector for net speciation estimates for each time bin}
#' \item{`net_speciation_2`}{A numeric vector for net speciation
#' estimates for each time bin}
#' \item{`net_speciation_3`}{A numeric vector
#' for net speciation estimates for each time bin}
#' \item{`net_speciation_4`}{A numeric vector for net speciation
#' estimates for each time bin}
#' \item{`relative_extinction_1`}{A numeric
#' vector for relative extinction estimates for each time bin}
#' \item{`relative_extinction_2`}{A numeric vector for relative
#' extinction estimates for each time bin}
#' \item{`relative_extinction_3`}{A numeric vector for relative
#' extinction estimates for each time bin}
#' \item{`relative_extinction_4`}{A numeric vector for relative
#' extinction estimates for each time bin}
#' \item{`relative_fossilization_1`}{A numeric vector for relative
#' fossilization estimates for each time bin}
#' \item{`relative_fossilization_2`}{A numeric vector for relative
#' fossilization estimates for each time bin}
#' \item{`relative_fossilization_3`}{A numeric vector for relative
#' fossilization estimates for each time bin}
#' \item{`relative_fossilization_4`}{A numeric vector for relative
#' fossilization estimates for each time bin}
#' \item{`tk02var`}{A numeric
#' vector for the variance on the base of the clock rate}
#' \item{`clockrate`}{A numeric vector for the base of the clock rate}
#' }
#' @seealso [`posterior3p`] for an example dataset of posterior parameter samples resulting from an analysis with 3 partitions rather than 1.
#'
#' @references
#' Simões, T. R. and S. E. Pierce (2021). Sustained High Rates of
#' Morphological Evolution During the Rise of Tetrapods. *Nature Ecology & Evolution* 5: 1403–1414.
#' @keywords datasets
NULL
