\name{posterior1p}
\alias{posterior1p}
\docType{data}
\title{
Posterior parameter samples (single clock)
}
\description{
An example dataset of posterior parameter samples resulting from a clock-based Bayesian inference analysis using the skyline fossilized birth–death process (FBD) tree model with Mr. Bayes after combining all parameter (.p) files into a single data frame with the \code{\link{combine_log}} function. This particular example was produced by analyzing the data set with a single morphological partition from Simões & Pierce (2021).

}
\usage{data("posterior1p")}
\format{
  A data frame with 4000 observations on several variables estimated for each generation during analysis:
  \describe{
    \item{Gen}{A numeric vector for the generation number}
    \item{\code{LnL}}{A numeric vector for the natural log likelihood of the cold chain}
    \item{\code{LnPr}}{A numeric vector for the natural log likelihood of the priors}
    \item{\code{TH}}{A numeric vector for the total tree height (sum of all branch durations, as chronological units)}
    \item{\code{TL}}{A numeric vector for total tree length (sum of all branch lengths, as accumulated substitutions/changes)}
    \item{\code{prop_ancfossil}}{A numeric vector indicating the proportion of fossils recovered as ancestors}
    \item{\code{sigma}}{A numeric vector for the standard deviation of the lognormal distribution governing how much rates vary across characters.}
    \item{\code{net_speciation_1}, \code{net_speciation_2}, \code{net_speciation_3}, \code{net_speciation_4}}{A numeric vector for net speciation estimates for each time bin}
    \item{\code{relative_extinction_1}, \code{relative_extinction_2}, \code{relative_extinction_3}, \code{relative_extinction_4}}{A numeric vector for relative extinction estimates for each time bin}
    \item{\code{relative_fossilization_1}, \code{relative_fossilization_2}, \code{relative_fossilization_3}, \code{relative_fossilization_4}}{A numeric vector for relative fossilization estimates for each time bin}
    \item{\code{tk02var}}{A numeric vector for the variance on the base of the clock rate}
    \item{\code{clockrate}}{A numeric vector for the base of the clock rate}
   }
}
\details{
Datasets like this one can be produced from parameter log (.p) files using \code{\link{combine_log}}. The number of variables depends on parameter set up, but for clock analyses with Mr. Bayes, will typically include the ones above, possibly also including \code{alpha}, which contains the shape of the gamma distribution governing how much rates vary across characters. When using the traditional FBD model rather than the skyline FBD model used to produce this dataset, there will be only one column for each of \code{net_speciation}, \code{relative_extinction} and \code{relative_fossilization}. When using more than one morphological partition, different columns may be present; see \code{\link{posterior3p}} for an example with 3 partitions.
}
\seealso{
\code{\link{posterior3p}} for an example dataset of posterior parameter samples resulting from an analysis with 3 partitions rather than 1.
}
\references{
Simões, T. R. and S. E. Pierce (2021). Sustained High Rates of Morphological Evolution During the Rise of Tetrapods. \emph{Nature Ecology & Evolution} 5: 1403–1414.
}
\keyword{datasets}



