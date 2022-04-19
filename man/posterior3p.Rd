\name{posterior3p}
\alias{posterior3p}
\docType{data}
\title{
Posterior parameter samples (3 clock partion analysis)
}
\description{
An example dataset of posterior parameter samples resulting from a clock-based Bayesian inference analysis using the skyline fossilized birth–death process (FBD) tree model with Mr. Bayes after combining all parameter (.p) files into a single data frame with the \code{\link{combine_log}} function. This particular example was produced by analyzing the data set with three morphological partitions from Simões & Pierce (2021).

}
\usage{data("posterior3p")}
\format{
  A data frame with 4000 observations on several variables estimated for each generation during analysis. The number of variables depends on parameter set up, but for clock analyses with Mr. Bayes, will typically include the following:
  \describe{
    \item{\code{Gen}}{A numeric vector for the generation number}
    \item{\code{LnL}}{A numeric vector for the natural log likelihood of the cold chain}
    \item{\code{LnPr}}{A numeric vector for the natural log likelihood of the priors}
    \item{\code{TH}}{A numeric vector for the total tree height (sum of all branch durations, as chronological units)}
    \item{\code{TL}}{A numeric vector for total tree length (sum of all branch lengths, as accumulated substitutions/changes)}
    \item{\code{prop_ancfossil}}{A numeric vector indicating the proportion of fossils recovered as ancestors (when using the SA tree model)}
    \item{\code{alpha.i}}{A numeric vector for the shape of the gamma distribution governing how much rates vary across characters for the\item{\code{i}}data partition (when shape of the distribution is unlinked across partitions).}
    \item{\code{sigma.i}}{A numeric vector for the standard deviation of the lognormal distribution governing how much rates vary across characters for the\item{\code{i}}data partition (when shape of the distribution is unlinked across partitions).}
    \item{\code{m.i}}{A numeric vector for the rate multiplier parameter for the \item{\code{i}} data partition}
    \item{\code{tk02var.i}}{A numeric vector for the variance on the base of the clock rate for the \item{\code{i}} clock partition}
    \item{\code{clockrate}}{A numeric vector for the base of the clock rate}
    \item{\code{net_speciation.i}}{A numeric vector for net speciation estimates for each time bin (a single value when using the traditional FBD model and multiple \item{\code{i}} values when using the skyline FBD model)}
    \item{\code{relative_extinction}}{A numeric vector for relative extinction estimates for each time bin (a single value when using the traditional FBD model and multiple \item{\code{i}} values when using the skyline FBD model)}
    \item{\code{relative_fossilization}} {A numeric vector for relative fossilization estimates for each time bin (a single value when using the traditional FBD model and multiple \item{\code{i}} values when using the skyline FBD model)}
  }
}
\seealso{
\code{\link{posterior3p_long}}
\code{\link{posterior1p}}
\code{\link{posterior1p_long}}
}
\keyword{datasets}

\references{
Simões, T. R. and S. E. Pierce (2021). Sustained High Rates of Morphological Evolution During the Rise of Tetrapods. \emph{Nature Ecology & Evolution} 5: 1403–1414.
}

