\name{FBD_normality_plot}
\alias{FBD_normality_plot}

\title{
Inspect FBD parameter distributions visually
}
\description{
\code{FBD_normality_plot()} produces plots of the distributions of fossilized birth–death process (FBD) parameters (net speciation, relative exinction, and relative fossilization) to facilitate the assessment of the assumptions of normality within time bins and homogeneity of variance across time bins.
}
\usage{
FBD_normality_plot(posterior)
}

\arguments{
  \item{posterior}{
a dataset of FBD parameter MCMC posterior samples imported using \code{\link{import_log}}. See the \code{\link{posterior1p}} dataset for an example of how the resulting file should look.
}
}
\details{
The plots produced include density plots for each parameter within each time bin (residualized to have a mean of zero), scaled so that the top of the density is at a value of one. Superimposed onto these densitys are the densities of a normal distribution with the same mean and variance (and scaled by the same amount). Deviations between the normal density and the density of the parameters indiciate deviations from normality. The standard deviation of each parameter is also displayed for each time bin to facilitate assessing homogenity of variance.
}
\value{
A \code{ggplot} object, which can be modified using \pkg{ggplot2} functions.
}

\seealso{
\code{\link{import_log}} for producing a data set of FBD parameter posterior samples from parameter log files.

\code{\link{FBD_tests1}} for statistical tests of normality and homogeneity of variance.

\code{\link{FBD_tests2}} for tests of differences in parameter means.
}
\examples{
data("posterior1p")

FBD_normality_plot(posterior1p)
}
