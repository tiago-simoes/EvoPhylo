\name{FBD_normality_plot}
\alias{FBD_normality_plot}

\title{
Inspect FBD parameter distributions visually
}
\description{
Produces plots of the distributions of fossilized birth–death process (FBD) parameters to facilitate the assessment of the assumptions of normality within time bins and homogeneity of variance across time bins.
}
\usage{
FBD_normality_plot(posterior)
}

\arguments{
  \item{posterior}{
A data frame of posterior parameter estimates containing a single "Time_bin" column and one column for each FBD parameter value (e.g.,  "net_speciation", "relative extinction", and "relative fossilization"). Such data frame can be imported using \code{\link{combine_log}} followed by \code{\link{FBD_reshape}}.
}
}
\details{
The plots produced include density plots for each parameter within each time bin (residualized to have a mean of zero), scaled so that the top of the density is at a value of one (in \emph{black}). Superimposed onto these densitys are the densities of a normal distribution with the same mean and variance (and scaled by the same amount) (in \emph{red}). Deviations between the normal density in \emph{red} and the density of the parameters in \emph{black} indiciate deviations from normality. The standard deviation of each parameter is also displayed for each time bin to facilitate assessing homogenity of variance.
}
\value{
A \code{ggplot} object, which can be modified using \pkg{ggplot2} functions.
}

\seealso{
\code{\link{combine_log}} for producing a single data set of parameter posterior samples from individual parameter log files.

\code{\link{FBD_reshape}} for converting posterior parameter table from wide to long format.

\code{\link{FBD_tests1}} for statistical tests of normality and homogeneity of variance.

\code{\link{FBD_tests2}} for tests of differences in parameter means.
}
\examples{
data("posterior3p")

posterior3p_long <- FBD_reshape(posterior3p)

FBD_normality_plot(posterior3p_long)
}
