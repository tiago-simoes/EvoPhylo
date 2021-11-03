\name{FBD_normality_plot}
\alias{FBD_normality_plot}

\title{
Inspect parameter distributions visually
}
\description{
\code{FBD_normality_plot()} produces plots of the distributions of parameters (net_speciation, relative_exinction, and relative_fossilization) to facilitate the assessment of the assumptions of normality within time bins and homogeneity of variance across time bins.
}
\usage{
FBD_normality_plot(AllRunsrMelted_MC)
}

\arguments{
  \item{AllRunsrMelted_MC}{
an AllRunsMelted_MC (long) dataset with a Time_bin column and columns for each parameter (net_speciation, relative_exinction, and relative_fossilization).
}
}
\details{
The plots produced include density plots for each parameter within each time bin (residualized to have a mean of zero), scaled so that the top of the density is at a value of one. Superimposed onto these densitys are the densities of a normal distribution with the same mean and variance. Deviations between the normal density and the density of the parameters indiciate deviations from normality. The standard deviation of each parameter is also displayed for each time bin to facilitate assessing homogenity of variance.
}
\value{
A \code{ggplot} object, which can be modified using \pkg{ggplot2} functions.
}

\seealso{
\code{\link{FBD_tests1}} for statistical tests of normality and homogeneity of variance.
}
\examples{

}
