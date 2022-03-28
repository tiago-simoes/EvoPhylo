\name{FBD_dens_plot}
\alias{FBD_dens_plot}

\title{
Density Plot of One FBD Parameter
}
\description{
Produces a density or violin plot displaying the distribution of MCMC fossilized birth-death process (FBD) parameter samples by time bin.
}
\usage{
FBD_dens_plot(posterior, parameter, type = "density",
              stack = FALSE, color = "red")
}

\arguments{
  \item{posterior}{
a dataset of FBD parameter MCMC posterior samples imported using \code{\link{import_log}}. See the \code{\link{posterior1p}} dataset for an example of how the resulting file should look.
}
  \item{parameter}{
a string containing the name of an FBD parameter (\code{"net_speciation"}, \code{"relative_extinction"}, or \code{"relative_fossilization"}); abbreviations allowed.
}
  \item{type}{
the type of plot; either \code{"density"} for a density plot or \code{"violin"} for violin plots. Abbreviations allowed.
}
  \item{stack}{
when \code{type = "density"}, whether to produce stacked densities (\code{TRUE}) or overlapping densities (\code{FALSE}, the default). Ignored otherwise.
}
  \item{color}{
when \code{type = "violin"}, the color of the plotted densities.
}
}
\details{
Density plots are produced using \code{\link[ggplot2:stat_density]{ggplot2::stat_density}}, and violin plots are produced using \code{\link[ggplot2:geom_violin]{ggplot2::geom_violin}}. On violin plots, a horizontal line indicates the median (of the density), and the black dot indicates the mean.
}
\value{
A \code{ggplot} object, which can be modified using \pkg{ggplot2} functions.
}
\note{
When setting \code{type = "violin"}, a warning may appear saying something like "In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) : collapsing to unique 'x' values". This warning can be ignored.
}
\seealso{
\code{\link[ggplot2:stat_density]{ggplot2::stat_density}}, \code{\link[ggplot2:geom_violin]{ggplot2::geom_violin}} for the underlying functions to produce the plots.

\code{\link{import_log}} for producing a data set of FBD parameter posterior samples from parameter log files.

\code{\link{FBD_summary}}, \code{\link{FBD_normality_plot}}, \code{\link{FBD_tests1}}, and \code{\link{FBD_tests2}} for other functions used to summarize and display the distributions of the parameters.
}
\examples{
data("posterior1p")

FBD_dens_plot(posterior1p, parameter = "net_speciation",
              type = "density")

FBD_dens_plot(posterior1p, parameter = "net_speciation",
              type = "violin", color = "blue")
}
