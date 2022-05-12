\name{FBD_dens_plot}
\alias{FBD_dens_plot}

\title{
Density plots for each FBD parameter
}
\description{
Produces a density or violin plot displaying the distribution of FBD parameter samples by time bin.
}
\usage{
FBD_dens_plot(posterior, parameter, type = "density",
              stack = FALSE, color = "red")
}

\arguments{
  \item{posterior}{
A data frame of posterior parameter estimates containing a single "Time_bin" column and one column for each FBD parameter value (e.g.,  "net_speciation", "relative extinction", and "relative fossilization"). Such data frame can be imported using \code{\link{combine_log}} followed by \code{\link{FBD_reshape}}.
}
  \item{parameter}{
A string containing the name of an FBD parameter "net_speciation", "relative_extinction", or "relative_fossilization"; abbreviations allowed.
}
  \item{type}{
The type of plot; either \code{"density"} for a density plot or \code{"violin"} for violin plots. Abbreviations allowed.
}
  \item{stack}{
When \code{type = "density"}, whether to produce stacked densities (\code{TRUE}) or overlapping densities (\code{FALSE}, the default). Ignored otherwise.
}
  \item{color}{
When \code{type = "violin"}, the color of the plotted densities.
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
\code{vignette("fbd-params")} for the use of this function as part of an analysis pipeline.

\code{\link[ggplot2:stat_density]{ggplot2::stat_density}}, \code{\link[ggplot2:geom_violin]{ggplot2::geom_violin}} for the underlying functions to produce the plots.

\code{\link{combine_log}} for producing a single data frame of FBD parameter posterior samples from multiple log files.

\code{\link{FBD_reshape}} for converting a single data frame of FBD parameter estimates, such as those imported using \code{\link{combine_log}}, from wide to long format.

\code{\link{FBD_summary}}, \code{\link{FBD_normality_plot}}, \code{\link{FBD_tests1}}, and \code{\link{FBD_tests2}} for other functions used to summarize and display the distributions of the parameters.
}
\examples{
# See vignette("fbd-params") for how to use this
# function as part of an analysis pipeline

data("posterior3p")

posterior3p_long <- FBD_reshape(posterior3p)

FBD_dens_plot(posterior3p_long, parameter = "net_speciation",
              type = "density", stack = FALSE)

FBD_dens_plot(posterior3p_long, parameter = "net_speciation",
              type = "density", stack = TRUE)

FBD_dens_plot(posterior3p_long, parameter = "net_speciation",
              type = "violin", color = "red")
}
