\name{FBD_dens_plot}
\alias{FBD_dens_plot}

\title{
Density Plot of One Parameter
}
\description{
Produces a density or violin plot of a single parameter by time bin.
}
\usage{
FBD_dens_plot(AllRunsrMelted_MC, parameter, type = "density",
              stack = FALSE, color = "red")
}

\arguments{
  \item{AllRunsrMelted_MC}{
(?) Dataset with analysis, Time_bin, and each parameter
}
  \item{parameter}{
a string containing the name of a parameter ("net_speciation", "relative_extinction", or "relative_fossilization"); abbreviations allowed.
}
  \item{type}{
the type of plot; either "density" for a density plot or "violin" for violins plots. Abbreviations allowed.
}
  \item{stack}{
when \code{type = "density"}, whether ot produce stacked densities (\code{TRUE}) or overlapping densities (\code{FALSE}, the defaults). Ignored otherwise.
}
  \item{color}{
when \code{type = "violin"}, the color of the plotted densities.
}
}
\details{
Density plots are produced using \code{\link[ggplot2:stat_density]{ggplot2::stat_density}}, and violin plots are produced using \code{\link[ggplot2:geom_violin]{ggplot2::geom_violin}}. On violin plots, a horizontal line indicates the median, and the black dot indicates the mean.
}
\value{
A \code{ggplot} object, which can be modified using \pkg{ggplot2} functions.
}

\author{

}
\note{

}

\seealso{
\code{\link[ggplot2:stat_density]{ggplot2::stat_density}}, \code{\link[ggplot2:geom_violin]{ggplot2::geom_violin}}
}
\examples{

}
