\name{clockrate_dens_plot}
\alias{clockrate_dens_plot}
\title{
Plot the density of clock rates
}
\description{
Plots the density of clock rates by clock and clade.
}
\usage{
clockrate_dens_plot(rate_table, clock = NULL,
                    stack = FALSE, nrow = 1,
                    scales = "fixed")
}
\arguments{
  \item{rate_table}{
a table of clock rates, such as from the output of \code{\link{get_clockrate_table}} with an extra \code{clade} column.
}
  \item{clock}{
for which clock(s) rates will be plotted. If unspecified, all clocks are plotted.
}
  \item{stack}{
whether to display stacked density plots (\code{TRUE}) or overlapping density plots (\code{FALSE}).
}
  \item{nrow}{
when plotting rates for more than one clock, how many rows should be filled by the plots. This is passed to \code{\link{facet_wrap}}.
}
  \item{scales}{
when plotting rates for more than one clock, whether the axis scales should be fixed across clocks or allowed to vary. This is passed to \code{\link{facet_wrap}}.
}
}
\details{

}
\value{
A \code{ggplot} object, which can be modified using \pkg{ggplot2} functions.
}

\note{

}

\seealso{
\code{\link{geom_density}}
}
\examples{

}
