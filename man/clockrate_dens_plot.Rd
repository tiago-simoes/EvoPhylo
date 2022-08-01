\name{clockrate_dens_plot}
\alias{clockrate_dens_plot}
\title{
Plot clock rate distributions
}
\description{
Plots the distribution density of clock rates by clock and clade. The input must have a "clade" column.
}
\usage{
clockrate_dens_plot(rate_table, clock = NULL,
                    stack = FALSE, nrow = 1,
                    scales = "fixed")
}
\arguments{
  \item{rate_table}{
A data frame of clock rates, such as from the output of \code{\link{get_clockrate_table_MrBayes}} with an extra "clade" column.
}
  \item{clock}{
Which clock rates will be plotted. If unspecified, all clocks are plotted.
}
  \item{stack}{
Whether to display stacked density plots (\code{TRUE}) or overlapping density plots (\code{FALSE}).
}
  \item{nrow}{
When plotting rates for more than one clock, how many rows should be filled by the plots. This is passed to \code{\link{facet_wrap}}.
}
  \item{scales}{
When plotting rates for more than one clock, whether the axis scales should be "fixed" (default) across clocks or allowed to vary ("free", "free_x", or "free_y"). This is passed to \code{\link{facet_wrap}}.
}
}
\details{
The user must manually add clades to the rate table produced by \code{\link{get_clockrate_table_MrBayes}} before it can be used with this function. This can be doen manually with in R, such as by using a graphical user interface for editing data like the \pkg{DataEditR} package, or by writing the rate table to a spreadsheet and reading it back in after adding the clades. The example below uses a table that has had the clades added.
}
\value{
A \code{ggplot} object, which can be modified using \pkg{ggplot2} functions.
}

\seealso{
\code{vignette("rates-selection")} for the use of this function as part of an analysis pipeline.

\code{\link{get_clockrate_table_MrBayes}}, \code{\link{geom_density}}
}
\examples{
# See vignette("rates-selection") for how to use this
# function as part of an analysis pipeline

data("rate_table_clades_means3")

# Overlapping plots
clockrate_dens_plot(rate_table_clades_means3, stack = FALSE,
                    nrow = 1, scales = "fixed")

# Stacked density for all three clocks, changing the color
# palette to viridis using ggplot2 functions
clockrate_dens_plot(rate_table_clades_means3,
                    clock = 1:3, nrow = 1, stack = TRUE,
                    scales = "fixed") +
  ggplot2::scale_color_viridis_d() +
  ggplot2::scale_fill_viridis_d()
}
