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
a data frame of clock rates, such as from the output of \code{\link{get_clockrate_table}} with an extra \code{clade} column.
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
The user must manually add clades to the rate table produced by \code{\link{get_clockrate_table}} before it can be used with this function. This can be doen manually with in R, such as by using a graphical user interface for editing data like the \pkg{DataEditR} package, or by writing the rate table to a spreadsheet and reading it back in after adding the clades. The example below uses a table that has had the clades added.
}
\value{
A \code{ggplot} object, which can be modified using \pkg{ggplot2} functions.
}

\seealso{
\code{\link{get_clockrate_table}}, \code{\link{geom_density}}
}
\examples{
data("tree1p")

rate_table <- get_clockrate_table(tree1p, summary = "mean")

head(rate_table)

#Now, manually add clades using, e.g., Excel:
\dontrun{
write.csv(rate_table, file = "rate_table.csv",
          row.names = FALSE)
## Manually edit rate_table.csv, adding a "clade" column, then
## read the file back in

rate_table_clades <- read.csv("rate_table.csv")
}
## The example dataset rate_table_means already has clades:
data("rate_table_means")

head(rate_table_means)

clockrate_dens_plot(rate_table_means, clock = 1)

# Stacked density for all three clocks, changing the color
# palette to viridis using ggplot2 functions
clockrate_dens_plot(rate_table_means, clock = 1:3, nrow = 1,
                    stack = TRUE) +
  ggplot2::scale_color_viridis_d() +
  ggplot2::scale_fill_viridis_d()
}
