\name{clockrate_reg_plot}
\alias{clockrate_reg_plot}
\title{
Plot the relationship between two sets of clock rates
}
\description{
Displays a scatterplot and fit regression line of one set of clock rates against another, optionally computing their correlation.
}
\usage{
clockrate_reg_plot(rate_table, clock_x, clock_y,
                   method = "lm", corr = TRUE,
                   ...)
}
\arguments{
  \item{rate_table}{
a table of clock rates, such as from the output of \code{\link{get_clockrate_table}}.
}
  \item{clock_x, clock_y}{
which clocks should go on the x- and y-axes, respectively.
}
  \item{method}{
the regression method used fit the regression of one clock on the other. Passed to \code{\link{geom_smooth}}. Default is \code{"lm"} for a linear regression. \code{"loess"} is another option.
}
  \item{corr}{
whether to display the correlation between the two sets of clock rates on the plot.
}
  \item{\dots}{
other arguments passed to \code{\link{geom_smooth}}.
}

}
\details{
\code{clockrate_reg_plot()} can only be used when multiple clocks are present in the clock rate table. Unlike \code{\link{clockrate_summary}} and \code{\link{clockrate_dens_plot}}, no clade column is required.
}
\value{
A \code{ggplot} object, which can be modified using \pkg{ggplot2} functions.
}


\seealso{
\code{\link{geom_point}}, \code{\link{geom_smooth}}
}
\examples{
data("rate_table_means")

#Plot correlations between clocks 1 and 3
clockrate_reg_plot(rate_table_means, clock_x = 1, clock_y = 3)

#Use arguments supplied to geom_smooth():
clockrate_reg_plot(rate_table_means, clock_x = 1, clock_y = 3,
                   color = "red", se = FALSE)
}
