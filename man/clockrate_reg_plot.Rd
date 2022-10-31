\name{clockrate_reg_plot}
\alias{clockrate_reg_plot}
\title{
Plot regression lines between sets of rates
}
\description{
Displays a scatterplot and fits regression line of one set of clock rates against another, optionally displaying their Pearson correlation coefficient (r) and R-squared values (R^2).
}
\usage{
clockrate_reg_plot(rate_table, clock_x, clock_y,
                   method = "lm", show_lm = TRUE,
                   ...)
}
\arguments{
  \item{rate_table}{
A table of clock rates, such as from the output of \code{\link{get_clockrate_table_MrBayes}}.
}
  \item{clock_x, clock_y}{
The clock rates that should go on the x- and y-axes, respectively.
}
  \item{method}{
The method (function) used fit the regression of one clock on the other. Check the \code{method} argument in the to \code{\link{geom_smooth}} function of \pkg{ggplot2} for all options. Default is \code{"lm"} for a linear regression model. \code{"glm"} and \code{"loess"} are alternative options.
}
  \item{show_lm}{
Whether to display the Pearson correlation coefficient (r) and R-squared values (R^2) between two sets of clock rates.
}
  \item{\dots}{
Other arguments passed to \code{\link{geom_smooth}}.
}

}
\details{
\code{clockrate_reg_plot()} can only be used when multiple clocks are present in the clock rate table. Unlike \code{\link{clockrate_summary}} and \code{\link{clockrate_dens_plot}}, no "clade" column is required.
}
\value{
A \code{ggplot} object, which can be modified using \pkg{ggplot2} functions.
}


\seealso{
\code{vignette("rates-selection")} for the use of this function as part of an analysis pipeline.

\code{\link{geom_point}}, \code{\link{geom_smooth}}
}
\examples{
# See vignette("rates-selection") for how to use this
# function as part of an analysis pipeline

data("RateTable_Means_3p_Clades")

#Plot correlations between clocks 1 and 3
clockrate_reg_plot(RateTable_Means_3p_Clades,
                   clock_x = 1, clock_y = 3)

#Use arguments supplied to geom_smooth():
clockrate_reg_plot(RateTable_Means_3p_Clades,
                   clock_x = 1, clock_y = 3,
                   color = "red", se = FALSE)
}
