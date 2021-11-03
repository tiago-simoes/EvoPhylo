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
                   method = "lm", corr = TRUE)
}
\arguments{
  \item{rate_table}{
a table of clock rates, such as from the output of \code{\link{get_clockrate_table}}.
}
  \item{clock_x, clock_y}{
which clocks should go on the x- and y-axes, respectively.
}
  \item{method}{
the regression method used fit the regression of one clock on the other. Passed to \code{\link{geom_smooth}}. Default is \code{"lm"} for a linear regression.
}
  \item{corr}{
whether to display the correlation between the two sets of clock rates on the plot.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A \code{ggplot} object, which can be modified using \pkg{ggplot2} functions.

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{geom_point}}, \code{\link{geom_smooth}}
}
\examples{

}
