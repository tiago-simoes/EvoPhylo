\name{get_pwt_rates}
\alias{get_pwt_rates}

\title{
Conduct pairwise t-tests between node rates and clock base rate
}
\description{
Produces a data frame containing the results of 1-sample t-tests for the mean of posterior clock rates against each node's absolute clock rate.
}
\usage{
get_pwt_rates(rate_table, posterior)
}

\arguments{
  \item{rate_table}{
A data frame containing a single "value" column (for all rate values) and one column for the "clock" variable (indicating to which clock partition each rate values refers to), such as from the output of \code{\link{get_clockrate_table}} with an extra \code{clade} column added, and followed by \code{\link{clock_reshape}}.
}
  \item{posterior}{
A data frame of posterior parameter estimates including a "clockrate" column indicating the base of the clock rate estimate for each generation that will be used for pairwise t-tests. Such data frame can be imported using \code{\link{combine_log}} (no need to reshape from wide to long). See the \code{\link{posterior1p}} or \code{\link{posterior3p}} datasets for an examples of how the input file should look.
}
}
\details{
\code{get_pwt_rates()} first transforms relative clock rates to absolute rate values for each node and each clock, by multiplying these by the mean posterior clock rate base value. Then, for each node and clock, a one-sample t-test is performed with the null hypothesis that the mean of the posterior clockrates is equal to that node and clock's absolute clock rate.
}
\value{
A long data frame with one row per node per clock and the following columns:
\item{clade}{The name of the clade, taken from the "clade" column of \code{rate_table}}
\item{nodes}{The node number, taken from the "node" column of \code{rate_table}}
\item{clock}{The clock partition number}
\item{relative rate}{The relative mean clock rate per node, taken from the "rates" columns of \code{rate_table}}
\item{absolute rate (mean)}{The absolute mean clock rate per node; the relative clock rate multiplied by the mean of the posterior clock rates}
\item{null}{The absolute clock rate used as the null value in the t-test}
\item{p.value}{The p-value of the test comparing the mean ofthe posterior clockrates to each absolute clockrate}
}

\seealso{
\code{vignette("rates-selection")} for the use of this function as part of an analysis pipeline.

\code{\link{combine_log}}
}
\examples{
# See vignette("rates-selection") for how to use this
# function as part of an analysis pipeline

# Load example rate table and posterior data sets
data("rate_table_clades_means3")
data("posterior3p")

get_pwt_rates(rate_table_clades_means3, posterior3p)
}
