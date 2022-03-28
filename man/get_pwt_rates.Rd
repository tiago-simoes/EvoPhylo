\name{get_pwt_rates}
\alias{get_pwt_rates}

\title{
TITLE
%Needs a title
}
\description{
Produces a data frame containing the results of 1-sample t-tests for the mean of posterior clockrates against each node's absolute clock rate.
}
\usage{
get_pwt_rates(rate_table, posterior)
}

\arguments{
  \item{rate_table}{
a data frame of clock rates, such as from the output of \code{\link{get_clockrate_table}} with an extra \code{clade} column.
}
  \item{posterior}{
a dataset of FBD parameter MCMC posterior samples imported using \code{\link{import_log}}. See the \code{\link{posterior1p}} dataset for an example of how the resulting file should look. The \code{clockrate} column from this dataset will be used in the t-tests.
}
}
\details{
\code{get_pwt_rates()} computes the absolute clockrate for each node and each clock as the corresponding rate multiplied by the mean of the posterior clockrates. Then, for each node and clock, a one-sample t-test is performed with the null hypothesis that the mean of the posterior clockrates is equal to that node and clock's absolute clock rate.
}
\value{
A long data frame with one row per node per clock and the following columns:
\item{clade}{the name of the clade, taken from the \code{clade} column of \code{rate_table}}
\item{nodes}{the node number, taken from the \code{node} column of \code{rate_table}}
\item{clock}{the clock number}
\item{relative rate}{the relative clockrate, taken from the \code{rates} columns of \code{rate_table}}
\item{absolute rate (mean)}{the absolute clockrate; the relative clock rate multiplied by the mean of the posterior clockrates}
\item{null}{the absolute clockrate used as the null value in the t-test}
\item{p.value}{the p-value of the test comparing the mean ofthe posterior clockrates to each absolute clockrate}
}

\seealso{

}
\examples{
data("rate_table_means")
data("posterior3p")

get_pwt_rates(rate_table_means, posterior3p)
}
