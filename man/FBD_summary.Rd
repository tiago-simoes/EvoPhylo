\name{FBD_summary}
\alias{FBD_summary}

\title{
Summarize FBD Parameter Values
}
\description{
Produces numerical summaries of each fossilized birth–death process (FBD) parameter (net speciation, relative extinction, and relative fossilization) by time bin.
}
\usage{
FBD_summary(posterior, file = NULL, digits = 3)
}

\arguments{
  \item{posterior}{
a dataset of FBD parameter MCMC posterior samples imported using \code{\link{import_log}}. See the \code{\link{posterior1p}} dataset for an example of how the resulting file should look.
}
  \item{file}{
an optional file path where the resulting table will be stored using \code{\link{write.csv}}.
}
  \item{digits}{
the number of digitis to round the summary results to. Default is 3. See \code{\link{round}}.
}
}

\value{
A data frame with a row for each paramater and time bin, and columns for different summary statistics. These include the number of data points (\code{n}) and the mean, standard deviation (\code{sd}), minimum (\code{min}), first quartile (\code{Q1}), median, third quartile (\code{Q3}), and maximum (\code{max}). When \code{file} is not \code{NULL}, a .csv file containing this data frame will be saved to the filepath specified in \code{file} and the output will be returned invisibly.
}

\seealso{
\code{\link{import_log}} for producing a data set of FBD parameter posterior samples from parameter log files.

\code{\link{FBD_dens_plot}}, \code{\link{FBD_normality_plot}}, \code{\link{FBD_tests1}}, and \code{\link{FBD_tests2}} for other functions used to summarize and display the distributions of the parameters.
}
\examples{
data("posterior1p")

FBD_summary(posterior1p)
}
