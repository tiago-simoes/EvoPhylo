\name{FBD_summary}
\alias{FBD_summary}

\title{
Summarize FBD posterior parameter estimate values
}
\description{
Produces numerical summaries of each fossilized birth–death process (FBD) parameter ("net diversification", "relative extinction", and "relative fossilization") by time bin.
}
\usage{
FBD_summary(posterior, file = NULL, digits = 3)
}

\arguments{
  \item{posterior}{
A data frame containing a single "value" column and one column for each FBD parameter (e.g.,  "net_speciation", "relative extinction", and "relative fossilization"), indicating to which FBD parameter each value refers to. Such data frame can be imported using \code{\link{combine_log}} followed by \code{\link{FBD_reshape}} . See the \code{\link{posterior1p_long}} or \code{\link{posterior3p_long}} datasets for an examples of how the input file should look.
}
  \item{file}{
An optional file path where the resulting table will be stored using \code{\link{write.csv}}.
}
  \item{digits}{
The number of digitis to round the summary results to. Default is 3. See \code{\link{round}}.
}
}

\value{
A data frame with a row for each paramater and time bin, and columns for different summary statistics. These include the number of data points (\code{n}) and the mean, standard deviation (\code{sd}), minimum (\code{min}), first quartile (\code{Q1}), median, third quartile (\code{Q3}), and maximum (\code{max}). When \code{file} is not \code{NULL}, a .csv file containing this data frame will be saved to the filepath specified in \code{file} and the output will be returned invisibly.
}

\seealso{
\code{\link{combine_log}} for producing a single data set of parameter posterior samples from individual parameter log files.

\code{\link{FBD_reshape}} for converting posterior parameter table from wide to long format.

\code{\link{FBD_dens_plot}}, \code{\link{FBD_normality_plot}}, \code{\link{FBD_tests1}}, and \code{\link{FBD_tests2}} for other functions used to summarize and display the distributions of the parameters.
}
\examples{
data("posterior3p_long")

FBD_summary(posterior3p_long)
}
