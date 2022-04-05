\name{FBD_reshape}
\alias{FBD_reshape}
\title{
Coverts FBD parameter estimates table from wide to long format
}
\description{
Coverts FBD parameter estimates table, such as those imported using \code{\link{import_log}}, from wide to long.
}
\usage{
FBD_reshape(samples)
}
\arguments{
  \item{samples}{
A data frame of parameter estimates from skyline birth-death tree models, such as the skyline fossilized birth-death (FBD) tree model. 
}
}
\details{
The log files of bayesian evolutionary analysis using skyline birth-death tree models, including the skyline FBD model, result into 2 or more estimates for each parameter (e.g., net diversification, relative extinction (turnover), and relative fossilization). This function will convert log files including FBD parameters from wide to long for downstream analyses, as similarly done by \code{\link{clock_reshape}} for clock rate tables.

The \code{samples} data frame can be obtained by combining several output log files (.p files) from Mr. Bayes using \code{\link{import_log}}, or using LogCombiner from the BEAST2 package. 
}
\value{
A data frame containing a single "value" column (for each FBD parameter values) and one column for each FBD parameter (e.g.,  "net_speciation", "relative extinction", and "relative fossilization"), indicating to which parameter partition each value refers to.
}

\seealso{
\code{\link{import_log}}, \code{\link{clock_reshape}} 
}
\examples{

## The example dataset rate_table_clades already has clades and 3 clock rate columns:
data("posterior1p")

head(posterior1p)

## Reshape FBD table to long format
posterior1p_long <- FBD_reshape(posterior1p)
}
