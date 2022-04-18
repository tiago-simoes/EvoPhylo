\name{FBD_reshape}
\alias{FBD_reshape}
\title{
Coverts FBD parameter estimates table from wide to long format
}
\description{
Coverts FBD parameter estimates table, such as those imported using \code{\link{import_log}}, from wide to long format.
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
The log files produced by bayesian evolutionary analyses using skyline birth-death tree models, including the skyline FBD model, result into 2 or more estimates for each FBD parameter (e.g., net diversification, relative extinction (turnover), and relative fossilization). This function will convert log files including FBD parameters from wide to long format, with one row per generation per time bin and a new column "Time_bin" containing the time bins as a factor. The long format is necessary for downstream analyses, as similarly done by \code{\link{clock_reshape}} for clock rate tables.

The "samples" data frame can be obtained by combining several output log files (.p files) from Mr. Bayes using \code{\link{combine_log}}.
}

\value{
A data frame containing a single "value" column and one column for each FBD parameter (e.g.,  "net_speciation", "relative extinction", and "relative fossilization"), indicating to which FBD parameter each value refers to.
}

\seealso{
\code{\link{combine_log}}, \code{\link{clock_reshape}}
}
\examples{

## The example dataset rate_table_clades already has clades and 3 clock rate columns:
data("posterior1p")

head(posterior1p)

## Reshape FBD table to long format
posterior1p_long <- FBD_reshape(posterior1p)
}


