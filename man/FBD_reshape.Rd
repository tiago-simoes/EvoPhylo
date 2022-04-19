\name{FBD_reshape}
\alias{FBD_reshape}
\title{
Converts an FBD posterior parameter table from wide to long format
}
\description{
Converts FBD posterior parameter table, such as those imported using \code{\link{combine_log}}, from wide to long format.
}
\usage{
FBD_reshape(posterior)
}
\arguments{
  \item{posterior}{
Single posterior parameter sample dataset with skyline FBD parameters produced with \code{\link{combine_log}}. See the \code{\link{posterior3p}} or \code{\link{posterior1p}} dataset for an example of how the resulting file should look.

}
}
\details{
The posterior parameters log files produced by Bayesian evolutionary analyses using skyline birth-death tree models, including the skyline FBD model, result into 2 or more estimates for each FBD parameter (e.g., "net diversification", "relative extinction", and "relative fossilization"). This function will convert parameters log files including FBD parameters from wide to long format, with one row per generation per time bin and a new column "Time_bin" containing the time bins as a factor. The long format is necessary for downstream analyses, as similarly done by \code{\link{clock_reshape}} for clock rate tables.

The "posterior" data frame can be obtained by combining several output log files (.p files) from Mr. Bayes using \code{\link{combine_log}}.
}

\value{
A data frame containing a single "value" column and one column for each FBD parameter (e.g.,  "net_speciation", "relative extinction", and "relative fossilization"), indicating to which FBD parameter each value refers to.
}

\seealso{
\code{\link{combine_log}}, \code{\link{clock_reshape}}
}
\examples{

## The example dataset rate_table_clades already has clades and 3 clock rate columns:
data("posterior3p")

head(posterior3p)

## Reshape FBD table to long format
posterior3p_long <- FBD_reshape(posterior3p)

head(posterior3p_long)
}


