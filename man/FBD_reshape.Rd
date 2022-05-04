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
Single posterior parameter sample dataset with skyline FBD parameters produced with \code{\link{combine_log}}.
}
}
\details{
The posterior parameters log files produced by Bayesian evolutionary analyses using skyline birth-death tree models, including the skyline FBD model, result into two or more estimates for each FBD parameter (e.g., "net diversification", "relative extinction", and "relative fossilization"), one for each time bin. This function will convert parameters (.p) log files with skyline FBD parameters from wide to long format, with one row per generation per time bin and a new column "Time_bin" containing the respective time bins as a factor. The long format is necessary for downstream analyses using \code{\link{FBD_summary}}, \code{\link{FBD_dens_plot}}, \code{\link{FBD_normality_plot}}, \code{\link{FBD_tests1}}, or \code{\link{FBD_tests2}}, as similarly done by \code{\link{clock_reshape}} for clock rate tables.

The "posterior" data frame can be obtained by combining several output log files (.p files) from Mr. Bayes using \code{\link{combine_log}}.
}

\value{
A data frame of posterior parameter estimates containing a single "Time_bin" column and one column for each FBD parameter value (e.g., "net_speciation", "relative extinction", and "relative fossilization").
}

\seealso{
\code{\link{combine_log}}, \code{\link{reshape}}
}
\examples{
data("posterior3p")

head(posterior3p)

## Reshape FBD table to long format
posterior3p_long <- FBD_reshape(posterior3p)

head(posterior3p_long)
}


