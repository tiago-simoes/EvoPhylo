\name{FBD_reshape}
\alias{FBD_reshape}
\title{
Convert an FBD posterior parameter table from wide to long format
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
  \item{variables}{
  Names of FBD rate variables in the log. If NULL (default), will attempt to auto-detect the names and log type.
  }
  \item{log.type} {
  Name of the software which produced the log (currently supported: MrBayes or BEAST2). Has to be set if \code{variables} is not NULL.
  }
}
\details{
The posterior parameters log files produced by Bayesian evolutionary analyses using skyline birth-death tree models, including the skyline FBD model, result into two or more estimates for each FBD parameter, one for each time bin. This function will convert a table of parameters with skyline FBD parameters from wide to long format, with one row per generation per time bin and a new column "Time_bin" containing the respective time bins as a factor. The long format is necessary for downstream analyses using \code{\link{FBD_summary}}, \code{\link{FBD_dens_plot}}, \code{\link{FBD_normality_plot}}, \code{\link{FBD_tests1}}, or \code{\link{FBD_tests2}}, as similarly done by \code{\link{clock_reshape}} for clock rate tables.

The format of the log files can either be specified using the \code{variables} and \code{log.type} or auto-detected by the function.
The "posterior" data frame can be obtained by reading in a log file directly (e.g. using the \code{read.table} function) or by combining several output log files from Mr. Bayes using \code{\link{combine_log}}.
}

\value{
A data frame of posterior parameter estimates containing a single "Time_bin" column and one column for each FBD parameter value.
}

\seealso{
\code{vignette("fbd-params")} for the use of this function as part of an analysis pipeline.

\code{\link{combine_log}}, \code{\link{reshape}}
}
\examples{
# See vignette("fbd-params") for how to use this
# function as part of an analysis pipeline

data("posterior3p")

head(posterior3p)

## Reshape FBD table to long format
posterior3p_long <- FBD_reshape(posterior3p)

head(posterior3p_long)
}


