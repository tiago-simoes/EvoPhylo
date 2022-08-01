\name{clock_reshape}
\alias{clock_reshape}
\title{
Convert clock rate tables from wide to long format
}
\description{
Converts clock rate tables, such as those produced by \code{\link{clockrate_summary}} and imported back after including clade names, from wide to long format.
}
\usage{
clock_reshape(rate_table)
}
\arguments{
  \item{rate_table}{
A data frame of clock rates, such as from the output of \code{\link{get_clockrate_table_MrBayes}} with an extra "clade" column.
}
}
\details{
This function will convert clock rate tables from wide to long format, with a new column "clock" containing the clock partition from where each rate estimate was obtained as a factor. The long format is necessary for downstream analyses of selection strength (mode), as similarly done by \code{\link{FBD_reshape}} for posterior parameter log files.
}
\value{
A data frame containing a single "value" column (for all rate values) and one column for the "clock" variable (indicating to which clock partition each rate values refers to)
}

\seealso{
\code{vignette("rates-selection")} for the use of this function as part of an analysis pipeline.

\code{\link{get_clockrate_table_MrBayes}}, \code{\link{summary}}, \code{\link{clockrate_summary}}, \code{\link{FBD_reshape}}
}
\examples{
# See vignette("rates-selection") for how to use this
# function as part of an analysis pipeline

## The example dataset rate_table_clades_means3 already
## has clades and 3 clock rate columns:
data("rate_table_clades_means3")

head(rate_table_clades_means3)

## Reshape a clock rate table with clade names to long format
rates_by_clade <- clock_reshape(rate_table_clades_means3)
}


