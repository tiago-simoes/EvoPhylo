\name{clockrate_summary}
\alias{clockrate_summary}
\title{
Compute rate summary statistics across clades and clocks
}
\description{
Computes summary statistics for each clade and/or each clock partition. The input must have a "clade" column.
}
\usage{
clockrate_summary(rate_table, file = NULL, digits = 3)
}
\arguments{
  \item{rate_table}{
A data frame of clock rates, such as from the output of \code{\link{get_clockrate_table_MrBayes}} with an extra \code{"clade"} column.
}
  \item{file}{
An optional file path where the resulting table will be stored using \code{\link{write.csv}}.
}
  \item{digits}{
The number of digits to round the summary results to. Default is 3. See \code{\link{round}}.
}
}
\details{
The user must manually add clades to the rate table produced by \code{\link{get_clockrate_table_MrBayes}} before it can be used with this function. This can be doen manually within R, such as by using a graphical user interface for editing data like the \pkg{DataEditR} package, or by writing the rate table to a spreadsheet and reading it back in after adding the clades. The example below uses a table that has had the clades added.
}
\value{
A data frame containing a row for each clade and each clock with summary statistics (n, mean, standard deviation, minimum, 1st quartile, median, third quartile, maximum).
}

\seealso{
\code{vignette("rates-selection")} for the use of this function as part of an analysis pipeline.

\code{\link{get_clockrate_table_MrBayes}}, \code{\link{summary}}
}
\examples{
# See vignette("rates-selection") for how to use this
# function as part of an analysis pipeline

data("RateTable_Means_3p_Clades")

clockrate_summary(RateTable_Means_3p_Clades)
}
