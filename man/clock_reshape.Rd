\name{clock_reshape}
\alias{clock_reshape}
\title{
Coverts clock rate tables from wide to long format
}
\description{
Coverts clock rate tables, such as those produced by \code{\link{clockrate_summary}} and imported back after including clade names, from wide to long format. The input must have a "clade" column for downstream analyses.
}
\usage{
clock_reshape(rate_table)
}
\arguments{
  \item{rate_table}{
a data frame of clock rates, such as from the output of \code{\link{get_clockrate_table}} with an extra \code{clade} column.
}
}
\details{
The user must manually add clades to the rate table produced by \code{\link{get_clockrate_table}} before it can be used with this function. This can be doen manually within R, such as by using a graphical user interface for editing data like the \pkg{DataEditR} package, or by writing the rate table to a spreadsheet and reading it back in after adding the clades. The example below uses a table that has had the clades added.
}
\value{
A data frame containing a single "value" column (for all rate values) and one column for the "clock" variable (indicating to which clock partition each rate values refers to)
}

\seealso{
\code{\link{get_clockrate_table}}, \code{\link{summary}}, \code{\link{clockrate_summary}}
}
\examples{

## The example dataset rate_table_clades already has clades and 3 clock rate columns:
data("rate_table_clades")

head(rate_table_clades)

## Reshape a clock rate table with clade names to long format
RatesByClade <- clock_reshape(RateTable_Medians)
}


