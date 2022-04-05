\name{clockrate_summary}
\alias{clockrate_summary}
\title{
Compute rate summary statistics across clades and clocks
}
\description{
Computes summary statistics for each clade and/or for each clock partition. The input must have a \code{clade} column.
}
\usage{
clockrate_summary(rate_table, file = NULL, digits = 3)
}
\arguments{
  \item{rate_table}{
a data frame of clock rates, such as from the output of \code{\link{get_clockrate_table}} with an extra \code{clade} column.
}
  \item{file}{
an optional file path where the resulting table will be stored using \code{\link{write.csv}}.
}
  \item{digits}{
the number of digits to round the summary results to. Default is 3. See \code{\link{round}}.
}
}
\details{
The user must manually add clades to the rate table produced by \code{\link{get_clockrate_table}} before it can be used with this function. This can be doen manually within R, such as by using a graphical user interface for editing data like the \pkg{DataEditR} package, or by writing the rate table to a spreadsheet and reading it back in after adding the clades. The example below uses a table that has had the clades added.
}
\value{
A data frame containing a row for each clade and each clock with summary statistics (n, mean, standard deviation, minimum, 1st quartile, median, third quartile, maximum).
}

\seealso{
\code{\link{get_clockrate_table}}, \code{\link{summary}}
}
\examples{
data("tree1p")

rate_table <- get_clockrate_table(tree1p, summary = "median")

head(rate_table)

#Now, manually add clades using, e.g., Excel:
\dontrun{
write.csv(rate_table, file = "rate_table.csv",
          row.names = FALSE)

##Then,          
## 1) Manually edit rate_table.csv, adding a "clade" column. S
## 2) Save the edited rate table with a different name to differentiate from the original output (e.g., rate_table_clades.csv). Then,
## 3) Read the file back in:

rate_table_clades <- read.csv("rate_table_clades.csv")
}
## The example dataset rate_table_clades already has clades:
data("rate_table_clades")

head(rate_table_clades)

clockrate_summary(rate_table_clades)
}
