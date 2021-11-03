\name{clockrate_summary}
\alias{clockrate_summary}
\title{
Compute rate summary statistics across clade and clock
}
\description{
Computes summary statistics for each clade for each clock. The input must have a \code{clade} column.
}
\usage{
clockrate_summary(rate_table, file = NULL, digits = 3)
}
\arguments{
  \item{rate_table}{
a table of clock rates, such as from the output of \code{\link{get_clockrate_table}} with an extra \code{clade} column.
}
  \item{file}{
an optional file path where the resulting table will be stored using \code{\link{write.csv}}.
}
  \item{digits}{
the number of digits to round the summary results to. Default is 3. See \code{\link{round}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame containing a row for each clade and each clock with summary statistics (n, mean, standard deviation, minimum, 1st quartile, median, third quartile, maximum).
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{get_clockrate_table}}, \code{\link{summary}}
}
\examples{

}
