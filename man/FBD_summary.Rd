\name{FBD_summary}
\alias{FBD_summary}

\title{
Summarize Parameter Values
}
\description{
Produces numerical summaries of each parameter (net speciation, relative extinction, and relative fossilization) by time bin and analysis type (tip, tip + node).
}
\usage{
FBD_summary(AllRunsrMelted_MC, file = NULL, digits = 3)
}

\arguments{
  \item{AllRunsrMelted_MC}{
(?) Dataset with analysis, Time_bin, and each parameter
}
  \item{file}{
an optional file path where the resulting table will be stored using \code{\link{write.csv}}.
}
  \item{digits}{
the number of digitis to round the summary results to. Default is 3. See \code{\link{round}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame with a row for each paramater, time bin, and analysis type, and columns for different summary statistics. These include the number of data points (\code{n}) and the mean, standard deviation (\code{sd}), minimum (\code{min}), first quartile (\code{Q1}), median, third quartile (\code{Q3}), and maximum (\code{max}). Rows with the analysis type "Combined" correspond to statistics collapsing across Tip and Tip + Node analyses for each parameter and time bin.
}

\author{

}
\note{

}

\seealso{

}
\examples{
}
