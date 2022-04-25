\name{posterior1p_long}
\alias{posterior1p_long}
\docType{data}
\title{
Posterior parameter samples (single clock) in long format
}
\description{
An example dataset similar to posterior1p, but reshaped from wide to long table format for some downstream analyses.

}
\usage{data("posterior1p_long")}
\format{
  Converting the posterior sample file from wide to long format will introduce a new column ("Time_bin") and change the other columns related to FBD variables as following:
  \describe{
    \item{\code{Time_bin}}{A factor with levels \code{1}, \code{2}, \code{3}, \code{4},...\emph{i} indicating the time bins from where FBD parameters have been sampled (only when using the skyline FBD model).}
    \item{\code{net_speciation}}{A numeric vector for net speciation estimates within each time bin.}
    \item{\code{relative_extinction}}{A numeric vector for relative extinction estimates within each time bin.}
    \item{\code{relative_fossilization}} {A numeric vector for relative fossilization estimates within each time bin.}
  }
}
\seealso{
\code{\link{posterior1p}}
\code{\link{posterior3p}}
\code{\link{posterior3p_long}}
}
\keyword{datasets}



