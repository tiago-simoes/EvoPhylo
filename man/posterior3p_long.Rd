\name{posterior3p_long}
\alias{posterior3p_long}
\docType{data}
\title{
Posterior parameter samples (3-partion analysis)
}
\description{
An example dataset similar to posterior3p, but reshaped from wide to long table format for some downstream analyses.

}
\usage{data("posterior3p_long")}
\format{
  Converting the posterior sample file from wide to long format will introduce a new column ("Time_bin") and change the other columns related to FBD variables as following:
  \describe{
    \item{\code{Time_bin}}{a factor with levels \code{1} \code{2} \code{3} \code{4}, etc. indicating the time bins from where FBD parameters have been sampled (only when using the skyline FBD model)}
    \item{\code{net_speciation}}{a numeric vector for net speciation estimates within each time bin}
    \item{\code{relative_extinction}}{a numeric vector for relative extinction estimates within each time bin}
    \item{\code{relative_fossilization}} {a numeric vector for relative fossilization estimates within each time bin}
  }
}
\seealso{
\code{\link{posterior3p}}
\code{\link{posterior1p}}
\code{\link{posterior1p_long}}
}
\keyword{datasets}



