\name{posterior3p}
\alias{posterior3p}
\docType{data}
\title{
FBD parameter posterior samples (3-partion analysis)
}
\description{
A dataset of MCMC parameter samples resulting from processing MrBayes fossilized birth–death process (FBD) parameter output log files to yield a row per time bin per generation. This is an example of the kind of output that would be produced by \code{\link{import_log}}.
}
\usage{data("posterior3p")}
\format{
  A data frame with 4000 observations on the following 20 variables.
  \describe{
    \item{\code{Gen}}{a numeric vector}
    \item{\code{LnL}}{a numeric vector}
    \item{\code{LnPr}}{a numeric vector}
    \item{\code{TH}}{a numeric vector}
    \item{\code{TL}}{a numeric vector}
    \item{\code{prop_ancfossil}}{a numeric vector}
    \item{\code{sigma.1}}{a numeric vector}
    \item{\code{sigma.2}}{a numeric vector}
    \item{\code{sigma.3}}{a numeric vector}
    \item{\code{m.1}}{a numeric vector}
    \item{\code{m.2}}{a numeric vector}
    \item{\code{m.3}}{a numeric vector}
    \item{\code{tk02var.1}}{a numeric vector}
    \item{\code{tk02var.2}}{a numeric vector}
    \item{\code{tk02var.3}}{a numeric vector}
    \item{\code{clockrate}}{a numeric vector}
    \item{\code{Time_bin}}{a factor with levels \code{1} \code{2} \code{3} \code{4}}
    \item{\code{net_speciation}}{a numeric vector}
    \item{\code{relative_extinction}}{a numeric vector}
    \item{\code{relative_fossilization}}{a numeric vector}
  }
}
\seealso{
\code{\link{posterior1p}}
}
\keyword{datasets}
