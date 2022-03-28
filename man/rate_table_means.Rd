\name{rate_table_means}
\alias{rate_table_means}
\docType{data}
\title{
Clockrates by clade
}
\description{
A dataset containing the mean clockrates for a tree with three partitions, such as the output of \code{\link{get_clockrate_table}} but with an additional \code{clade} column added, which is required for use in \code{\link{clockrate_summary}} and \code{\link{clockrate_dens_plot}}.
}
\usage{data("rate_table_means")}
\format{
  A data frame with 79 observations on the following 5 variables.
  \describe{
    \item{\code{clade}}{a character vector containing the clade correspondign to each node}
    \item{\code{nodes}}{a numeric vector}
    \item{\code{rates1}}{a numeric vector containing the mean posterior clock rate for each node for the first partition}
    \item{\code{rates2}}{a numeric vector containing the mean posterior clock rate for each node for the second partition}
    \item{\code{rates3}}{a numeric vector containing the mean posterior clock rate for each node for the third partition}
  }
}
\details{
\code{rate_table_means} was formed by running \code{get_clockrate_table(tree3p)} and then adding a \code{close} column.
}
\seealso{
\code{\link{tree3p}} for the tree from which the clockrates were extracted.

\code{\link{get_clockrate_table}} for extracting a clockrate table from a tree.

\code{\link{clockrate_summary}}, \code{\link{clockrate_dens_plot}}, and \code{\link{clockrate_reg_plot}} for examples of using a clockrate table.
}
\keyword{datasets}
