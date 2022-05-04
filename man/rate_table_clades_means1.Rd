\name{rate_table_clades_means1}
\alias{rate_table_clades_means1}
\docType{data}
\title{
Mean clock rates by node and clade (single clock)
}
\description{
A data set containing the mean clock rates for a tree with 1 clock partition, such as the output of \code{\link{get_clockrate_table}} but with an additional "clade" column added, which is required for use in \code{\link{clockrate_summary}} and \code{\link{clockrate_dens_plot}}.
}
\usage{data("rate_table_clades_means1")}
\format{
  A data frame with 79 observations on the following 3 variables.
  \describe{
    \item{\code{clade}}{A character vector containing the clade names for each corresponding node}
    \item{\code{nodes}}{A numeric vector for the node numbers in the summary tree}
    \item{\code{rates}}{A numeric vector containing the mean posterior clock rate for each node}
  }
}
\details{
\code{rate_table_clades_means1} was created by running \code{get_clockrate_table(tree1p)} and then adding a "clade" column. It can be produced by using the following procedure:

1) Import tree file:

\preformatted{data("tree1p")}

2) Produce clock rate table with, for instance, mean rate values from each branch in the tree:

\preformatted{rate_table <- get_clockrate_table(tree1p, summary = "mean")}

\preformatted{write.csv(rate_table, file = "rate_table.csv", row.names = FALSE)}

3) Now, manually add clades using, e.g., Excel:

3.1) Manually edit rate_table.csv, adding a "clade" column. This introduces customized clade names to individual nodes in the tree.

3.2) Save the edited rate table with a different name to differentiate from the original output (e.g., rate_table_clades_means.csv).

4) Read the file back in:

\preformatted{rate_table_clades_means1 <- read.csv("rate_table_clades_means.csv")}

\preformatted{head(rate_table_clades_means1)}

}
\seealso{
\code{\link{tree1p}} for the tree from which the clock rates were extracted.

\code{\link{get_clockrate_table}} for extracting a clock rate table from a tree.

\code{\link{clockrate_summary}}, \code{\link{clockrate_dens_plot}}, and \code{\link{clockrate_reg_plot}} for examples of using a clockrate table.
}
\keyword{datasets}
