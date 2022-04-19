\name{rate_table_clades_means}
\alias{rate_table_clades_means}
\docType{data}
\title{
Table with mean clock rates by node and clade
}
\description{
A dataset containing the mean clock rates for a tree with 3 clock partitions, such as the output of \code{\link{get_clockrate_table}} but with an additional "clade'" column added, which is required for use in \code{\link{clockrate_summary}} and \code{\link{clockrate_dens_plot}}.
}
\usage{data("rate_table_clades_means")}
\format{
  A data frame with 79 observations on the following 5 variables.
  \describe{
    \item{\code{clade}}{a character vector containing the clade correspondign to each node}
    \item{\code{nodes}}{a numeric vector for the node numbers in the summary tree}
    \item{\code{rates1}}{a numeric vector containing the mean posterior clock rate for each node for the first partition}
    \item{\code{rates2}}{a numeric vector containing the mean posterior clock rate for each node for the second partition}
    \item{\code{rates3}}{a numeric vector containing the mean posterior clock rate for each node for the third partition}
  }
}
\details{
\code{rate_table_clades_means} was formed by running \code{get_clockrate_table(tree3p)} and then adding a "clade" column. It can be produced by using the following procedure:

#Import tree file
data("tree3p")

#Produce clock rate table with, for instance, median rate values from each branch in the tree
rate_table <- get_clockrate_table(tree3p, summary = "median")

#Now, manually add clades using, e.g., Excel:
\dontrun{
write.csv(rate_table, file = "rate_table.csv",
          row.names = FALSE)

##Then,
# 1) Manually edit rate_table.csv, adding a "clade" column.
#    This introduces customized clade names to individual nodes in the tree.
# 2) Save the edited rate table with a different name to differentiate from the original output
#    (e.g., rate_table_clades_means.csv)
# 3) Read the file back in:

rate_table_clades_means <- read.csv("rate_table_clades_means.csv")
}

head(rate_table_clades_means)


}
\seealso{
\code{\link{tree3p}} for the tree from which the clock rates were extracted.

\code{\link{get_clockrate_table}} for extracting a clock rate table from a tree.

\code{\link{clockrate_summary}}, \code{\link{clockrate_dens_plot}}, and \code{\link{clockrate_reg_plot}} for examples of using a clockrate table.
}
\keyword{datasets}
