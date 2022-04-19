\name{get_clockrate_table}
\alias{get_clockrate_table}

\title{
Extract evolutionary rates from a Bayesian clock tree
}
\description{
Extract evolutionary rate summary statistics for each node from a Bayesian clock (time) summary tree produced by Mr. Bayes and stores them in a data frame.
}
\usage{
get_clockrate_table(tree, summary = "median",
                    drop_dummyextant = TRUE)
}
\arguments{
  \item{tree}{
An S4 class object of type \code{treedata}; a Bayesian clock tree imported using \code{\link[treeio:read.mrbayes]{treeio::read.mrbayes}} for Mr. Bayes summary trees.
}
  \item{summary}{
The name of the rate summary. Should be one of \code{"mean"} or \code{"median"}.
}
  \item{drop_dummyextant}{
\code{logical}; whether to drop a "dummy" extant tip labeled as "Dummyextant" from the Mr. Bayes summary tree prior to extracting the clock rates (when present). Default is \code{TRUE}.
}
}

\value{
A data frame with a column containing the node identifer (\code{node}) and one column for each relaxed clock partition in the tree object containing clock rates.
}

\seealso{
\code{\link{clockrate_summary}} for summarizing and examining properties of the resulting rate table. Note that clade membership each node must be customized (manually added) before these functions can be used, since this is tree and dataset dependent.
}
\examples{
#Load example tree file
data("tree3p")

#OR

#Import example tree (.t or .tre) file from example directory
tree3p<-treeio::read.mrbayes("_dev/Examples/MultiClockTree/Tree_3p.t")

rate_table <- get_clockrate_table(tree3p)

head(rate_table)
}
