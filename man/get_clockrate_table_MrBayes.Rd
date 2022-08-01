\name{get_clockrate_table_MrBayes}
\alias{get_clockrate_table_MrBayes}

\title{
Extract evolutionary rates from a Bayesian clock tree produced by Mr. Bayes
}
\description{
Extract evolutionary rate summary statistics for each node from a Bayesian clock summary tree produced by Mr. Bayes and stores them in a data frame.
}
\usage{
get_clockrate_table_MrBayes(tree, summary = "median",
                    drop_dummy = NULL)
}
\arguments{
  \item{tree}{
An S4 class object of type \code{treedata}; a Bayesian clock tree imported using \code{\link[treeio:read.mrbayes]{treeio::read.mrbayes}} for Mr. Bayes summary trees.
}
  \item{summary}{
The name of the rate summary. Should be one of \code{"mean"} or \code{"median"}.
}
  \item{drop_dummy}{
if not \code{NULL}, will drop the dummy extant tip with the given label from the Mr. Bayes summary tree prior to extracting the clock rates (when present). Default is \code{NULL}.
}
}

\value{
A data frame with a column containing the node identifier (\code{node}) and one column for each relaxed clock partition in the tree object containing clock rates.
}

\seealso{
\code{vignette("rates-selection")} for the use of this function as part of an analysis pipeline.
\code{\link{get_clockrate_table_BEAST2}} for the equivalent function for BEAST2 output files.
\code{\link{clockrate_summary}} for summarizing and examining properties of the resulting rate table. Note that clade membership for each node must be customized (manually added) before these functions can be used, since this is tree and dataset dependent.
}
\examples{
# See vignette("rates-selection") for how to use this
# function as part of an analysis pipeline

## Import summary tree with three clock partitions produced by
## Mr. Bayes (.t or .tre files) from your local directory
\dontrun{
tree3p <- treeio::read.mrbayes("Tree3p.t")
}

#Or use the example Mr.Bayes multi-clock tree file (\code{tree3p})
data("tree3p")

# obtain the rate table from MrBayes tree
rate_table <- get_clockrate_table_MrBayes(tree3p)

head(rate_table)
}
