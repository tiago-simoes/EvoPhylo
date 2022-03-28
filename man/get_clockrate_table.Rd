\name{get_clockrate_table}
\alias{get_clockrate_table}

\title{
Extract clock rate summaries from a Tree
}
\description{
Extracts the clock rate summaries for each node from a tree and stores them in a data frame.
}
\usage{
get_clockrate_table(tree, summary = "mean",
                    drop_dummyextent = TRUE)
}
\arguments{
  \item{tree}{
a \code{treedata} object; the output of a call to \code{\link[treeio:read.mrbayes]{treeio::read.mrbayes}}.
}
  \item{summary}{
the name of the rate summary. Should be one of \code{"mean"} or \code{"median"}.
}
  \item{drop_dummyextent}{
\code{logical}; whether to drop the "DummyExtent" tip of the tree prior to extracting the clock rates. Default is \code{TRUE}.
}
}

\value{
A data frame containing a column containing the node identifer (\code{node}) and one column for each clock in the tree object containing the rates.
}

\seealso{
\code{\link{clockrate_summary}}, \code{\link{clockrate_summary}}, and \code{\link{clockrate_summary}} for summarizing and examining properties of the resulting rate table. Note that clade must be manually added before these functions can be used.
}
\examples{
data("tree1p")

rate_table <- get_clockrate_table(tree1p)

head(rate_table)
}
