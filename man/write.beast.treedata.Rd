\name{write.beast.treedata }
\alias{write.beast.treedata }

\title{
Export multiple treedata objects (S4 class tree files) to BEAST NEXUS file
}
\description{
This function was adopted and modified from treeio::write.beast to export a list of trees instead of a single tree.
}
\usage{
write.beast.treedata = function(treedata, file = "", translate = TRUE, tree.name = "STATE")

}
\arguments{
  \item{treedata}{
An S4 class object of type \code{treedata} containing multiple trees; a Bayesian clock tree imported using \code{\link[treeio:read.beast]{treeio::read.beast}} or \code{\link[treeio:read.mrbayes]{treeio::read.mrbayes}}.
}
  \item{file}{
Output file. If file = "", print the output content on screen
}
  \item{translate}{
Whether translate taxa labels
}
  \item{tree.name}{
Name of the tree
}
}

\value{
Writes object type \code{treedata} containing multiple trees to a file or file content on screen
}

\examples{
#Load file with multiple trees
trees_file = system.file("_dev/extdata/BEAST2/OffsetTrees", "ex_offset.trees", package = "EvoPhylo")

#Write multiple trees to a file
write.beast.treedata(posterior_trees_offset, file = "posterior_trees_offset.t")
}

