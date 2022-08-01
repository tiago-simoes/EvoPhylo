\name{write.beast.treedata }
\alias{write.beast.treedata }

\title{
Export multiple treedata objects (S4 class tree files) to BEAST NEXUS file
}
\description{
This function was adopted and modified from treeio::write.beast to export a list of trees instead of a single tree.
}
\usage{
write.beast.treedata(treedata, file = "", translate = TRUE, tree.name = "STATE")

}
\arguments{
  \item{treedata}{
An S4 class object of type \code{treedata} containing multiple trees; e.g. a Bayesian clock tree distribution imported using \code{\link[treeio:read.beast]{treeio::read.beast}} or \code{\link[treeio:read.mrbayes]{treeio::read.mrbayes}}.
}
  \item{file}{
Output file. If \code{file = ""}, prints the output content on screen.
}
  \item{translate}{
Whether to translate taxa labels.
}
  \item{tree.name}{
Name of the trees, default \code{"STATE"}.
}
}

\value{
Writes object type \code{treedata} containing multiple trees to a file or file content on screen
}

\examples{
#Load file with multiple trees
trees_file = system.file("extdata", "ex_offset.trees", package = "EvoPhylo")
posterior_trees_offset = treeio::read.beast(trees_file)

#Write multiple trees to screen
write.beast.treedata(posterior_trees_offset)
}

