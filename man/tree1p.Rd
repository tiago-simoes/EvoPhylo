\name{tree1p}
\alias{tree1p}
\docType{data}
\title{
Phylogenetic tree with a single clock partition
}
\description{
A clock Bayesian phylogenetic tree, imported as an S4 class object using \code{treeio::read.mrbayes()}.

}
\usage{data("tree1p")}
\format{
A \code{tidytree} object.
}
\details{
This example tree file was produced by analyzing the data set with a single morphological partition from Simões & Pierce (2021).
}
\seealso{
\code{\link{tree3p}} for another tree object with 3 clock partitions.

\code{\link{get_clockrate_table_MrBayes}} for extratcing the poserior clockrates from a tree object.
}
\keyword{datasets}

\references{
Simões, T. R. and S. E. Pierce (2021). Sustained High Rates of Morphological Evolution During the Rise of Tetrapods. \emph{Nature Ecology & Evolution} 5: 1403–1414.
}

