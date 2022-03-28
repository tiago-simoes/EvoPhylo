\name{get_gower_dist}
\alias{get_gower_dist}

\title{
Compute Gower distances
}
\description{
Computes Gower distance between features from a phylogenetic data matrix
}
\usage{
get_gower_dist(x, numeric = FALSE)
}
\arguments{
  \item{x}{
a data frame or matrix with a column for each character containing feature qualities in the rows, or the path of a .nex file containing such a phylogenetic data matrix, which will be read using \code{\link[ape:read.nexus.data]{ape::read.nexus.data}}. The data cannot include polymorphisms.
}
  \item{numeric}{
whether to treat the values contained in the \code{x} as numeric or categorical. If \code{FALSE} (default), features will be considered categorical; if \code{TRUE}, they will be considered numeric.
}
}
\author{
This function uses code adapted from \code{StatMatch::gower.dist()} written by Marcello D'Orazio.
}
\value{
The Gower distance matrix.
}

\seealso{
\code{\link[StatMatch:gower.dist]{StatMatch::gower.dist}}
}
\examples{
data("characters")

dist_matrix <- get_gower_dist(characters)
}
