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
A phylogenetic data matrix in Nexus (.nex) format, or in any other data frame or matrix format with a column for each character and terminal taxa as rows, which will be read using \code{\link[ape:read.nexus.data]{ape::read.nexus.data}}. The data cannot include polymorphisms.
}
  \item{numeric}{
Whether to treat the values contained in the \code{x} as numeric or categorical. If \code{FALSE} (default), features will be considered categorical; if \code{TRUE}, they will be considered numeric. For further details, see "Vignette_ImportData".
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

#Reading example file as categorical data
Dmatrix <- get_gower_dist(characters)

#Reading tetrapod data matrix as categorical data (recommended)
Dmatrix <- get_gower_dist(DataMatrix.nex)

#Reading tetrapod data matrix as numeric data
Dmatrix <- get_gower_dist(characters, numeric = FALSE)
Dmatrix <- get_gower_dist(DataMatrix.nex, numeric = FALSE)
}
