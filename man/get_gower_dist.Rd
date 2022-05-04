\name{get_gower_dist}
\alias{get_gower_dist}

\title{
Compute Gower distances between characters
}
\description{
Computes Gower distance between characters from a phylogenetic data matrix.
}
\usage{
get_gower_dist(x, numeric = FALSE)
}
\arguments{
  \item{x}{
A phylogenetic data matrix in Nexus (.nex) format, or in any other data frame or matrix format with a column for each character and terminal taxa as rows, which will be read using \code{\link[ape:read.nexus.data]{ape::read.nexus.data}}. The data cannot include polymorphisms.
}
  \item{numeric}{
Whether to treat the values contained in the \code{x} as numeric or categorical. If \code{FALSE} (default), features will be considered categorical; if \code{TRUE}, they will be considered numeric.
}
}
\author{
This function uses code adapted from \code{StatMatch::gower.dist()} written by Marcello D'Orazio.
}
\value{
The Gower distance matrix.
}

\examples{
# Load phylogenetic data matrix as categorical data and
# create distance matrix
\dontrun{
Dmatrix <- get_gower_dist("DataMatrix.nex")}

# Or load example phylogenetic data matrix
data("characters")

# Create distance matrix
Dmatrix <- get_gower_dist(characters)

#Reading data matrix as numeric data
Dmatrix <- get_gower_dist(characters, numeric = TRUE)
# OR
\dontrun{
Dmatrix <- get_gower_dist("DataMatrix.nex", numeric = TRUE)}
}
