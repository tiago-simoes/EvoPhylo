\name{get_gower_dist}
\alias{get_gower_dist}

\title{
Compute Gower distances
}
\description{
Computes Gower distance between features from a phylogenetic data matrix
}
\usage{
get_gower_dist(file, numeric = FALSE)
}
\arguments{
  \item{file}{
a data frame or matrix with a row for each character containing feature qualities in the columns, or the path of a CSV file containing such a matrix, which will be read using \code{\link{read.csv}}.
}
  \item{numeric}{
whether to treat the values contained in the \code{file} as numeric or categorical. If \code{FALSE} 9default), features will be considered categorical; if \code{TRUE}, they will be considered numeric.
}
}
\details{
The Gower distance matrix is computed using \code{\link[StatMatch:gower.dist]{StatMatch::gower.dist}}.
}
\value{
The Gower distance matrix.
}
%\references{
%
%}
%\author{
%
%}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[StatMatch:gower.dist]{StatMatch::gower.dist}}
}
\examples{

}
