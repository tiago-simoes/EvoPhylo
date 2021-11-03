\name{cluster_to_nexus}
\alias{cluster_to_nexus}

\title{
Create a NEXUS text file of cluster membership
}
\description{
Creates a text file readable by NEXUS containing cluster membership.
}
\usage{
cluster_to_nexus(clusters, file = "")
}

\arguments{
  \item{clusters}{
a \code{cluster_df}; the output of a call to \code{\link{make_clusters}}.
}
  \item{file}{
the path of the text file to be created that is to contain the cluster information in a NEXUS-readblae format. If \code{""} (the default), the text will be printed to the console. Passed directly to the \code{file} argument of \code{\link{cat}}.
}
}
\details{

}
\value{
The text as a string, returned invisibly. Use \code{\link{cat}} on the resulting output to format it correctly (i.e., to turn \code{"\\n"} into line breaks).
}
\references{

}
\author{

}
\note{

}

\seealso{
\code{\link{cat}}
}
\examples{

}
