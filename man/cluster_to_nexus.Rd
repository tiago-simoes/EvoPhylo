\name{cluster_to_nexus}
\alias{cluster_to_nexus}

\title{
Create a NEXUS text file of cluster membership
}
\description{
Creates a text file readable by NEXUS containing cluster membership.
}
\usage{
cluster_to_nexus(cluster_df, file = NULL)
}

\arguments{
  \item{cluster_df}{
A \code{cluster_df} object; the output of a call to \code{\link{make_clusters}}.
}
  \item{file}{
The path of the text file to be created containing the cluster information in Nexus format (can be copied and pasted into a Mr. Bayes command block). If \code{NULL} (the default), no file will be written and the output will be returned as a string. If \code{""}, the text will be printed to the console. Passed directly to the \code{file} argument of \code{\link{cat}}.
}
}
\value{
The text as a string, returned invisibly if \code{file} is not \code{NULL}. Use \code{\link{cat}} on the resulting output to format it correctly (i.e., to turn \code{"\\n"} into line breaks).
}

\seealso{
\code{\link{cat}}
}
\examples{
data("characters")

#Reading example file as categorical data
Dmatrix <- get_gower_dist(characters)

#Create clusters with PAM under k=3 partitions
cluster_df <- make_clusters(Dmatrix, k = 3)

#Name of the desired .nex file, e.g., "clusters.nex"
file <- tempfile(fileext = ".nex")

#Write to nexus file and store text in cluster_file variable
cluster_file <- cluster_to_nexus(cluster_df, file = file)

cat(cluster_file)
}
