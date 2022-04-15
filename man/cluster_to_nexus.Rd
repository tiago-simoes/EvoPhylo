\name{cluster_to_nexus}
\alias{cluster_to_nexus}

\title{
Exports character clusters to a nexus file
}
\description{
Creates and exports a Nexus file with a list of characters and their respective partitions as inferred by the \code{make_clusters} function. The contents can be copied and pasted directly into a Mr. Bayes commands block for a partitioned clock Bayesian inference analysis.
}
\usage{
cluster_to_nexus(cluster_df, file = NULL)
}

\arguments{
  \item{cluster_df}{
A \code{cluster_df} object; the output of a call to \code{\link{make_clusters}}.
}
  \item{file}{
The path of the text file to be created containing the partitioning information in Nexus format. If \code{NULL} (the default), no file will be written and the output will be returned as a string. If \code{""}, the text will be printed to the console. Passed directly to the \code{file} argument of \code{\link{cat}}.
}
}
\value{
The text as a string, returned invisibly if \code{file} is not \code{NULL}. Use \code{\link{cat}} on the resulting output to format it correctly (i.e., to turn \code{"\\n"} into line breaks).
}

\seealso{
\code{\link{cat}}
}
\examples{
#Load example phylogenetic data matrix
data("characters")

#Create distance matrix
Dmatrix <- get_gower_dist(characters)

#Find optimal partitioning scheme using PAM under k=3 partitions
cluster_df <- make_clusters(Dmatrix, k = 3)

#Name of the desired .nex file, e.g., "clusters.nex"
file <- tempfile(fileext = ".nex")

#Write to nexus file and store text in cluster_file variable
cluster_file <- cluster_to_nexus(cluster_df, file = file)

cat(cluster_file)
}
