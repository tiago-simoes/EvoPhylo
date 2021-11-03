\name{make_clusters}
\alias{make_clusters}

\title{
Produce Character Clusters
}
\description{
Determines cluster membership using the supplied Gower distance matrix and requested number of clusters, optionally plotting the result.
}
\usage{
make_clusters(Dmatrix, k, tsne = FALSE, plot = tsne,
              ...)
}
\arguments{
  \item{Dmatrix}{
a Gower distance matrix, the output of a call to \code{\link{get_gower_dist}}.
}
  \item{k}{
the desired number of clusters.
}
  \item{tsne}{
whether to perform Barnes-Hut t-distributed stochastic neighbor embedding to produce a two-dimensional representation of the distance matrix using \code{\link[Rtsne:Rtsne]{Rtsne::Rtsne}}. See Details. Default is \code{FALSE}.
}
  \item{plot}{
whether to plot cluster membership. If \code{TRUE}, \code{\link{plot.cluster_df}} is called on the resulting output. See \code{\link{plot.cluster_df}} for details. Default is \code{TRUE} when \code{tsne = TRUE} and \code{FALSE} otherwise.
}
  \item{\dots}{
further argument to \code{\link{plot.cluster_df}} when \code{plot = TRUE}.
}
}
\details{
\code{make_clusters} calls \code{\link[cluster:pam]{cluster::pam}} on the supplied Gower distance matrix with the specified number of clusters to determine cluster membership for each character. When \code{tsne = TRUE}, Barnes-Hut t-distributed stochastic neighbor embedding is used to compute a two-dimensional embedding of the distance matrix, which can help in visualizing the resulting clusters. \code{\link[Rtsne:Rtsne]{Rtsne::Rtsne}} is used to do this, with \code{theta = 0} and \code{dims = 2}. The resulting dimensions will be included in the output; see Value below.
}
\value{
A data frame, inheriting from class \code{"cluster_df"}, with a row for each character with its number (\code{character_number}) and cluster membership (\code{cluster}). When \code{tsne = TRUE}, two additional columns, labeled \code{tSNE_Dim1} and \code{tSNE_Dim2}, will be included, containing the values on the two dimensions computed using \code{Rtsne}.

The \code{pam} fit resulting from \code{cluster::pam} is returned in the \code{"pam.fit"} attribute of the outut object.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
When \code{plot = TRUE}, a warning may appear from \code{ggrepel}, saying something along the lines of "unlabeled data points (too many overlaps). Consider increasing max.overlaps". See \code{\link[ggrepel:geom_text_repel]{ggrepel::geom_text_repel}} for details; the \code{max.overlaps} argument can be supplied to \code{\dots}.
}

\seealso{
\code{\link{get_gower_dist}}, \code{\link{get_sil_widths}}, \code{\link{plot.cluster_df}}

\code{\link[cluster:pam]{cluster::pam}}, \code{\link[Rtsne:Rtsne]{Rtsne::Rtsne}}
}
\examples{

}
