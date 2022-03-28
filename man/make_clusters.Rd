\name{make_clusters}
\alias{make_clusters}
\alias{plot.cluster_df}

\title{
Produce and Plot Character Clusters
}
\description{
Determines cluster membership using the supplied Gower distance matrix and requested number of clusters. The results can be plotted using an estimated multi-dimensional representation.
}
\usage{
make_clusters(dist_mat, k, tsne = FALSE,
              tsne_dim = 2, tsne_theta = 0,
              ...)

\method{plot}{cluster_df}(x, seed = NA, nrow = 1,
              ...)
}
\arguments{
  \item{dist_mat}{
a Gower distance matrix, the output of a call to \code{\link{get_gower_dist}}.
}
  \item{k}{
the desired number of clusters.
}
  \item{tsne}{
whether to perform Barnes-Hut t-distributed stochastic neighbor embedding (tSNE) to produce a multi-dimensional representation of the distance matrix using \code{\link[Rtsne:Rtsne]{Rtsne::Rtsne}}. The number of dimensions is controlled b ythe \code{tsne_dim} argument. See Details. Default is \code{FALSE}.
}
  \item{tsne_dim}{
when \code{tsne = TRUE}, the number of dimensions to use to represent the distance matrix. This is passed to the \code{dims} argument of \code{\link[Rtsne:Rtsne]{Rtsne::Rtsne}}. Default is 2.
}
  \item{tsne_theta}{
when \code{tsne = TRUE}, a parameter controlling the speed/accuracy trade-off (increase for faster but less accurate results). This is passed to the \code{theta} argument of \code{\link[Rtsne:Rtsne]{Rtsne::Rtsne}}. Default is 0 for exact tSNE.
}
  \item{\dots}{
for \code{make_clusters()}, other arguments passed to \code{\link[Rtsne:Rtsne]{Rtsne::Rtsne}} when \code{tsne = TRUE}.

For \code{plot.cluster_df()}, other arguments passed to \code{\link[ggrepel:geom_label_repel]{ggrepel::geom_text_repel}} to control display of the observation labels.
}
  \item{x}{
a \code{cluster_df} object; the output of a call to \code{make_clusters()}.
}
  \item{seed}{
the seed used to control the placement of the labels and the jittering of the points. Jittering only occurs when \code{tsne = FALSE} in the call to \code{make_clusters()}. Using a non-\code{NA} seed ensure replicability across uses.
}
  \item{nrow}{
when \code{tsne = TRUE} in the call to \code{make_clusters()} and \code{tsne_dim} is greater than 2, the number of rows used to display the resulting 2-dimensional plots. Default is 1 for side-by-side plots.
}
}
\details{
\code{make_clusters} calls \code{\link[cluster:pam]{cluster::pam}} on the supplied Gower distance matrix with the specified number of clusters to determine cluster membership for each character. When \code{tsne = TRUE}, Barnes-Hut t-distributed stochastic neighbor embedding is used to compute a multi-dimensional embedding of the distance matrix, which can help in visualizing the resulting clusters. \code{\link[Rtsne:Rtsne]{Rtsne::Rtsne}} is used to do this. The resulting dimensions will be included in the output; see Value below.

\code{plot.cluster_df()} plots cluster membership ina scatterplot with points colored based on cluster membership. When \code{tsne = TRUE} in the call to \code{make_clusters()}, the x- and y-axes will correspond to requested tSNE dimensions. With more than 2 dimensions, several plots will be produced, one for each pair of tSNE dimensions. These are displayed together using \code{\link[patchwork:plot_layout]{patchwork::plot_layout}}. When \code{tsne = FALSE}, the points will be arrange horizontally by cluster membership and randomly placed vertically.
}
\value{
A data frame, inheriting from class \code{"cluster_df"}, with a row for each character with its number (\code{character_number}) and cluster membership (\code{cluster}). When \code{tsne = TRUE}, additional columns  will be included, one for each requested tSNE dimension, labeled \code{tSNE_Dim1}, \code{tSNE_Dim2}, etc., containing the values on the dimensions computed using \code{Rtsne()}.

The \code{pam} fit resulting from \code{cluster::pam} is returned in the \code{"pam.fit"} attribute of the outut object.
}

\note{
When using \code{plot.cluster_df}, warnings may appear from \code{ggrepel}, saying something along the lines of "unlabeled data points (too many overlaps). Consider increasing max.overlaps". See \code{\link[ggrepel:geom_label_repel]{ggrepel::geom_text_repel}} for details; the \code{max.overlaps} argument can be supplied to \code{\dots}. This warning can generally be ignored, though.
}

\seealso{
\code{\link{get_gower_dist}}, \code{\link{get_sil_widths}}

\code{\link[cluster:pam]{cluster::pam}}, \code{\link[Rtsne:Rtsne]{Rtsne::Rtsne}}
}
\examples{
data("characters")

dist_matrix <- get_gower_dist(characters)

sil_widths <- get_sil_widths(dist_matrix, max.k = 7)

sil_widths
# 3 clusters yields the highest silhouette width

#Create clusters
cluster_df <- make_clusters(dist_matrix, k = 3)

#Plot clusters
plot(cluster_df, seed = 12345)

#Create clusters and estimate 3 tNSE dimensions
cluster_df_tsne <- make_clusters(dist_matrix, k = 3,
                                 tsne = TRUE, tsne_dim = 3)

#Plot clusters
plot(cluster_df_tsne, seed = 12345)
}
