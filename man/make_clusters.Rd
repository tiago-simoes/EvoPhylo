\name{make_clusters}
\alias{make_clusters}
\alias{plot.cluster_df}

\title{
Estimates and plots character partitions
}
\description{
Determines cluster(partition) membership for phylogenetic morphological characters from the supplied Gower distance matrix and requested number of clusters using partitioning around medoids (PAM, or K-medoids). For further and independently testing the quality of the chosen partitioning scheme, users may also poduce graphic clustering (tSNEs), coloring data points according to PAM clusters, to verify PAM clustering results.
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
A Gower distance matrix, the output of a call to \code{\link{get_gower_dist}}.
}
  \item{k}{
The desired number of clusters (or character partitions), the output from \code{\link{get_sil_widths}}.
}
  \item{tsne}{
Whether to perform Barnes-Hut t-distributed stochastic neighbor embedding (tSNE) to produce a multi-dimensional representation of the distance matrix using \code{\link[Rtsne:Rtsne]{Rtsne::Rtsne}}. The number of dimensions is controlled by the \code{tsne_dim} argument. See Details. Default is \code{FALSE}.
}
  \item{tsne_dim}{
When \code{tsne = TRUE}, the number of dimensions for the tSNE multidimensional scaling plots. This is passed to the \code{dims} argument of \code{\link[Rtsne:Rtsne]{Rtsne::Rtsne}}. Default is 2.
}
  \item{tsne_theta}{
When \code{tsne = TRUE}, a parameter controlling the speed/accuracy trade-off (increase for faster but less accurate results). This is passed to the \code{theta} argument of \code{\link[Rtsne:Rtsne]{Rtsne::Rtsne}}. Default is 0 for exact tSNE.
}
  \item{\dots}{
For \code{make_clusters()}, other arguments passed to \code{\link[Rtsne:Rtsne]{Rtsne::Rtsne}} when \code{tsne = TRUE}.
}
  \item{\dots}{
For \code{plot()}, when plotting a \code{cluster_df} object, other arguments passed to \code{\link[ggrepel:geom_label_repel]{ggrepel::geom_text_repel}} to control display of the observation labels.
}
  \item{x}{
A \code{cluster_df} object; the output of a call to \code{make_clusters()}.
}
  \item{seed}{
The seed used to control the placement of the labels and the jittering of the points. Jittering only occurs when \code{tsne = FALSE} in the call to \code{make_clusters()}. Using a non-\code{NA} seed ensure replicability across uses.
}
  \item{nrow}{
When \code{tsne = TRUE} in the call to \code{make_clusters()} and \code{tsne_dim} is greater than 2, the number of rows used to display the resulting 2-dimensional plots. Default is 1 for side-by-side plots.
}
}
\details{
\code{make_clusters} calls \code{\link[cluster:pam]{cluster::pam}} on the supplied Gower distance matrix with the specified number of clusters to determine cluster membership for each character. PAM is analogous to K-means, but it has its clusters centered around medoids instead of centered around centroids, which are less prone to the impact from outliers and heterogeneous cluster sizes. PAM also has the advantage over k-means of utilizing Gower distance matrices instead of Euclidean distance matrices only.

When \code{tsne = TRUE}, a Barnes-Hut t-distributed stochastic neighbor embedding is used to compute a multi-dimensional embedding of the distance matrix, coloring data points according to the PAM-defined clusters, as estimated by the function \code{make_clusters}. This graphic clustering allows users to independently test the quality of the chosen partitioning scheme from PAM, and can help in visualizing the resulting clusters. \code{\link[Rtsne:Rtsne]{Rtsne::Rtsne}} is used to do this. The resulting dimensions will be included in the output; see Value below.

\code{plot()} plots all morphological characters in a scatterplot with points colored based on cluster membership. When \code{tsne = TRUE} in the call to \code{make_clusters()}, the x- and y-axes will correspond to requested tSNE dimensions. With more than 2 dimensions, several plots will be produced, one for each pair of tSNE dimensions. These are displayed together using \code{\link[patchwork:plot_layout]{patchwork::plot_layout}}. When \code{tsne = FALSE}, the points will be arrange horizontally by cluster membership and randomly placed vertically.
}

\value{
A data frame, inheriting from class \code{"cluster_df"}, with a row for each character with its number (\code{character_number}) and cluster membership (\code{cluster}). When \code{tsne = TRUE}, additional columns  will be included, one for each requested tSNE dimension, labeled \code{tSNE_Dim1}, \code{tSNE_Dim2}, etc., containing the values on the dimensions computed using \code{Rtsne()}.

The \code{pam} fit resulting from \code{cluster::pam} is returned in the \code{"pam.fit"} attribute of the outut object.
}

\note{
When using \code{plot()} on a \code{cluster_df} object, warnings may appear from \code{ggrepel} saying something along the lines of "unlabeled data points (too many overlaps). Consider increasing max.overlaps". See \code{\link[ggrepel:geom_label_repel]{ggrepel::geom_text_repel}} for details; the \code{max.overlaps} argument can be supplied to \code{plot()} to increase the maximum number of element overlap in the plot. Alternatively, users can increase the size of the plot when exporting it, as it will increase the plot area and reduce the number of elements overlap. This warning can generally be ignored, though.
}

\seealso{
\code{\link{get_gower_dist}}, \code{\link{get_sil_widths}}

\code{\link[cluster:pam]{cluster::pam}}, \code{\link[Rtsne:Rtsne]{Rtsne::Rtsne}}
}
\examples{
data("characters")

#Reading example file as categorical data
Dmatrix <- get_gower_dist(characters)

sil_widths <- get_sil_widths(Dmatrix, max.k = 7)

sil_widths
# 3 clusters yields the highest silhouette width

#Create clusters with PAM under k=3 partitions
cluster_df <- make_clusters(Dmatrix, k = 3)

#Simple plot of clusters
plot(cluster_df, seed = 12345)

#Create clusters with PAM under k=3 partitions and perform tSNE (default = 2 dimensions)
cluster_df_tsne1 <- make_clusters(Dmatrix, k = 3, tsne = TRUE)

#Plot clusters
plot(cluster_df_tsne1)

#Create clusters with PAM under k=3 partitions and perform tSNE (3 dimensions)
cluster_df_tsne2 <- make_clusters(Dmatrix, k = 3, tsne = TRUE, tsne_dim = 3)

#Plot clusters, plots divided into 2 rows, and increasing overlap of text labels (default = 10)
plot(cluster_df_tsne2, nrow = 2, max.overlaps=20)
}

