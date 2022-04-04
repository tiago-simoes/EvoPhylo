\name{get_sil_widths}
\alias{get_sil_widths}
\alias{plot.sil_width_df}

\title{
Compute Silhouette Widths for Various Numbers of Clusters
}
\description{
Computes silhouette widths for several possible numbers of clusters \code{k}, which determines how well an object falls within their cluster compared to other clusters. The best number of clusters \code{k} is the one with the highest silhouette width.
}
\usage{
get_sil_widths(dist_mat, max.k = 10)

\method{plot}{sil_width_df}(x, ...)
}
\arguments{
  \item{dist_mat}{
a Gower distance matrix, the output of a call to \code{\link{get_gower_dist}}.
}
  \item{max.k}{
the maximum number of clusters to search across.
}
  \item{x}{
a \code{sil_width_df} object; the output of a call to \code{get_sil_widths()}.
}
  \item{\dots}{
further arguments passed to \code{\link[ggplot2:geom_path]{ggplot2::geom_line}} to control the appearance of the plot.
}
}
\details{
\code{get_sil_widths} calls \code{\link[cluster:pam]{cluster::pam}} on the supplied Gower distance matrix with each number of clusters up to \code{mak.k} and stores the average silhouette widths across the clustered characters. When \code{plot = TRUE}, a plot of the sillhouette widths against the number of clusters is produced, though this can also be produced seperately on the resulting data frame using \code{plot.sil_width_df()}. The number of clusters with the greatest silhouette width should be selected for use in the final clustering specification. See Examples.
}
\value{
For \code{get_sil_widths()}, a data frame, inheriting from class \code{"sil_width_df"}, with two columns: \code{k} is the number of clusters, and \code{sil_width} is the silhouette widths for each number of clusters. If \code{plot = TRUE}, the output is returned invisibly.

For \code{plot.sil_widths_df()}, a \code{ggplot} object that can be manipulated using \pkg{ggplot2} syntax (e.g., to change the \code{theme} or labels).
}

\seealso{
\code{\link{get_gower_dist}}, \code{\link[cluster:pam]{cluster::pam}}
}
\examples{
data("characters")

#Reading example file as categorical data
Dmatrix <- get_gower_dist(characters)

#Get silhouette widths for k=7
sw <- get_sil_widths(Dmatrix, max.k = 7)

sw

plot(sw, color = "red", size =2)
}
