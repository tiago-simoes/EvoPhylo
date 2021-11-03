\name{get_sil_widths}
\alias{get_sil_widths}
\title{
Compute Silhouette Widths for Various Numbers of Clusters
}
\description{
Computes silhouette widths for several possible numbers of clusters to help determine the best number, usually as that with the highest silhouette width.
}
\usage{
get_sil_widths(Dmatrix, max.k = 10, plot = TRUE)
}
\arguments{
  \item{Dmatrix}{
a Gower distance matrix, the output of a call to \code{\link{get_gower_dist}}.
}
  \item{max.k}{
the maximum number of clusters to search across.
}
  \item{plot}{
whether to plot the silhouette widths against the the number of clusters. If \code{TRUE}, \code{\link{plot.sil_width_df}} is called on the resulting output.
}
}
\details{
\code{get_sil_widths} calls \code{\link[cluster:pam]{cluster::pam}} on the supplied Gower distance matrix with each number of clusters up to \code{mak.k} and stores the average silhouette widths across the clustered characters. When \code{plot = TRUE}, a plot of the sillhouette widths against the number of clusters is produced, though this can also be produced seperately on the resulting data frame using \code{\link{plot.sil_width_df}}. The number of clusters with the greatest silhouette width should be selected for use in the final clustering specification. See Examples.
}
\value{
A data frame, inheriting from class \code{"sil_width_df"}, with two columns: \code{k} is the number of clusters, and \code{sil_width} is the silhouette widths for each number of clusters. If \code{plot = TRUE}, the output is returned invisibly.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{get_gower_dist}}, \code{\link{plot.sil_width_df}}, \code{\link[cluster:pam]{cluster::pam}}
}
\examples{

}
