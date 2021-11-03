\name{plot.cluster_df}
\alias{plot.cluster_df}

\title{
Plot Cluster Membership
}
\description{
Plots cluster membership along two dimensions if present.
}
\usage{
\method{plot}{cluster_df}(x, seed = 100, ...)
}
\arguments{
  \item{x}{
a \code{cluster_df} object; the output of a call to \code{\link{make_clusters}}.
}
  \item{seed}{
the seed passed to \code{\link[ggplot2:position_jitter]{ggplot2::position_jitter}}. Can generally be ignored, but different values can produce slightly different plots. Ignored if \code{tsne = FALSE} in the original call to \code{\link{make_clusters}}.
}
  \item{\dots}{
additional arguments passed to \code{\link[ggrepel:geom_text_repel]{ggrepel::geom_text_repel}} to control label position.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A \code{ggplot} object that can be manipulated using \pkg{ggplot2} syntax (e.g., to change the \code{theme} or labels).
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
