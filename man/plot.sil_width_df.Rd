\name{plot.sil_width_df}
\alias{plot.sil_width_df}

\title{
Plot Silhouette Widths against Number of Clusters
}
\description{
Makes a plot of silhouette widths for various number of clusters for use in selecting the clustering specification to use.
}
\usage{
\method{plot}{sil_width_df}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a \code{sil_width_df} object; the output of a call to \code{\link{get_sil_widths}}.
}
  \item{\dots}{
ignored.
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
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{get_sil_widths}}
}
\examples{

}
