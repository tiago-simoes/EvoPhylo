\name{plot_treerates_sgn}
\alias{plot_treerates_sgn}

\title{
TITLE
%Needs a title
}
\description{
%Needs a description
}
\usage{
plot_treerates_sgn(tree, posterior, clock = 1,
   summary = "mean", threshold = c("1 SD", "2 SD"),
   drop.dummyextent = TRUE,
   low = "blue", mid = "gray90", high = "red", size = 2,
   xlim = NULL, geo_skip = c("Quaternary", "Holocene",
                             "Late Pleistocene"))
}
\arguments{
  \item{tree}{
a \code{tidytree} object; the output of a call to \code{\link[treeio:read.mrbayes]{treeio::read.mrbayes}}.
}
  \item{posterior}{
a dataset of FBD parameter MCMC posterior samples imported using \code{\link{import_log}}. See the \code{\link{posterior1p}} dataset for an example of how the resulting file should look. The \code{clockrate} column from this dataset (divided by its mean to form relative clockrates) will be used along with the supplied thresholds to color the tree.
}
  \item{clock}{
the clock number to plot, (i.e., the number "rateTK02Brlens\{#\}_mean"). Ignored if only one clock is available.
}
  \item{summary}{
the rate summary to plot, (i.e., the value "rateTK02Brlens1_\{summary\}"). Only "mean" and "median" are allowed. Default is "mean".
}
  \item{threshold}{
a vector of threshold values. Default is to display thresholds of ±1 relative standard deviation of the relative posterior clockrates. Should be specified as a number of standard deviations (e.g., \code{"1 SD"}) or the confidence level for a confidence interal around the mean relative posterior clockrate (e.g., \code{"95\%"}). Multiple values are allowed to produce a plot with multiple thresholds. Set to \code{NULL} to omit thresholds.
}
  \item{drop.dummyextent}{
\code{logical}; whether to drop the "DummyExtent" tip of the tree prior to plotting the tree. Default is \code{TRUE}.
}
  \item{low, mid, high}{
colors passed to \code{\link{scale_color_steps2}} to control the colors of the branches based on which thresholds are exceeded. When no thresholds are supplied, use \code{mid} to control the color of the tree.
}
  \item{size}{
the thickness of the lines that form the tree.
}
  \item{xlim}{
the x-axis limits. Should be two negative numbers (though the axis labels will be in absolute value, i.e., Ma).
}
 \item{geo_skip}{
a vector of interval names indicating which intervals should not be labeled. Passed directly to the \code{skip} argument of \code{\link[deeptime:coord_geo]{deeptime::coord_geo}}.
}
}
\details{
Plots the phylogentic tree contained in \code{tree} using \code{\link[ggtree:ggtree]{ggtree::ggtree}}.
}
\value{
A \code{ggtree} object, which inherits from \code{ggplot}.
}

\seealso{
\code{\link[ggtree:ggtree]{ggtree::ggtree}}, \code{\link[deeptime:coord_geo]{deeptime::coord_geo}}
}
\examples{
data("tree3p")
data("posterior3p")

plot_treerates_sgn(tree3p, posterior3p, clock = 1,
                   summary = "mean",
                   threshold = c("1 SD", "95\%"))
}

