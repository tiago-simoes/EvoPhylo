\name{f2}
\alias{f2}

\title{
Plots a tree with thresholds(?)
}
\description{
Plots trees with user-specified thresholds
}
\usage{
plot_treerates_sgn(tree, posterior.clockrate.all., clock = 1,
   summary = "mean", threshold = c("1 SD", "2 SD"),
   drop.dummyextent = TRUE,
   low = "blue", mid = "gray90", high = "red",
   xlim = NULL, geo_skip = c("Quaternary", "Holocene",
                             "Late Pleistocene"))
}
\arguments{
  \item{tree}{
a \code{tidytree} object; the output of a call to \code{\link[treeio:read.mrbayes]{treeio::read.mrbayes}}.
}
  \item{posterior.clockrate.all.}{
a vector of posterior clock rates (the \code{"clockrate.all."} column of a posterior object(?))
}
  \item{clock}{
the clock number to plot, (i.e., the number "rateTK02Brlens\{#\}_mean"). Ignored if only one clock is available.
}
  \item{summary}{
the rate summary to plot, (i.e., the value "rateTK02Brlens1_\{summary\}"). Only "mean" and "median" are allowed. Default is "mean".
}
  \item{threshold}{
a vector of threshold values. Default is to display thresholds of ±1 relative standard deviation of the posterior clockrates. Should be specified as a number of standard deviations (e.g., \code{"1 SD"}) or the confidence level for a confidence interal around the mean posterior clockrate (e.g., \code{"95\%"}). Multiple values are allowed to produce a plot with multiple thresholds.
}
  \item{drop.dummyextent}{
\code{logical}; whether to drop the "DummyExtent" tip of the tree prior to plotting the tree. Default is \code{TRUE}.
}
  \item{low, mid, high}{
colors passed to \code{\link{scale_color_steps2}} to control the colors of the branches based on which thresholds are exceeded.
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

}

