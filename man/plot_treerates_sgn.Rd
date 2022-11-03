\name{plot_treerates_sgn}
\alias{plot_treerates_sgn}

\title{
Plot Bayesian evolutionary tree with rate thresholds for selection mode
}
\description{
Plots the summary Bayesian evolutionary tree with branches, according to user-defined thresholds (in units of standard deviations) used to infer the strength and mode of selection.
}
\usage{
plot_treerates_sgn(type = c("MrBayes", "BEAST2"),
                  tree, posterior,
                  trans = c("none", "log", "log10"),
                  summary = "mean", drop.dummyextant = TRUE,
                  clock = 1, threshold = c("1 SD", "2 SD"),
                  low = "blue", mid = "gray90", high = "red",
                  branch_size = 2, tip_size = 2,
                  xlim = NULL, nbreaks = 10, geo_size = list(2, 3),
                  geo_skip = c("Quaternary", "Holocene", "Late Pleistocene"))
}
\arguments{
  \item{type}{
Whether to use data output from "Mr.Bayes" or "BEAST2".
}
  \item{tree}{
A \code{tidytree} object; the output of a call to \code{\link[treeio:read.beast]{treeio::read.beast}}. Summary trees from Mr. Bayes will include branch specific rates for all clock partitions, and the partition to be plotted will be specified using the "clock" argument. On the other hand, BEAST2 will output one separate summary tree file for each clock partition. For the latter, the tree file for the partition of interest should be provided for plotting.
}
  \item{posterior}{
A data frame of posterior parameter estimates (log file). From Mr.Bayes, it includes a "clockrate" column indicating the mean (background) clock rate estimate for each generation that will be used for pairwise t-tests. Such data frame can be imported using \code{\link{combine_log}} (no need to reshape from wide to long). See the \code{\link{posterior1p}} or \code{\link{posterior3p}} datasets for an examples of how the input file should look.
     From BEAST2, it will include at least one "rate<filename>.mean" column indicating the mean (background) clock rate estimate for each generation. If there are "P" unlinked clock partitions in BEAST2, there will be P x "rate<filename>.mean" columns (one for each partition) in the posterior log file.
}
  \item{trans}{
Type of data transformation to perform on background rates extracted from the posterior log file from Mr. Bayes or BEAST2. Options include "none" (if rates are normally distributed), natural log transformation "log", and log of base 10 transformation "log10". The necessity of using data transformation can be tested using the function \code{\link{plot_back_rates}}.
}
  \item{summary}{
Only when using Mr. Bayes trees. The rate summary stats chosen to calculate selection mode. Only rates "mean" and "median" are allowed. Default is "mean".
}
  \item{drop.dummyextant}{
\code{logical}; Only when using Mr. Bayes trees. Whether to drop the "Dummyextant" tip (if present) from the tree before plotting the tree. Default is \code{TRUE}.
}
  \item{clock}{
The clock partition number to calculate selection mode. Ignored if only one clock is available.
}
  \item{threshold}{
A vector of threshold values. Default is to display thresholds of ±1 relative standard deviation (SD) of the relative posterior clock rates. Should be specified as a number of standard deviations (e.g., \code{"1 SD"}) or the confidence level for a confidence interal around the mean relative posterior clockrate (e.g., \code{"95\%"}). Multiple values are allowed to produce a plot with multiple thresholds. Set to \code{NULL} to omit thresholds.
}
  \item{low, mid, high}{
Colors passed to \code{\link{scale_color_steps2}} to control the colors of the branches based on which thresholds are exceeded. When no thresholds are supplied, use \code{mid} to control the color of the tree.
}
  \item{branch_size}{
The thickness of the lines that form the tree.
}
  \item{tip_size}{
The font size for the tips of the tree.
}
  \item{xlim}{
The x-axis limits. Should be two negative numbers (though the axis labels will be in absolute value, i.e., Ma).
}
  \item{nbreaks}{
The number of interval breaks in the geological timescale.
}
  \item{geo_size}{
The font size for the labels in the geological scale. The first value in \code{list()} is the font size for geological epochs and the second value is for geological periods. Passed directly to the \code{size} argument of \code{\link[deeptime:coord_geo]{deeptime::coord_geo}}.
}
 \item{geo_skip}{
A vector of interval names indicating which intervals should not be labeled. Passed directly to the \code{skip} argument of \code{\link[deeptime:coord_geo]{deeptime::coord_geo}}.
}
}
\details{
Plots the phylogentic tree contained in \code{tree} using \code{\link[ggtree:ggtree]{ggtree::ggtree}}. Branches undergoing accelerating evolutionary rates (e.g., >\code{"1 SD"}, \code{"3 SD"}, or \code{"5 SD"} relative to the background rate) for each morphological clock partition suggest directional (or positive) selection for that morphological partition in that branch of the tree. Branches undergoing decelerating evolutionary rates (e.g., <\code{"1 SD"}, \code{"3 SD"}, or \code{"5 SD"} relative to the background rate) for each morphological clock partition suggest stabilizing selection for that morphological partition in that branch of the tree. For details on rationale, see Simões & Pierce (2021).

Please double check that the distribution of background rates (mean rates for the tree) sampled from the posterior follow the assumptions of a normal distribution (e.g., check for normality of distribution in Tracer). Otherwise, displayed results may not have a valid interpretation.
}

\value{
A \code{ggtree} object, which inherits from \code{ggplot}.
}
\references{
Simões, T. R. and S. E. Pierce (2021). Sustained High Rates of Morphological Evolution During the Rise of Tetrapods. \emph{Nature Ecology & Evolution} 5: 1403–1414.
}
\seealso{
\code{vignette("rates-selection")} for the use of this function as part of an analysis pipeline.
\code{\link[ggtree:ggtree]{ggtree::ggtree}}, \code{\link[deeptime:coord_geo]{deeptime::coord_geo}}
}
\examples{
# See vignette("rates-selection") for how to use this
# function as part of an analysis pipeline

## MrBayes example
# Load example tree and posterior
data("tree3p")
data("posterior3p")

plot_treerates_sgn(
  type = "MrBayes",
  tree3p, posterior3p,          #MrBayes tree file with data for all partitions
  trans = "none",
  summary = "mean",             #MrBayes specific argument
  drop.dummyextant = TRUE,      #MrBayes specific argument
  clock = 1,                           #Show rates for clock partition 1
  threshold = c("1 SD", "3 SD"),       #sets background rate threshold for selection mode
  branch_size = 1.5, tip_size = 3,                          #sets size for tree elements
  xlim = c(-450, -260), nbreaks = 8, geo_size = list(3, 3)) #sets limits and breaks for geoscale

\dontrun{
## BEAST2 example
tree_clock1 <- system.file("extdata", "Penguins_MCC_morpho_part1", package = "EvoPhylo")
tree_clock1 <- treeio::read.beast(tree_clock1)
posterior <- system.file("extdata", "Penguins_log.log", package = "EvoPhylo")
posterior <- read.table(posterior, header = TRUE)

plot_treerates_sgn(
  type = "BEAST2",
  tree_clock1, posterior,                 #BEAST2 tree file with data for partition 1
  trans = "log10",
  clock = 1,                              #Show rates for clock partition 1
  threshold = c("1 SD", "3 SD"),          #sets background rate threshold for selection mode
  branch_size = 1.5, tip_size = 3,                        #sets size for tree elements
  xlim = c(-70, 30), nbreaks = 8, geo_size = list(3, 3))  #sets limits and breaks for geoscale
}}

