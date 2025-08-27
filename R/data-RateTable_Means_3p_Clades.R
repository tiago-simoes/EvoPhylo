#' Mean clock rates by node and clade (3 clock partitions)
#'
#' A data set containing the mean clock rates for a tree with 3 clock
#' partitions, such as the output of [get_clockrate_table_MrBayes()]
#' but with an additional "clade" column added, which is required for use in
#' [clockrate_summary()] and [clockrate_dens_plot()].
#'
#' `RateTable_Means_3p_Clades` was created by running `get_clockrate_table_MrBayes(tree3p)` and then adding a "clade" column. It can be produced by using the following procedure:
#'
#' \enumerate{
#' \item{Import tree file:
#'
#' \preformatted{data("tree3p")}
#' }
#' \item{Produce clock rate table with, for instance, mean rate values from each
#' branch in the tree:
#'
#' \preformatted{rate_table <- get_clockrate_table_MrBayes(tree3p, summary = "mean")
#' write.csv(rate_table, file = "rate_table.csv", row.names = FALSE)}
#' }
#' \item{Now, manually add clades using, e.g., Excel:
#'
#' 3.1. Manually edit rate_table.csv, adding a "clade" column. This introduces
#' customized clade names to individual nodes in the tree.
#'
#' 3.2. Save the edited rate table with a different name to differentiate from
#' the original output (e.g., rate_table_clades_means.csv).
#' }
#' \item{
#' Read the file back in:
#'
#' \preformatted{RateTable_Means_3p_Clades <- read.csv("rate_table_clades_means.csv")
#'
#' head(RateTable_Means_3p_Clades)}
#' }
#' }
#' @name RateTable_Means_3p_Clades
#' @docType data
#' @usage data("RateTable_Means_3p_Clades")
#' @format A data frame with 79 observations on the following 5 variables.
#' \describe{
#' \item{`clade`}{A character vector containing the clade
#' names for each corresponding node}
#' \item{`nodes`}{A numeric vector for
#' the node numbers in the summary tree}
#' \item{`rates1`}{A numeric vector
#' containing the mean posterior clock rate for each node for the first
#' partition}
#' \item{`rates2`}{A numeric vector containing the mean
#' posterior clock rate for each node for the second partition}
#' \item{`rates3`}{A numeric vector containing the mean posterior clock
#' rate for each node for the third partition} }
#' @seealso [`tree3p`] for the tree from which the clock rates were
#' extracted.
#'
#' [get_clockrate_table_MrBayes()] for extracting a clock rate table
#' from a tree.
#'
#' [clockrate_summary()], [clockrate_dens_plot()], and
#' [clockrate_reg_plot()] for examples of using a clockrate table.
#' @keywords datasets
NULL
