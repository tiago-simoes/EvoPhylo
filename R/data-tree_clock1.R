#' BEAST2 phylogenetic tree with clock rates from partition 1
#'
#' A clock Bayesian phylogenetic tree with clock rates from a single clock
#' partition (partition 1 here), imported as an S4 class object using
#' [treeio::read.beast()].
#'
#' This example tree file was produced by analyzing the data set with a single
#' morphological partition from
#'
#' @name tree_clock1
#' @docType data
#' @usage data("tree_clock1")
#' @format A `tidytree` object.
#' @seealso
#'
#' [`tree_clock2`] for another BEAST2 tree object with clock rates
#' from partition 2 for this same dataset.
#'
#' [`tree3p`] for another tree object with 3 clock partitions from
#' Mr.Bayes.
#'
#' [`tree1p`] for another tree object with a single clock from
#' Mr.Bayes.
#'
#' [get_clockrate_table_BEAST2()] for extracting the posterior clock
#' rates from BEAST2 tree objects.
#' @keywords datasets
NULL
