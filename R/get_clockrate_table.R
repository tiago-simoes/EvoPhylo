#' Extract evolutionary rates from Bayesian clock trees produced by BEAST2 or Mr. Bayes
#'
#' Extract evolutionary rate summary statistics for each node from a Bayesian clock summary tree and stores them in a data frame. BEAST2 stores the rates for each clock in a separate file; all such trees need to be loaded using [treeio::read.beast()].
#'
#' @name get_clockrate_table
#'
#' @param ... `treedata` objects containing the summary trees from BEAST2 with associated data on the rates for each separate clock.
#' @param tree An S4 class object of type `treedata`; a Bayesian clock tree imported using [treeio::read.mrbayes()] for Mr. Bayes summary trees.
#' @param summary summary metric used for the rates. Currently supported: `"mean"` or `"median"`, default `"median"`.
#' @param drop_dummy if not `NULL`, will drop the dummy extant tip with the given label from the summary trees prior to extracting the clock rates (when present). Default is `NULL`.
#'
#' @returns
#' A data frame with a column containing the node identifier (`node`) and one column containing the clock rates for each tree provided, in the same order as the trees.
#'
#' @seealso
#' [clockrate_summary()] for summarizing and examining properties of the resulting rate table. Note that clade membership for each node must be customized (manually added) before these functions can be used, since this is tree and dataset dependent.
#'
#' See `vignette("rates-selection")` for how to use this function as part of an analysis pipeline.
#'
#' @examples
#' #Import all clock summary trees produced by BEAST2 from your local directory
#' \dontrun{
#' tree_clock1 <- treeio::read.beast("tree_file_clock1.tre")
#' tree_clock2 <- treeio::read.beast("tree_file_clock2.tre")
#' }
#'
#' #Or use the example BEAST2 multiple clock trees that accompany EvoPhylo.
#' data(tree_clock1)
#' data(tree_clock2)
#'
#' # obtain the rate table from BEAST2 trees
#' rate_table <- get_clockrate_table_BEAST2(tree_clock1, tree_clock2, summary = "mean")
#'
#' head(rate_table)
#'
#' #Import summary tree with three clock partitions produced by
#' #Mr. Bayes (.t or .tre files) from your local directory
#' \dontrun{
#'   tree3p <- treeio::read.mrbayes("Tree3p.t")
#' }
#'
#' #Or use the example Mr.Bayes multi-clock tree file (`tree3p`)
#' data("tree3p")
#'
#' # obtain the rate table from MrBayes tree
#' rate_table <- get_clockrate_table_MrBayes(tree3p)
#'
#' head(rate_table)

#' @export
#' @rdname get_clockrate_table
get_clockrate_table_BEAST2 <- function(..., summary = "median", drop_dummy = NULL) {
  trees <- list(...)

  if (length(trees) == 0L) {
    stop("No trees provided")
  }

  summary <- match.arg(summary, c("mean", "median"))
  name <- if (summary == "mean") "rate" else "rate_median"

  if (!is.null(drop_dummy)) {
    trees <- lapply(trees, treeio::drop.tip, drop_dummy)
  }

  rate_table <- data.frame(nodes = as.integer(trees[[1]]@data$node))

  for (tr in trees) {
    data <- tr@data[match(rate_table$nodes, as.integer(tr@data$node)), ] #get rates in same order as 1st column
    rate_table <- cbind(rate_table, data[[name]])
  }

  colnames(rate_table)[2:ncol(rate_table)] <- paste0("rates", 1:(ncol(rate_table) - 1))
  row.names(rate_table) <- NULL

  return(rate_table)
}

#' @export
#' @rdname get_clockrate_table
get_clockrate_table_MrBayes <- function(tree, summary = "median", drop_dummy = NULL) {

  if (!is.null(drop_dummy)) {
    tree <- treeio::drop.tip(tree, drop_dummy)
  }

  nodes <- as.integer(tree@data$node)

  p <- unglue::unglue_data(names(tree@data), "rate<model>rlens<clock>_<summary>",
                           open = "<", close = ">")
  rownames(p) <- names(tree@data)

  p <- p[rowSums(is.na(p)) < ncol(p),,drop=FALSE]

  p$clock <- gsub("\\{|\\}", "", p$clock)

  summary <- match.arg(summary, c("mean", "median"))

  rates <- rownames(p)[p$summary == summary]

  rate_table <- setNames(data.frame(nodes, tree@data[rates]),
                         c("nodes", paste0("rates", p[rates, "clock"])))

  for (i in seq_len(ncol(rate_table))[-1]) {
    rate_table[[i]] <- as.numeric(rate_table[[i]])
  }

  rate_table
}
