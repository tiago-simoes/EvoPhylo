#' Remove dummy tip from summary trees, accounting for metadata on the
#' tips
#'
#' This method is designed to remove the dummy tip added on offset trees once
#' postprocessing is complete (for instance once the summary tree has been
#' built using TreeAnnotator). `drop.dummy.beast()` is designed to remove the dummy tip added to a dataset before
#' running with BEAST2, and `drop.dummy.mb()` is for running with Mr. Bayes.
#'
#' @name drop.dummy
#'
#' @param tree.file path to file containing the tree with dummy tip
#' @param output.file path to file to write converted tree. If `NULL`
#' (default), the tree is simply returned.
#' @param dummy.name name of the added dummy tip, default `dummy`.
#' @param convert.heights,convert.ages whether height or age metadata should be converted to height
#' - offset or age - offset (required to plot e.g. HPD intervals correctly). Default `TRUE`.
#'
#' @returns
#' A list of `tree` converted tree (as treedata) and `offset`
#' height/age of the youngest tip in the final tree
#'
#' @examples
#' # Analyze the trees with dummy tips - for instance, calculate the MCC summary tree
#' # Then remove the dummy tip from the MCC tree
#' final_tree1 <- drop.dummy.beast(system.file("extdata", "ex_offset.MCC.tre", package = "EvoPhylo"))
#'
#' # Remove the dummy tip from the summary tree
#' final_tree2 <- drop.dummy.mb(system.file("extdata", "tree_mb_dummy.tre", package = "EvoPhylo"))

#' @export
#' @rdname drop.dummy
drop.dummy.beast <- function(tree.file, output.file = NULL, dummy.name = "dummy", convert.heights = TRUE) {
  tmp <- treeio::read.beast(tree.file)
  tip <- which(tmp@phylo$tip.label == dummy.name)
  tipn <- which(tmp@data$node == as.character(tip))
  node <- tmp@phylo$edge[which(tmp@phylo$edge[,2] == tip),1]
  noden <- which(tmp@data$node == as.character(node))
  offset <- min(tmp@data$height_median[which(as.numeric(tmp@data$node) <= length(tmp@phylo$tip.label) &
                                               as.numeric(tmp@data$node) != tip)])

  tmp@phylo <- ape::drop.tip(tmp@phylo, tip)
  tmp@data <- tmp@data[-c(tipn,noden),]
  tmp@data$node <- as.numeric(tmp@data$node)
  tmp@data$node[(tmp@data$node) > tip] <- tmp@data$node[(tmp@data$node) > tip] - 1
  tmp@data$node[(tmp@data$node) > node] <- tmp@data$node[(tmp@data$node) > node] - 1
  tmp@data$node <- as.character(tmp@data$node)

  if (convert.heights) {
    for (m in c("height_0.95_HPD", "height_range", "height_median", "height")) {
      tmp@data[[m]] <- lapply(tmp@data[[m]], function(x) x - offset)
    }
  }

  if (!is.null(output.file)) {
    write.beast.treedata(list(tmp), file = output.file)
  }

  list(tree = tmp, offset = offset)
}

#' @export
#' @rdname drop.dummy
drop.dummy.mb <- function(tree.file, output.file = NULL, dummy.name = "dummy", convert.ages = TRUE) {
  tmp <- treeio::read.mrbayes(tree.file)

  tmp@data$`prob+-sd` <- as.factor(tmp@data$`prob+-sd`)
  for (i in seq_along(tmp@data)) {
    if (is.character(tmp@data[[i]])) {
      tmp@data[[i]] <- as.numeric(tmp@data[[i]])
    }
  }

  tip <- which(tmp@phylo$tip.label == dummy.name)
  tipn <- which(tmp@data$node == as.character(tip))
  node <- tmp@phylo$edge[which(tmp@phylo$edge[,2] == tip),1]
  noden <- which(tmp@data$node == as.character(node))
  offset <- min(tmp@data$age_median[which(as.numeric(tmp@data$node) <= length(tmp@phylo$tip.label) &
                                            as.numeric(tmp@data$node) != tip)])

  tmp@phylo <- ape::drop.tip(tmp@phylo, tip)
  tmp@data <- tmp@data[-c(tipn,noden),]
  tmp@data$node <- as.numeric(tmp@data$node)
  tmp@data$node[(tmp@data$node) > tip] <- tmp@data$node[(tmp@data$node) > tip] - 1
  tmp@data$node[(tmp@data$node) > node] <- tmp@data$node[(tmp@data$node) > node] - 1
  tmp@data$node <- as.character(tmp@data$node)

  if (convert.ages) {
    for (m in c("age_median", "age_mean")) {
      tmp@data[[m]] <- lapply(tmp@data[[m]], "-", offset)
    }

    tmp@data$age_0.95HPD <- lapply(tmp@data$age_0.95HPD,"-", offset)
  }

  if (!is.null(output.file)) {
    write.beast.treedata(list(tmp), file = output.file)
  }

  list(tree = tmp, offset = offset)
}
