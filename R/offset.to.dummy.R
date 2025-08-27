#' Convert trees produced by a BEAST2 FBD analysis with offset to trees with correct ages
#'
#' This method adds a dummy tip at the present (t = 0) to fully extinct trees with offsets,
#' in order to have correct ages (otherwise the most recent tip is assumed to be at 0).
#' This is a workaround to get the proper ages of the trees into other tools such as TreeAnnotator. `offset.to.dummy()` drops metadata on the tips, while `offset.to.dummy.metadata()` retains it (and is therefore a bit slower).
#'
#' @param trees.file path to BEAST2 output file containing posterior trees
#' @param log.file path to BEAST2 trace log file containing offset values
#' @param output.file path to file to write converted trees. If `NULL` (default), trees are simply returned.
#' @param dummy.name name of the added dummy tip, default `"dummy"`.
#'
#' @returns A list of converted trees (as treedata)
#'
#' @examples
#' # Convert trees with offset to trees with dummy tip
#' trees_file <- system.file("extdata", "ex_offset.trees", package = "EvoPhylo")
#' log_file <- system.file("extdata", "ex_offset.log", package = "EvoPhylo")
#' converted_trees <- offset.to.dummy.metadata(trees_file, log_file)
#'
#' # Do something with the converted trees - for instance, calculate the MCC summary tree
#' # Then remove the dummy tip from the MCC tree
#' final_tree <- drop.dummy.beast(system.file("extdata", "ex_offset.MCC.tre", package = "EvoPhylo"))

#' @export
offset.to.dummy <- function(trees.file, log.file, output.file = NULL,
                            dummy.name = "dummy") {
  trees <- ape::read.nexus(trees.file)
  log <- read.table(log.file, header = TRUE)
  log <- log$offset

  for (i in seq_len(length(trees))) {
    trees[[i]]$offset <- log[i]
  }

  presenttrees <- lapply(trees, function(t) {
    offset.to.dummy.phylo(t, dummy.name = dummy.name)
  })

  if (!is.null(output.file)) {
    ape::write.nexus(presenttrees, file = output.file)
  }

  presenttrees
}

#' @export
#' @rdname offset.to.dummy
offset.to.dummy.metadata <- function(trees.file, log.file, output.file = NULL,
                                     dummy.name = "dummy") {
  trees <- treeio::read.beast(trees.file)
  log <- read.table(log.file, header = TRUE)
  log <- log$offset

  print(stats::median(log))

  for (i in seq_along(trees)) {
    trees[[i]]@phylo$offset <- log[i]
  }

  presenttrees <- lapply(trees, function(tmp) {
    ntips <- length(tmp@phylo$tip.label)
    root <- ntips + 2

    t <- offset.to.dummy.phylo(tmp@phylo, dummy.name = dummy.name)

    tmp@data$node <- as.numeric(tmp@data$node)
    tmp@data$node[which(tmp@data$node > ntips)] <- tmp@data$node[which(tmp@data$node > ntips)]+1
    tmp@data$node[which(tmp@data$node >= root)] <- tmp@data$node[which(tmp@data$node >= root)]+1
    tmp@data$node <- as.character(tmp@data$node)
    metas <- colnames(tmp@data)
    metas <- metas[metas != "node"]
    new.data <- list(node = as.character(ntips + 1))
    for (m in metas) new.data[[m]] <- 0
    tmp@data <- rbind(tmp@data, new.data)
    tmp@phylo <- t
    tmp
  })

  if (!is.null(output.file)) {
    write.beast.treedata(presenttrees, file = output.file)
  }

  presenttrees
}


# helper function for offset.to.dummy and offset.to.dummy.metadata
# adds dummy tip to a phylo object
offset.to.dummy.phylo <- function(t, dummy.name = "dummy") {
  ntips <- length(t$tip.label)
  totn <- ntips + t$Nnode
  times <- ape::node.depth.edgelength(t)
  root_time <- max(times) + t$offset
  root <- ntips + 2

  t$edge[which(t$edge > ntips)] <- t$edge[which(t$edge > ntips)]+1
  t$edge[which(t$edge >= root)] <- t$edge[which(t$edge >= root)]+1
  t$edge <- rbind(c(root, root+1), c(root, ntips+1),t$edge)
  t$edge.length <- c(0.1 , root_time, t$edge.length)
  t$tip.label <- c(t$tip.label, dummy.name)
  t$Nnode <- t$Nnode + 1

  ntips <- length(t$tip.label)
  totn <- ntips + t$Nnode

  t
}
