#' Convert trees produced by a BEAST2 FBD analysis with offset to trees with correct ages,
#' accounting for possible metadata on the tips.
#'
#' This method adds a dummy tip at the present (t = 0) to fully extinct trees with offsets,
#' in order to have correct ages (otherwise the most recent tip is assumed to be at 0).
#' This is a workaround to get the proper ages of the trees into other tools such as TreeAnnotator.
#'
#' @param trees.file path to BEAST2 output file containing posterior trees
#' @param log.file path to BEAST2 trace log file containing offset values
#' @param output.file path to file to write converted trees. If \code{NULL} (default), trees are simply returned.
#' @param dummy.name name of the added dummy tip, default \code{dummy}.
#'
#' @return list of converted trees (as treedata)
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
#'
#' @export
#' @importFrom stats median
#' @seealso [offset.to.dummy()] (faster version discarding metadata)
#' @md
offset.to.dummy.metadata <- function(trees.file, log.file, output.file = NULL,
                                     dummy.name = "dummy") {
  trees <- treeio::read.beast(trees.file)
  log <- read.table(log.file, header = TRUE)
  log <- log$offset
  print(median(log))
  
  for(i in 1:length(trees)) trees[[i]]@phylo$offset <- log[i]
  
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
    for(m in metas) new.data[[m]] <- 0
    tmp@data <- rbind(tmp@data, new.data)
    
    tmp@phylo <- t
    return(tmp)
  })
  
  if(!is.null(output.file)) write.beast.treedata(presenttrees, file = output.file)
  presenttrees
}

#' Convert trees produced by a BEAST2 FBD analysis with offset to trees with correct ages.
#'
#' This method adds a dummy tip at the present (t = 0) to fully extinct trees with offsets,
#' in order to have correct ages (otherwise the most recent tip is assumed to be at 0).
#' This is a workaround to get the proper ages of the trees into other tools such as TreeAnnotator.
#'
#' **NB:** Any metadata present on the tips will be discarded. If you want to keep metadata (such as clock rate values),
#' use \code{offset.to.dummy.metadata} instead.
#'
#' @param trees.file path to BEAST2 output file containing posterior trees
#' @param log.file path to BEAST2 trace log file containing offset values
#' @param output.file path to file to write converted trees. If \code{NULL} (default), trees are simply returned.
#' @param dummy.name name of the added dummy tip, default \code{dummy}.
#'
#' @return list of converted trees (as treedata)
#'
#' @seealso [offset.to.dummy.metadata()] (slower version, keeping metadata)
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
#'
#' @export
#' @importFrom stats median
#' @md
offset.to.dummy <- function(trees.file, log.file, output.file = NULL,
                            dummy.name = "dummy") {
  trees <- ape::read.nexus(trees.file)
  log <- read.table(log.file, header = TRUE)
  log <- log$offset
  
  for(i in 1:length(trees)) trees[[i]]$offset <- log[i]
  
  presenttrees <- lapply(trees, function(t) {
    offset.to.dummy.phylo(t, dummy.name = dummy.name)
  })
  
  if(!is.null(output.file)) ape::write.nexus(presenttrees, file = output.file)
  presenttrees
}


#' Remove dummy tip from beast summary trees, accounting for metadata on the tips
#'
#' This method is designed to remove the dummy tip added on offset trees once postprocessing is
#' complete (for instance once the summary tree has been built using TreeAnnotator).
#'
#' @param tree.file path to file containing the tree with dummy tip
#' @param output.file path to file to write converted tree. If \code{NULL} (default), the tree is simply returned.
#' @param dummy.name name of the added dummy tip, default \code{dummy}.
#' @param convert.heights whether height metadata should be converted to height - offset (required to plot e.g. HPD intervals correctly). Default TRUE.
#'
#' @return list of \code{tree} converted tree (as treedata) ; and \code{offset} age of the youngest tip in the final tree
#'
#' @seealso [drop.dummy.mb()] for the same function using summary trees with a "dummy" extant from Mr. Bayes
#'
#' @examples
#' # Analyze the trees with dummy tips - for instance, calculate the MCC summary tree
#' # Then remove the dummy tip from the MCC tree
#' final_tree <- drop.dummy.beast(system.file("extdata", "ex_offset.MCC.tre", package = "EvoPhylo"))
#'
#' @export
#'
#' @md
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
  
  if(convert.heights) {
    for(m in c("height_0.95_HPD", "height_range", "height_median", "height")) {
      tmp@data[[m]] <- lapply(tmp@data[[m]], function(x) x - offset)
    }
  }
  
  if(!is.null(output.file)) write.beast.treedata(list(tmp), file = output.file)
  
  list(tree = tmp, offset = offset)
}



#' Remove dummy tip from Mr. Bayes summary trees, accounting for metadata on the tips
#'
#' This method is designed to remove the dummy tip added to a dataset before running with Mr. Bayes.
#'
#' @param tree.file path to file containing the tree with dummy tip
#' @param output.file path to file to write converted tree. If \code{NULL} (default), the tree is simply returned.
#' @param dummy.name name of the added dummy tip, default \code{dummy}.
#' @param convert.ages whether height metadata should be converted to height - offset (required to plot e.g. HPD intervals correctly). Default TRUE.
#'
#' @return list of \code{tree} converted tree (as treedata) ; and \code{offset} age of the youngest tip in the final tree
#'
#' @seealso [drop.dummy.beast()] for the same function using summary trees with a "dummy" extant from BEAST2
#'
#' @examples
#' # Remove the dummy tip from the summary tree
#' final_tree <- drop.dummy.mb(system.file("extdata", "tree_mb_dummy.tre", package = "EvoPhylo"))
#'
#' @export
#'
#' @md
drop.dummy.mb <- function(tree.file, output.file = NULL, dummy.name = "dummy", convert.ages = TRUE) {
  
  tmp <- treeio::read.mrbayes(tree.file)
  
  tmp@data$`prob+-sd` <- as.factor(tmp@data$`prob+-sd`)
  tmp@data <- dplyr::mutate_if(tmp@data, is.character,as.numeric)
  
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
  
  if(convert.ages) {
    for(m in c("age_median", "age_mean")) {
      tmp@data[[m]] <- lapply(tmp@data[[m]], "-", offset)
    }
    tmp@data$age_0.95HPD <- lapply(tmp@data$age_0.95HPD,"-", offset)
  }
  
  if(!is.null(output.file)) write.beast.treedata(list(tmp), file = output.file)
  
  list(tree = tmp, offset = offset)
}

# helper function for offset.to.dummy and offset.to.dummy.metadata
# adds dummy tip to a phylo object
offset.to.dummy.phylo <- function(t, dummy.name = "dummy") {
  ntips <- length(t$tip.label)
  totn <- ntips + t$Nnode
  
  times <- ape::node.depth.edgelength(t)
  root_time <- max(times) + t$offset
  
  root <- ntips+2
  t$edge[which(t$edge > ntips)] <- t$edge[which(t$edge > ntips)]+1
  t$edge[which(t$edge >= root)] <- t$edge[which(t$edge >= root)]+1
  
  t$edge <- rbind(c(root, root+1), c(root, ntips+1),t$edge)
  t$edge.length <- c(0.1 , root_time, t$edge.length)
  t$tip.label <- c(t$tip.label, dummy.name)
  t$Nnode <- t$Nnode + 1
  
  ntips <- length(t$tip.label)
  totn <- ntips + t$Nnode
  
  return(t)
}

# adapted from treeio to handle list of trees instead of single trees
write.beast.treedata <- function(treedata, file = "",
                                 translate = TRUE, tree.name = "STATE"){
  cat("#NEXUS\n", file = file)
  cat(paste("[R-package treeio, ", date(), "]\n\n", sep = ""),
      file = file, append = TRUE)
  N <- treeio::Ntip(treedata[[1]])
  
  obj <- lapply(treedata, treeio::as.phylo)
  ntree <- length(obj)
  cat("BEGIN TAXA;\n", file = file, append = TRUE)
  cat(paste("\tDIMENSIONS NTAX = ", N, ";\n", sep = ""), file = file,
      append = TRUE)
  cat("\tTAXLABELS\n", file = file, append = TRUE)
  cat(paste("\t\t", obj[[1]]$tip.label, sep = ""), sep = "\n",
      file = file, append = TRUE)
  cat("\t;\n", file = file, append = TRUE)
  cat("END;\n", file = file, append = TRUE)
  cat("BEGIN TREES;\n", file = file, append = TRUE)
  
  if (translate) {
    cat("\tTRANSLATE\n", file = file, append = TRUE)
    obj <- ape::.compressTipLabel(obj)
    X <- paste("\t\t", 1:N, "\t", attr(obj, "TipLabel"),
               ",", sep = "")
    X[length(X)] <- gsub(",", "", X[length(X)])
    cat(X, file = file, append = TRUE, sep = "\n")
    cat("\t;\n", file = file, append = TRUE)
    class(obj) <- NULL
    for (i in 1:ntree) obj[[i]]$tip.label <- as.character(1:N)
  }
  else {
    if (is.null(attr(obj, "TipLabel"))) {
      for (i in 1:ntree) obj[[i]]$tip.label <- ape::checkLabel(obj[[i]]$tip.label)
    }
    else {
      attr(obj, "TipLabel") <- ape::checkLabel(attr(obj, "TipLabel"))
      obj <- ape::.uncompressTipLabel(obj)
    }
  }
  
  for(i in 1:ntree) {
    treedata[[i]]@phylo <- obj[[i]]
    root.tag <- if (treeio::is.rooted(treedata[[i]])) "= [&R] " else "= [&U] "
    
    cat("\tTREE *", paste0(tree.name,"_",i-1), root.tag, file = file, append = TRUE)
    cat(treeio::write.beast.newick(treedata[[i]], file = ""), "\n", sep = "",
        file = file, append = TRUE)
  }
  
  cat("END;\n", file = file, append = TRUE)
}
