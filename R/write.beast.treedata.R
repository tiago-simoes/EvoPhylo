#' Export multiple treedata objects (S4 class tree files) to BEAST NEXUS file
#'
#' This function was adopted and modified from [treeio::write.beast()] to export a
#' list of trees instead of a single tree.
#'
#' @param treedata An S4 class object of type `treedata` containing
#' multiple trees; e.g. a Bayesian clock tree distribution imported using
#' [treeio::read.beast()] or [treeio::read.mrbayes()].
#' @param file Output file. If `file = ""`, prints the output content on
#' screen.
#' @param translate Whether to translate taxa labels.
#' @param tree.name Name of the trees, default `"STATE"`.
#'
#' @returns
#' Writes object type `treedata` containing multiple trees to a
#' file or file content on screen.
#'
#' @examples
#' #Load file with multiple trees
#' \dontrun{
#' trees_file = system.file("extdata", "ex_offset.trees", package = "EvoPhylo")
#' posterior_trees_offset = treeio::read.beast(trees_file)
#'
#' #Write multiple trees to screen
#' write.beast.treedata(posterior_trees_offset)
#' }

#' @export
write.beast.treedata <- function(treedata, file = "", translate = TRUE, tree.name = "STATE") {

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
    for (i in seq_len(ntree)) {
      obj[[i]]$tip.label <- as.character(seq_len(N))
    }
  }
  else if (is.null(attr(obj, "TipLabel"))) {
    for (i in seq_len(ntree)) {
      obj[[i]]$tip.label <- ape::checkLabel(obj[[i]]$tip.label)
    }
  }
  else {
    attr(obj, "TipLabel") <- ape::checkLabel(attr(obj, "TipLabel"))
    obj <- ape::.uncompressTipLabel(obj)
  }

  for (i in seq_len(ntree)) {
    treedata[[i]]@phylo <- obj[[i]]
    root.tag <- if (treeio::is.rooted(treedata[[i]])) "= [&R] " else "= [&U] "

    cat("\tTREE *", paste0(tree.name, "_", i - 1), root.tag, file = file, append = TRUE)
    cat(treeio::write.beast.newick(treedata[[i]], file = ""), "\n", sep = "",
        file = file, append = TRUE)
  }

  cat("END;\n", file = file, append = TRUE)
}
