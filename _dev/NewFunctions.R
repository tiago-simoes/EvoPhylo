#New functions developed by Joëlle Barido-Sottani (see references therein) or by
#us to import tree files from BEAST2 and Mr. Bayes with metadata embedded (e.g.,
#evolutionary rates), and to remove "dummy taxa" introduced by the "treeWoffset"
#operator from BEAST2.

#NOTE: Needs examples to test, needs packages from Bioconductor

## Function 1: Add dummy extant to posterior trees from BEAST2 keeping tree metadata (modified by J. Barido-Sottani)
# Barido-Sottani, J., Aguirre-Fern?ndez, G., Hopkins, M. J., Stadler, T., Warnock, R. 2019. Ignoring stratigraphic age uncertainty leads to erroneous estimates of species divergence times under the fossilized birth-death process. Proc. R. Soc. Lond., Ser. B: Biol. Sci. 286, 20190685.

transform.woffset.root.metadata <- function(treesfile, logfile, outfile) {
  trees <- treeio::read.beast(treesfile)
  log <- read.table(logfile, header = TRUE)
  log <- log$offset
  # print(median(log))

  for(i in seq_along(trees)) trees[[i]]@phylo$offset <- log[i]

  #Add dummy extant taxa (and update ages)
  presenttrees <- lapply(trees, function(tmp) {
    t <- tmp@phylo # get the phylo class of the S4 class tmp tree and save as t (S3 class)
    ntips <- length(t$tip.label)
    totn <- ntips + t$Nnode

    times <- ape::node.depth.edgelength(t)
    root_time <- max(times) + t$offset

    root <- ntips+2
    t$edge[t$edge > ntips] <- t$edge[t$edge > ntips] + 1
    t$edge[t$edge >= root] <- t$edge[t$edge >= root] + 1

    t$edge <- rbind(c(root, root+1), c(root, ntips+1), t$edge)
    t$edge.length <- c(0.1, root_time, t$edge.length)
    t$tip.label <- c(t$tip.label, "dummy")
    t$Nnode <- t$Nnode + 1

    tmp@data$node <- as.numeric(tmp@data$node)
    tmp@data$node[tmp@data$node > ntips] <- tmp@data$node[tmp@data$node > ntips] + 1
    tmp@data$node[tmp@data$node >= root] <- tmp@data$node[tmp@data$node >= root] + 1
    tmp@data$node <- as.character(tmp@data$node)
    tmp@data <- rbind(tmp@data, list(rate = 0, node = as.character(ntips + 1)))

    tmp@phylo <- t # save the modified S3 class tree file (t) as the phylo component of the S4 class tmp t
    return(tmp) # save that as presenttrees
  })

  write.beast.treedatas(presenttrees, file = outfile) # save the final tree with metadata using function described below
  invisible(presenttrees)
}

##NEW Function- Write presenttrees from above as trees with metadata

write.beast.treedatas <- function(treedatas, file = "", translate = TRUE, tree.name = "STATE"){
  cat("#NEXUS\n", file = file)
  cat(paste("[R-package treeio, ", date(), "]\n\n", sep = ""),
      file = file, append = TRUE)
  N <- treeio::Ntip(treedatas[[1]])

  obj <- lapply(treedatas, treeio::as.phylo)
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

  #allow or not translate taxon names data block
  if (translate) {
    cat("\tTRANSLATE\n", file = file, append = TRUE)
    obj <- compressTipLabel(obj)
    X <- paste("\t\t", seq_len(N), "\t", attr(obj, "TipLabel"),
               ",", sep = "")
    X[length(X)] <- gsub(",", "", X[length(X)])
    cat(X, file = file, append = TRUE, sep = "\n")
    cat("\t;\n", file = file, append = TRUE)
    class(obj) <- NULL
    for (i in seq_len(ntree)) obj[[i]]$tip.label <- as.character(seq_len(N))
  }
  else {
    if (is.null(attr(obj, "TipLabel"))) {
      for (i in seq_len(ntree)) obj[[i]]$tip.label <- treeio::checkLabel(obj[[i]]$tip.label)
    }
    else {
      attr(obj, "TipLabel") <- treeio::checkLabel(attr(obj, "TipLabel"))
      obj <- uncompressTipLabel(obj)
    }
  }


  #replace @phylo from  presenttrees with new @phylo (obj)
  for (i in seq_len(ntree)) {
    treedatas[[i]]@phylo <- obj[[i]]
    root.tag <- if (treeio::is.rooted(treedatas[[i]])) "= [&R] " else "= [&U] "

    cat("\tTREE *", paste0(tree.name,"_",i-1), root.tag, file = file, append = TRUE)
    cat(treeio::write_beast_newick(treedatas[[i]], file = ""), "\n", sep = "",
        file = file, append = TRUE)
  }

  cat("END;\n", file = file, append = TRUE)
}

uncompressTipLabel <- function(x) {
  #Taken from ape:::.uncompressTipLabel
  Lab <- attr(x, "TipLabel")
  if (is.null(Lab)) return(x)
  class(x) <- NULL
  for (i in seq_along(x)) x[[i]]$tip.label <- Lab
  class(x) <- "multiPhylo"
  attr(x, "TipLabel") <- NULL
  x
}

compressTipLabel <- function (x, ref = NULL) {
  #Taken from ape:::.compressTipLabel
  if (!is.null(attr(x, "TipLabel")))
    return(x)
  if (is.null(ref))
    ref <- x[[1]]$tip.label
  n <- length(ref)
  if (length(unique(ref)) != n)
    stop("Some tip labels are duplicated in tree no. 1.", call. = FALSE)
  relabel <- function(y) {
    label <- y$tip.label
    if (!identical(label, ref)) {
      if (length(label) != length(ref))
        stop("One tree has a different number of tips.", call. = FALSE)
      ilab <- match(label, ref)
      if (anyNA(ilab))
        stop("One tree has different tip labels.", call. = FALSE)
      ie <- match(seq_len(n), y$edge[, 2])
      y$edge[ie, 2] <- ilab
    }
    y$tip.label <- NULL
    y
  }
  x <- unclass(x)
  x <- lapply(x, relabel)
  attr(x, "TipLabel") <- ref
  class(x) <- "multiPhylo"
  x
}


#### FUNCTION 2: REMOVE "DUMMY" extant FROM BEAST MCC TREE ######
#Barido-Sottani, J., Aguirre-Fern?ndez, G., Hopkins, M. J., Stadler, T., Warnock, R. 2019. Ignoring stratigraphic age uncertainty leads to erroneous estimates of species divergence times under the fossilized birth-death process. Proc. R. Soc. Lond., Ser. B: Biol. Sci. 286, 20190685.

convert.MCC.tree = function(treefile, outfile) {
  tmp <- treeio::read.beast(treefile)

  tip <- which(tmp@phylo$tip.label == "dummy")
  tipn <- which(tmp@data$node == as.character(tip))
  node <- tmp@phylo$edge[tmp@phylo$edge[,2] == tip,1]
  noden <- which(tmp@data$node == as.character(node))

  # offset <- min(tmp@data$height_median[as.numeric(tmp@data$node) <= length(tmp@phylo$tip.label) &
  #                                             as.numeric(tmp@data$node) != tip])
  # print(offset)

  tmp@phylo <- treeio::drop.tip(tmp@phylo, tip)
  tmp@data <- tmp@data[-c(tipn, noden),]
  tmp@data$node <- as.numeric(tmp@data$node)
  tmp@data$node[tmp@data$node > tip] <- tmp@data$node[tmp@data$node > tip] - 1
  tmp@data$node[tmp@data$node > node] <- tmp@data$node[tmp@data$node > node] - 1
  tmp@data$node <- as.character(tmp@data$node)

  tmp@data$height_0.95_HPD[] <- lapply(tmp@data$height_0.95_HPD, function(x) x-offset)
  treeio::write.beast(tmp,outfile)
  system(paste0("gsed -i 's%* UNTITLED%TREE1%' ", outfile))
}



#### FUNCTION 3: REMOVE "DUMMY" FROM Mr. Bayes MRC/MCT TREES  (modified by T. Sim?es from convert.MCC.tree) ######

convert.MrB.trees = function(treefile, outfile) {

  tmp <- treeio::read.mrbayes(treefile)

  tmp@data$`prob+-sd` <- as.factor(tmp@data$`prob+-sd`)

  char_col <- vapply(tmp@data, is.character, logical(1L))
  tmp@data[char_col] <- lapply(tmp@data[char_col], as.numeric)

  tip <- which(tmp@phylo$tip.label == "Dummyextant")
  tipn <- which(tmp@data$node == as.character(tip))
  node <- tmp@phylo$edge[tmp@phylo$edge[,2] == tip,1]
  noden <- which(tmp@data$node == as.character(node))

  # offset <- min(tmp@data$age_median[which(as.numeric(tmp@data$node) <= length(tmp@phylo$tip.label) &
  #                                             as.numeric(tmp@data$node) != tip)])

  tmp@phylo <- treeio::drop.tip(tmp@phylo, tip)
  tmp@data <- tmp@data[-c(tipn,noden),]
  tmp@data$node <- as.numeric(tmp@data$node)
  tmp@data$node[tmp@data$node > tip] <- tmp@data$node[tmp@data$node > tip] - 1
  tmp@data$node[tmp@data$node > node] <- tmp@data$node[tmp@data$node > node] - 1
  tmp@data$node <- as.character(tmp@data$node)

  treeio::write.beast(tmp,outfile)
  # system(paste0("gsed -i 's%* UNTITLED%TREE1%' ", outfile))
}



