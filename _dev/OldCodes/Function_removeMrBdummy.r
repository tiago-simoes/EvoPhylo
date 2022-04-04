# library(ape)
# library(treeio)

#### FUNCTION 1: REMOVE "DUMMY" FROM Mr. Bayes MRC/MCT TREES  (modified from convert.MCC.tree) ######

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

######## EXAMPLE RUNS

convert.MrB.trees(treefile = "_dev/Examples/NewFunctions/TetrapodMRC(dummy).t",
                  outfile = "_dev/Examples/NewFunctions/TetrapodMRC(NOdummy).t")

