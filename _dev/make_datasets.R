# setwd ("E:/Git/EvoPhylo/")

######### Morphologicl data matrix example file
characters <- as.data.frame(ape::read.nexus.data("_dev/extdata/MrBayes/DataMatrix.nex"))
save(characters, file = "data/characters.rda")

######### MrBayes example files

#Making posterior3p
posterior3p <- combine_log("_dev/extdata/MrBayes/MrBayes/MultiClockTree/LogFiles3p/",
                        burnin = .25, downsample = 10000)
save(posterior3p, file = "data/posterior3p.rda")

#Making posterior1p
posterior1p <- combine_log("_dev/extdata/MrBayes/SingleClockTree/LogFiles1p/",
                        burnin = .25, downsample = 10000)
save(posterior1p, file = "data/posterior1p.rda")

#Making tree3p
tree3p <- treeio::read.mrbayes("_dev/extdata/MrBayes/MultiClockTree/Tree3p.t")
save(tree3p, file = "data/tree3p.rda")

#Making tree1p
tree1p <- treeio::read.mrbayes("_dev/extdata/MrBayes/SingleClockTree/Tree1p.t")
save(tree1p, file = "data/tree1p.rda")

#making rate_table_clades_means3
#same as get_clockrate_table(tree3p) with clades added
rate_table_clades_means3 <- read.csv("_dev/extdata/MrBayes/MultiClockTree/RateTable_Means_Clades.csv")
save(rate_table_clades_means3, file = "data/rate_table_clades_means3.rda")

#making rate_table_clades_means1
#same as get_clockrate_table(tree1p) with clades added
rate_table_clades_means1 <- read.csv("_dev/extdata/MrBayes/SingleClockTree/RateTable_Means_Clades.csv")
save(rate_table_clades_means1, file = "data/rate_table_clades_means1.rda")

######### BEAST2 example files

#Handling BEAST2 trees with offsets
posterior_trees_offset <- treeio::read.beast("_dev/extdata/BEAST2/OffsetTrees/ex_offset.trees")
save(posterior_trees_offset, file = "data/posterior_trees_offset.rda")

posterior_log_offset <- read.table("_dev/extdata/BEAST2/OffsetTrees/ex_offset.log", header = TRUE)
save(posterior_log_offset, file = "data/posterior_log_offset.rda")

mcc_dummy <- treeio::read.beast("_dev/extdata/BEAST2/OffsetTrees/ex_offset.MCC.tre")
save(mcc_dummy, file = "data/mcc_dummy.rda")


