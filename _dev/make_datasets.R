setwd ("E:/Git/EvoPhylo/")

######### Morphologicl data matrix example file

characters <- as.matrix(characters)
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

#making RateTable_Means_3p_Clades
#same as get_clockrate_table(tree3p) with clades added
RateTable_Means_3p_Clades <- read.csv("_dev/extdata/MrBayes/MultiClockTree/RateTable_Means_Clades.csv")
save(RateTable_Means_3p_Clades, file = "data/RateTable_Means_3p_Clades.rda")

#making RateTable_Means_1p_Clades
#same as get_clockrate_table(tree1p) with clades added
RateTable_Means_1p_Clades <- read.csv("_dev/extdata/MrBayes/SingleClockTree/RateTable_Means_Clades.csv")
save(RateTable_Means_1p_Clades, file = "data/RateTable_Means_1p_Clades.rda")

######### BEAST2 example files

#multiple post trees file
post_trees_path = system.file("extdata", "ex_offset.trees", package = "EvoPhylo")
post_trees <- treeio::read.beast(post_trees_path)
save(post_trees, file = "data/post_trees.rda")

#multiple clocks tree files
tree_file_clock1 = system.file("extdata", "ex_multiclock.clock1.MCC.tre", package = "EvoPhylo")
tree_clock1 = treeio::read.beast(tree_file_clock1)
save(tree_clock1, file = "data/tree_clock1.rda")

tree_file_clock2 = system.file("extdata", "ex_multiclock.clock2.MCC.tre", package = "EvoPhylo")
tree_clock2 = treeio::read.beast(tree_file_clock2)
save(tree_clock2, file = "data/tree_clock2.rda")


#Handling BEAST2 trees with offsets
posterior_trees_offset <- treeio::read.beast("_dev/extdata/BEAST2/OffsetTrees/ex_offset.trees")
save(posterior_trees_offset, file = "data/posterior_trees_offset.rda")

posterior_log_offset <- read.table("_dev/extdata/BEAST2/OffsetTrees/ex_offset.log", header = TRUE)
save(posterior_log_offset, file = "data/posterior_log_offset.rda")

mcc_dummy <- treeio::read.beast("_dev/extdata/BEAST2/OffsetTrees/ex_offset.MCC.tre")
save(mcc_dummy, file = "data/mcc_dummy.rda")


