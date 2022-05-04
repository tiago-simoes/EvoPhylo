# setwd ("E:/Git/EvoPhylo/")

#Making posterior3p
posterior3p <- combine_log("_dev/extdata/MultiClockTree/LogFiles3p/",
                        burnin = .25, downsample = 10000)
save(posterior3p, file = "data/posterior3p.rda")

#Making posterior1p
posterior1p <- combine_log("_dev/extdata/SingleClockTree/LogFiles1p/",
                        burnin = .25, downsample = 10000)
save(posterior1p, file = "data/posterior1p.rda")

#Making characters
characters <- as.data.frame(ape::read.nexus.data("_dev/extdata/DataMatrix.nex"))
save(characters, file = "data/characters.rda")

#Making tree3p
tree3p <- treeio::read.mrbayes("_dev/extdata/MultiClockTree/Tree3p.t")
save(tree3p, file = "data/tree3p.rda")

#Making tree1p
tree1p <- treeio::read.mrbayes("_dev/extdata/SingleClockTree/Tree1p.t")
save(tree1p, file = "data/tree1p.rda")

#making rate_table_clades_means3
#same as get_clockrate_table(tree3p) with clades added
rate_table_clades_means3 <- read.csv("_dev/extdata/MultiClockTree/RateTable_Means_Clades.csv")
save(rate_table_clades_means3, file = "data/rate_table_clades_means3.rda")

#making rate_table_clades_means1
#same as get_clockrate_table(tree1p) with clades added
rate_table_clades_means1 <- read.csv("_dev/extdata/SingleClockTree/RateTable_Means_Clades.csv")
save(rate_table_clades_means1, file = "data/rate_table_clades_means1.rda")
