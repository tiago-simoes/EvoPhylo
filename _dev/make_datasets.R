setwd ("E:/Git/EvoPhylo/")

#Making posterior3p
posterior3p <- combine_log("examples/MultiClockTree/LogFiles3p/",
                        burnin = .25, downsample = 10000)
save(posterior3p, file = "data/posterior3p.rda")

#Making posterior3p_long
posterior3p_long <- FBD_reshape(posterior3p)
save(posterior3p_long, file = "data/posterior3p_long.rda")

#Making posterior1p
posterior1p <- combine_log("examples/SingleClockTree/LogFiles1p/",
                        burnin = .25, downsample = 10000)
save(posterior1p, file = "data/posterior1p.rda")

#Making posterior1p_long
posterior1p_long <- FBD_reshape(posterior1p)
save(posterior1p_long, file = "data/posterior1p_long.rda")

#Making characters
characters <- as.data.frame(ape::read.nexus.data("examples/DataMatrix.nex"))
save(characters, file = "data/characters.rda")

#Making tree3p
tree3p <- treeio::read.mrbayes("examples/MultiClockTree/Tree3p.t")
save(tree3p, file = "data/tree3p.rda")

#Making tree1p
tree1p <- treeio::read.mrbayes("examples/SingleClockTree/Tree1p.t")
save(tree1p, file = "data/tree1p.rda")

#making rate_table_clades_means
#same as get_clockrate_table(tree3p) with clades added
rate_table_clades_means <- read.csv("examples/MultiClockTree/RateTable_Means_Clades.csv")
save(rate_table_clades_means, file = "data/rate_table_clades_means.rda")

#making rate_table_clades_medians
#same as get_clockrate_table(tree3p) with clades added
rate_table_clades_medians <- read.csv("examples/MultiClockTree/RateTable_Medians_Clades.csv")
save(rate_table_clades_medians, file = "data/rate_table_clades_medians")
