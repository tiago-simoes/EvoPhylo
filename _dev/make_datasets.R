#Making posterior3p
posterior <- read.table("_dev/Examples/SelectionStrength/3p_TipNodeMC_CombLog_8k.p", header = TRUE)
names(posterior) <- gsub(".all.", "", names(posterior), fixed = T)
names(posterior)[endsWith(names(posterior), ".")] <- substring(names(posterior)[endsWith(names(posterior), ".")], 1, nchar(names(posterior)[endsWith(names(posterior), ".")]) - 1)
burnin <- quantile(posterior$Gen, .25)
posterior <- posterior[posterior$Gen > burnin,]
posterior <- posterior[round(seq(1, nrow(posterior), length.out = 1000)),]
posterior3p <- FBD_reshape(posterior)
save(posterior3p, file = "data/posterior3p.rda")

#Making posterior1p
posterior1p <- import_log(paste0("_dev/new files/Input files/1p_run", 1:4, ".p"),
                        burnin = .25, downsample = 1000)
save(posterior1p, file = "data/posterior1p.rda")

#Making characters
characters <- as.data.frame(ape::read.nexus.data("_dev/new files/Input files/DataMatrix_noPolys.nex"))
save(characters, file = "data/characters.rda")

#Making tree3p
tree3p <- treeio::read.mrbayes("_dev/new files/Input files/Tree_3p(raw).t")
save(tree3p, file = "data/tree3p.rda")

#Making tree1p
tree1p <- treeio::read.mrbayes("_dev/new files/Input files/Tree_1p.t")
save(tree1p, file = "data/tree1p.rda")

#making rate_table_means
#same as get_clockrate_table(tree3p) with clades added
rate_table_means <- read.csv("_dev/new files/Input files/RateTable_Means_Clades.csv")
save(rate_table_means, file = "data/rate_table_means.rda")
