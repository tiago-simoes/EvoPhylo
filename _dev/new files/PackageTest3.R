devtools::install_github("ngreifer/tiagos-package")

library(tiagospackage)
library(ggtree)
library(treeio)
library(ape)


help(, "tiagospackage")


#### CharPart
(Documented)

#Load character data and produce Gower distance matrix
Data_M<-ape::read.nexus.data("DataMatrix_noPolys.nex") ##T: Data matrix cannot include polymorphisms (users will know what it means)
Data_M <- as.data.frame(Data_M)
Dmatrix <- get_gower_dist(Data_M, numeric = FALSE) ##T: USE the UPDATED get_gower_dist function (Without transposing)

#Save D matrix
write.csv(Dmatrix, file="Dmatrix.csv")

#Plot number of cluster against silhouette width
get_sil_widths(Dmatrix, max.k = 10, plot = TRUE)

#Decide on number of clusters based on plot; here, k = 3
#Generate and visualize clusters
clusters <- make_clusters(Dmatrix, k = 3, tsne = TRUE, plot = TRUE)

#Write clusters to Nexus file
cluster_to_nexus(clusters, file = "Clusters_Nexus.txt")


#### ST&P_Rates
(Documented)

#get evol tree with one (1p) or several (e.g., 3p) clock partitions (either in .t or .tre formats)
tree <- treeio::read.beast("Tree_3p(raw).t")
tree@data$`rateTK02Brlens{3}_median`

#Get table of clock rates with given summary stat
RateTable_Means <- get_clockrate_table(tree, summary = "mean")
RateTable_Medians <- get_clockrate_table(tree, summary = "median")

#Export the rate_table
write.csv(RateTable_Means, file="RateTable_Means.csv")
write.csv(RateTable_Medians, file="RateTable_Medians.csv")

#Get and export Node labels
tree_nodes<-ggtree(tree, layout = "rectangular", ladderize=TRUE, right=TRUE, 
                 branch.length="none", size = 0.05)+
  geom_tiplab(size=2, linesize = 0.01, color="black",  offset = 0.2)+
  geom_label(aes(label=node), size=2, color="purple", position = "dodge")
tree_nodes
ggsave("tree_nodes.pdf")


#### Get summary statistic table for each clade by clock partition
#Import customized RateTable_Medians or RateTable_Means with clade names for each node (Example provided)
rate_table_clade <- read.csv("RateTable_Medians_Clades.csv")
#Get summary stats
clockrate_summary(rate_table_clade)
#Export the summary stats for rate_table_clade
write.csv(rate_table_clade, file="Sum_RateTable_Medians.csv")


#Plot distributions of rates by clock and clade
clockrate_dens_plot(rate_table_clade, stack = TRUE, nrow = 3) ##T: Can we add  other color palettes (e.g., viridis?)

#Plot regressions of rates from two clocks
p12 <- clockrate_reg_plot(rate_table_clade, clock_x = 1, clock_y = 2)
p13 <- clockrate_reg_plot(rate_table_clade, clock_x = 1, clock_y = 3)
p23 <- clockrate_reg_plot(rate_table_clade, clock_x = 2, clock_y = 3)
gridExtra::grid.arrange(p12, p13, p23, nrow = 2)


#### SelectionStrength
(Documented)

#Import customized RateTable_Means with clade names for each node (Example provided)
RatesByClade <- read.csv("RateTable_Means_Clades.csv")
RatesByClade <- clock_reshape(RatesByClade)

posterior <- read.table("Log(4runs)10k_3p.p", header = TRUE)
tree <- treeio::read.beast("Tree_3p.tre")

#Get matrix of pairwise t-tests for sign difference between the posterior base of clock rates (background rate) and
# each rate means for every branch in summary tree
f1(RatesByClade, posterior$clockrate.all.)

#Plot tree using one sign threshold
f2(tree, posterior$clockrate.all., clock = 1,
   threshold = c("1 SD"))

#Plot tree using various sign thresholds
f2(tree, posterior$clockrate.all., clock = 1,
   threshold = c("1 SD", "2 SD", "95%"))




#### ST&P_FBD
(Documented)

#Get FBD posterior from log file and Reshape Allruns data from wide to long
AllRuns <- read.table("Log(4runs)10k_1p.p", header = TRUE)
AllRuns_Melted <- FBD_reshape(AllRuns)

#Summarize parameters by time bin and analysis
Sum_FBD <- FBD_summary(AllRuns_Melted)
#Export the summary stats for rate_table_clade
write.csv(Sum_FBD, file="Sum_FBD.csv")

#Plot density of parameter by time bin
FBD_dens_plot(AllRuns_Melted, parameter = "net_speciation", type = "density")
FBD_dens_plot(AllRuns_Melted, parameter = "net_speciation", type = "violin")
FBD_dens_plot(AllRuns_Melted, parameter = "relative_extinction", type = "violin")
FBD_dens_plot(AllRuns_Melted, parameter = "relative_fossilization", type = "violin")

#Plot density of clockrate estimates by analysis type
f3(AllRuns_Melted) #+ coord_cartesian(xlim = c(0, .08)) ##T: Not necessary

#Tests for normality and homoscedasticity for each parameter across time bins
#and analyses
FBD_tests1(AllRuns_Melted) 

#Visualize deviations from normality and similarity of variances
FBD_normality_plot(AllRuns_Melted) 

#Test differences in location for each parameter between time bins 
#for each analysis type
FBD_tests2(AllRuns_Melted, contrast = "Time_bin") 

#Test differences in location for each parameter between analysis types
#for each time bin
#FBD_tests2(AllRuns_Melted, contrast = "analysis") ##T: Remove From vignette
