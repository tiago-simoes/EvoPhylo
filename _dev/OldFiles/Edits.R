######Additions and changes

library(tiagospackage)


#-1) Log File Data Import
# #Get FBD parameter estimates from collection of log files (.p)
# files <- list.files(pattern='\\.p$')
# AllRuns <- do.call(rbind, lapply(files, function(x)
#   read.table(x, skip=1, header = TRUE))) #I am using skip=1 to remove the first line of each log file that contains an unwanted ID string/Beyond that, I also wanted to remove the first 25% of all rows in each file
#
# #Filtering down to 10k rows by skipping every 14 rows to further reduce file size
# AllRuns<- AllRuns[seq(1, nrow(AllRuns), 14),] #this should create a final downsampled log file equivalent to Log(4runs)10k_3p or Log(4runs)10k_1p
#
# #Export reduced log file
# write.csv(AllRuns, file="AllRuns.csv")

#-2) Reshape Allruns data from wide to long

# AllRuns_Melted <- FBD_reshape(AllRuns) #Equivalent to AllRuns_COMB_Melted.csv

#BUT, adapt FBD_tests2 and FBD_summary functions. E.g.:


#CURRENT:
# FBD_summary <- function(AllRunsrMelted_MC, file = NULL, digits = 3) {
#   time.bins <- sort(unique(AllRunsrMelted_MC$Time_bin))
#   parameters <- c("net_speciation", "relative_extinction", "relative_fossilization")
#   analyses <- c("Tip", "TipNode", "Combined")
#
#   #NEW:
#   FBD_summary <- function(AllRunsrMelted_MC, file = NULL, digits = 3) {
#     time.bins <- sort(unique(AllRunsrMelted_MC$Time_bin))
#     parameters <- c("net_speciation", "relative_extinction", "relative_fossilization")
#     analyses <- sort(unique(AllRunsrMelted_MC$Analysis)) #For putative new variables created by users or simply remove analyses altogether
#


#-3)Importing phylogenetic data matrices directly from Nexus files (instead of csv files)

#### CharPart

# #CURRENT: Load character data and produce Gower distance matrix
#
# Dmatrix <- get_gower_dist("Tetra_Characters.csv", numeric = FALSE)
#
# #NEW: Incorporate ape function read.nexus.data into this function
# Data_M<-ape::read.nexus.data("DataMatrix_43t_178ch(Cl20)_UPDATED_NoPolys.nex")
# Data_M <- as.data.frame(Data_M)
# Dmatrix <- get_gower_dist(Data_M, numeric = FALSE) #No transposition of the matrix needed in this case (ape seems to be directly importing the character variables as observations as needed)
#


#-4) New argument to choose between 2 or 3 tSNE dimensions to plot
# #CURRENT:
# clusters <- make_clusters(Dmatrix, k = 3, tsne = TRUE, plot = TRUE)
#
# #I would suggest the addition of another argument where users can determine how many t-SNE dimensions they would like to plot (only first 2 or also 3).
# #The reason is that in my experience important insights may come from comparing tSNE_Dim1 to tSNE_Dim2, or tSNE_Dim2 to tSNE_Dim3.
# #tsne_dim = 2 would be the default option.
#
# #NEW:
# clusters <- make_clusters(Dmatrix, k = 3, tsne = TRUE, plot = TRUE, tsne_dim = 2)
#
# #The output in the case of 3 tSNE dimensions would be three graphs next to each other


-#5) NEW FUNCTION (document)
#
# clock_reshape <- function(RatesByClade) {
#   RatesByClade_long <- reshape(RatesByClade, direction = "long",
#                                varying = names(RatesByClade)[startsWith(names(RatesByClade), "rates")],
#                                v.names = c("rates"),
#                                timevar = "clock",
#                                idvar = "nodes",
#                                sep = "_")
#   RatesByClade_long[["clock"]] <- factor(RatesByClade_long[["clock"]])
#   RatesByClade_long
# }



