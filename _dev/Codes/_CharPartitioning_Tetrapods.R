library(ape)
library(cluster)
library(dplyr) 
library(ISLR) 
library(Rtsne)
library(factoextra)
library(naniar)
library(psych)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(StatMatch)


####### IMPORT Raw Characters
Dataraw <- read.csv("Tetra_Characters.csv")

#Transpose rows and columns 

DataDF <- as.data.frame(t(Dataraw))
Data1 <- DataDF %>%  mutate_all(as.ordered) #If multistates treated as ordered variables for char distances
Data2 <- DataDF %>% mutate_all(as.numeric) #If multistates treated as unordered variables for char distances

############# DISTANCE MATRIX  ####################


Data_M <- as.matrix(Data1)
Dmatrix <- gower.dist(Data_M, KR.corr=TRUE)
Dmatrix <- as.matrix(Dmatrix)
summary(Dmatrix)
write.csv(Dmatrix, file="Dmatrix.csv")

############ Cluster analysis using partitioning around medoids (PAM)
### Use Dmatrix above or import

Dmatrix <- read.csv("Dmatrix.csv", header=TRUE, row.names = 1)
Dmatrix <- as.matrix(Dmatrix)


############ k selection using silhouette index
sil_width <- c(NA)

for(i in 2:10){
  
  pam_fit <- pam(Dmatrix,
                 diss = TRUE,
                 k = i)
  
  sil_width[i] <- pam_fit$silinfo$avg.width
  
}

sil_widthDF<-as.data.frame(sil_width)
write.table(sil_widthDF, file="ClusterSizes_K.csv")

# Plot sihouette width (higher is better)

plot(1:10, sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(1:10, sil_width)

##############Clustering with PAM under the chosen K value

pam_fit <- pam(Dmatrix, diss = TRUE, k = 3)

summary(pam_fit)


#ClusterStats
ClusterStats<-pam_fit$clusinfo
write.table(ClusterStats, file="ClusterStats.csv")

#Variable-clusterings
Clusters <- as.data.frame(pam_fit$clustering)
df.Clusters = data.frame(character.name = row.names(Clusters), character.number = seq(1,nrow(Clusters)), cluster = Clusters$`pam_fit$clustering`)
write.csv(df.Clusters, file="ClusterTable.csv")

# Group rows by cluster name 
s = df.Clusters %>% group_by(cluster) %>% summarise(vector=paste(character.number, collapse=" "))

# Write clusters as character partitions in NEXUS format

file<-file("Clusters_Nexus.txt")
cat("#NEXUS\n", file = file)
cat(paste("charset morph_p", s$cluster," = ",s$vector, ";", sep=""), sep = "\n", file = "Clusters_Nexus.txt", append = TRUE)
cat("partition all","=" ,nrow(s1), ":", paste("morph_p", s1$cluster, sep="", collapse = ","),  ";\n", file = "ClustersHom_3p_Nexus.txt", append = TRUE)
cat("set partition=all;\n", file = "Clusters_Nexus.txt", append = TRUE)
readLines(file)
close(file)



########## tSNE graphic cluster analysis
#
#

tsne <- Rtsne(Dmatrix, is_distance = TRUE, theta=0.0, dims=2)

tsne_data <- tsne$Y %>% data.frame() %>% setNames(c("tSNE_Dim1", "tSNE_Dim2")) %>%
  mutate(cluster = factor(pam_fit$clustering),
         Character_number = df.Clusters$character.number,
         Character_name = df.Clusters$character.name)

ggplot(tsne_data, aes(x = tSNE_Dim1, y = tSNE_Dim2)) +
  geom_point()+
  geom_text_repel(aes(label=Character_number, color = cluster))