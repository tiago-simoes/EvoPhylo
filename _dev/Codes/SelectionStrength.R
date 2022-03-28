library(ape)
library(phylotate)
library(rlist)
library(ggtree)
library(phytools)
library(strap)
library(devtools)
library(ggplot2)
library(deeptime)
library(treeio)
library(viridisLite)
library(viridis)
library(cowplot)
library(plotly)
library(ggrepel)
library(plyr)
library(lattice)
library(reshape2)
library(FSA)
library(car)
library(Rmisc)
library(dplyr)
library(tidyr)
library(scales)
library(rstatix)

#Deeptime
library(devtools)
install_github("willgearty/deeptime")

#BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")

#Treeio
BiocManager::install(c("treeio"))

#Latest GGTREE
remotes::install_github("YuLab-SMU/ggtree")
##################### TIP+Node MultiCons Plots




######################### SECTION 1: STATISTICAL TESTS#########################


# Import melted from RWTY script
AllRunsrMelted_MC <- read.csv("_dev/Examples/SelectionStrength/AllRuns_COMB_Melted.csv")
names(AllRunsrMelted_MC)

#Downsample posterior sample to 1k (necessary for some stats tests)
Down_AllRunsrMelted_MC<- AllRunsrMelted_MC[round(seq(1, nrow(AllRunsrMelted_MC), length.out = 1000)),]

write.csv(Down_AllRunsrMelted_MC, file= "Down1k_1p_Tip_Par_Melted.csv")

#Order groups by Time_bin

AllRunsrMelted_MC$Time_bin = factor(AllRunsrMelted_MC$Time_bin,
                                    levels=unique(AllRunsrMelted_MC$Time_bin))

summary(AllRunsrMelted_MC)



####### Density plot of clockrate estimates from posterior from TIP vs TIP+NODE analyses
#clockrate tip vs tip+node

Pl1a<- ggplot(AllRunsrMelted_MC, aes(clockrate, fill = analysis)) +
  geom_density(position="dodge",  alpha = 0.6)+
  scale_x_continuous(limits=c(0,0.075))+
  theme(legend.position="top")
Pl1a



########## Tests for Relative rates diff from clockrate (Tip+Node analysis) #################

#Import mean relative rates for each clock partition for all clades from summary tree
RatesByClade <- read.csv("RateByClockrMelted_means.csv", header = TRUE)
names(RatesByClade)

#Import all clockrate values from full posterior output
posterior <- read.table("_dev/Examples/SelectionStrength/3p_TipNodeMC_CombLog_8k.p", header = TRUE)
names(posterior)

#Summary stats of clockrate parameter from the posterior (absolute rate)

Sum_clockrate<-Summarize(posterior$clockrate.all.,
                         data=posterior,
                         digits=9)
Sum_clockrate<- as.data.frame(t(Sum_clockrate))
Sum_clockrate

write.csv(Sum_clockrate, file="Sum_posterior_clockrate.csv")


#Get absolute rate per branch (branch relative means x tree background clockrate mean)
RatesByClade_Abs <- RatesByClade
RatesByClade_Abs$rate <-RatesByClade_Abs$rate*Sum_clockrate$mean


#Summary stats of MCT absolute branch rates

Sum_AbsBrRate<-Summarize(RatesByClade_Abs$rate,
                         data=RatesByClade_Abs,
                         digits=9)
Sum_AbsBrRate<- as.data.frame(t(Sum_AbsBrRate))
Sum_AbsBrRate

write.csv(Sum_AbsBrRate, file="Sum_MCT_AbsolRate.csv")



########### t-tests: use Absolute Rate Means
###### NOT RECOMMENDED - TEST TOO WEAK - very small data dispersion (see clockrate TipNode plot)
###### Use instead at least +-1 SD (see below)

# One-sample t-test for all mean post estimates of absolute branch rates

#H0: Is each absolute branch rate (mu) signif equal to the mean of the background clockrate?
#Test for first branch (rate = mu = 0.0118122079)

posterior$clockrate.all. #All background clockrates (mean = 0.0125)

X2 <- t.test(posterior$clockrate.all., mu = 0.0118122079, conf.int = 0.99 )
X2
X2$conf.int[2]

#Reject H0 (they are sign diff). 95 percent confidence interval: 0.01244794 0.01258599


# Now Loop for all MCT branches mean rate values
test <- vector(mode = "list", length = 7)
out<-data.frame(matrix(as.numeric(), nrow=length(RatesByClade_Abs$rate), ncol=2))


for(i in 1:length(RatesByClade_Abs$rate)){
  test[[i]] <- t.test(posterior$clockrate.all., mu = RatesByClade_Abs$rate[i], conf.int = 0.95)
  out[i, 1] <- as.numeric(test[[i]]$null.value)
  out[i, 2] <- as.numeric(test[[i]]$p.value)
}

out <- cbind(RatesByClade_Abs$clade,
             RatesByClade_Abs$nodes,
             RatesByClade_Abs$clock,
             RatesByClade$rate,
             RatesByClade_Abs$rate,
             out)
names(out) <-c("clade", "nodes", "clock", "Relative rate","Absolute rate (mean)", "null", "p.value")

write.csv(out, file="RateSign_T.test(Tip+Node MCT).csv")


#
#



######################### SECTION 2: PLOTTING RESULTS ON  MrBayes tree MCT  #########################
# For plotting trees without dummy extant, use the following only (not remove.MrB.extant) - SINGLE & MULTI CLOCK


MCT_TRE2<- read.mrbayes("_dev/Examples/SelectionStrength/3p_TK02_FT_SFBD4l_TipNode_MultiTopCons_AllCom.t")
names(MCT_TRE2@data)

#Drop extant "dummy" tip
MCT_TRE2 <- drop.tip(MCT_TRE2, "Dummyextant")

#MULTI clock tree (3 partitions file)

MCT_TRE2@data$`prob+-sd` <- as.factor(MCT_TRE2@data$`prob+-sd`)
MCT_TRE2@data <- MCT_TRE2@data %>% mutate_if(is.character,as.numeric)

offset = min(MCT_TRE2@data$age_median)
offset

center <- MCT_TRE2@data$age_median

MCT_TRE2@data$age_0.95HPD


# #POsterior
#
# Pl1b<-ggtree(MCT_TRE2, position = position_nudge(x = -offset),
#              layout = "rectangular", ladderize=TRUE, right=TRUE, aes(color=prob))+
#   geom_tiplab(size=2.1, linesize = 0.01, fontface = "italic", color="black", offset = -offset)+
#   geom_tiplab(aes(x=branch, label=round(prob, digits = 2)), vjust=-.5, size=1.5, offset = -offset)+
#   scale_colour_gradient(low = "blue", high = "red",
#                         na.value = "grey50",
#                         guide = "colourbar", aesthetics = "colour")+
#   geom_nodelab(aes(x=branch, label=round(prob, digits = 2)), vjust=-.5, size=2, nudge_x = -offset)+
#   coord_geo(xlim=c(-450,-250), ylim=c(0,Ntip(MCT_TRE2)+2), expand=FALSE,
#             dat = list("epochs", "periods"), abbrv = list(TRUE, FALSE),
#             skip = c("Pliocene", "Pleistocene", "Holocene","Quaternary"),
#             pos = list("bottom", "bottom"),alpha = 1, height = unit(1, "line"),
#             rot = 0, size = list(2.5,3), neg = TRUE) +
#   scale_x_continuous(breaks=seq(-450,-250,10), labels=abs(seq(-450,-250,10)))+
#   theme_tree2()+
#   theme(legend.position=c(.8, .2))
# Pl1b
# Pl1b<- revts(Pl1b)
# Pl1b
#
#
# #Ages
# MCT_raw2<-read.mrbayes("3p_TK02_FT_SFBD4l_TipNode_MultiTopCons_AllCom.t")
#
# MCT_raw2@data$`prob+-sd` <- as.factor(MCT_raw2@data$`prob+-sd`)
# MCT_raw2@data <- MCT_raw2@data %>% mutate_if(is.character,as.numeric)
#
#
# Pl2b<-ggtree(MCT_raw2, layout = "rectangular", ladderize=TRUE, right=TRUE, color="black",
#              size=I(MCT_raw2@data$prob))%>% flip(65, 75)+
#   #geom_tiplab(aes(x=branch, label=round(prob, digits = 2)), vjust=-.5, size=1.5, offset = -offset)+
#   geom_range(range='age_0.95HPD', color='purple', alpha=.4, size=1.5)+
#   geom_nodelab(aes(x=branch, label=round(age_median, digits = 1)), vjust=-.2, hjust=-1.2, size=2)+
#   geom_tiplab(size=3, linesize = 0.01, fontface = "italic", color="black")+
#   coord_geo(xlim=c(-450,-260), ylim=c(0,Ntip(MCT_raw2)+2), expand=FALSE,
#             dat = list("epochs", "periods"), abbrv = list(TRUE, FALSE),
#             skip = c("Pliocene", "Pleistocene", "Holocene","Quaternary"),
#             pos = list("bottom", "bottom"),alpha = 1, height = unit(1, "line"),
#             rot = 0, size = list(2.5,3), neg = TRUE) +
#   scale_x_continuous(breaks=seq(-450,-260,10), labels=abs(seq(-450,-260,10)))+
#   theme_tree2()+
#   theme(legend.position=c(.8, .2))
# Pl2b
# Pl2b <- revts(Pl2b)
# Pl2b



############## VISUALIZING SIGNIFICANT RATE CHANGES ON TREES (Rate Means of final tip+node calibrated tree)
#Use upper (clockrate + iSD) and lower (clockrate - iSD) bounds

#Get mean absolute background clockrate +- 95%CI (from t-test): sign branch rate dev from background
Upper95CI <- X2$conf.int[2]
Lower95CI <- X2$conf.int[1]

#Get mean absolute background clockrate +- 1SD: sign branch rate dev from background
Upper1SD <- Sum_clockrate$mean+Sum_clockrate$sd
Lower1SD <- Sum_clockrate$mean-Sum_clockrate$sd

#Get mean absolute background clockrate +- 1SD: sign branch rate dev from background
Upper2SD <- Sum_clockrate$mean+(2*Sum_clockrate$sd)
Lower2SD <- Sum_clockrate$mean-(2*Sum_clockrate$sd)


#RELATIVE THRESHOLDS
#Get mean absolute background clockrate +- 95%CI (from t-test): sign branch rate dev from background
RelUpper95CI <- Upper95CI/Sum_clockrate$mean
RelLower95CI <- Lower95CI/Sum_clockrate$mean

#Convert to relative background clockrate +- 1SD (to plot on MCT): sign branch rate dev from background
RelUpper1SD <- Upper1SD/Sum_clockrate$mean
RelLower1SD <- Lower1SD/Sum_clockrate$mean

#Convert to relative background clockrate +- 2SD (to plot on MCT): strong branch rate dev from background
RelUpper2SD <- Upper2SD/Sum_clockrate$mean
RelLower2SD <- Lower2SD/Sum_clockrate$mean


############ Rates1

#Rates1 (Sign shift = Outside 95%CI)

Pl4R1a<-ggtree(MCT_TRE2, layout = "rectangular", ladderize=TRUE, right=TRUE, position = position_nudge(x = -offset),
               aes(color=rateTK02Brlens1_mean), size=2)%>% flip(63, 73)+
  geom_tiplab(size=3, linesize = 0.01, fontface = "italic", color="black", offset = -285)+
  #geom_tiplab(aes(x=branch, label=round(rateTK02Brlens1_mean, digits = 2)), color="black", vjust=-0.5, size=2, offset = -288)+
  scale_colour_steps2("Background Rate Threshold",
                      low ="blue", mid= "gray", high ="red",
                      breaks=c(RelLower95CI, RelUpper95CI),
                      labels =c("Lower 95%CI","Upper 95%CI"),
                      limits = c(0,2),
                      midpoint = 1)+
  #geom_range(range='rateTK02Brlens_0.95HPD', color='red', alpha=.6, size=2) +
  #geom_nodelab(aes(x=branch, label=round(rateTK02Brlens1_mean, digits = 2)), color="black", vjust=-.5, hjust=+.5, size=2, nudge_x = -offset)+
  coord_geo(xlim=c(-450,-260), ylim=c(0,Ntip(MCT_TRE2)+2), expand=FALSE,
            dat = list("epochs", "periods"), abbrv = list(TRUE, FALSE),
            skip = c("Pliocene", "Pleistocene", "Holocene","Quaternary"),
            pos = list("bottom", "bottom"),alpha = 1, height = unit(1, "line"),
            rot = 0, size = list(2.5,3), neg = TRUE) +
  scale_x_continuous(breaks=seq(-450,-260,20), labels=abs(seq(-450,-260,20)))+
  theme_tree2()+
  labs(title = "Significant (weak) Selection")+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position=c(.9, .2),
        legend.title = element_text(size = 10, face = "bold"))
# Pl4R1a
Pl4R1a <- revts(Pl4R1a)
Pl4R1a


#Rates1 (Robust Sign shift => 1SD)
Pl4R1b<-ggtree(MCT_TRE2, layout = "rectangular", ladderize=TRUE, right=TRUE, position = position_nudge(x = -offset),
               aes(color=rateTK02Brlens1_mean), size=2)%>% flip(63, 73)+
  geom_tiplab(size=3, linesize = 0.01, fontface = "italic", color="black", offset = -285)+
  #geom_tiplab(aes(x=branch, label=round(rateTK02Brlens1_mean, digits = 2)), color="black", vjust=-0.5, size=2, offset = -288)+
  scale_colour_steps2("Background Rate Threshold",
                      low ="blue", mid= "gray", high ="red",
                      breaks=c(RelLower1SD, RelUpper1SD),
                      labels =c( "-1 SD", "+1 SD"),
                      limits = c(0,2),
                      midpoint = 1)+
  #geom_range(range='rateTK02Brlens_0.95HPD', color='red', alpha=.6, size=2) +
  #geom_nodelab(aes(x=branch, label=round(rateTK02Brlens1_mean, digits = 2)), color="black", vjust=-.5, hjust=+.5, size=2, nudge_x = -offset)+
  coord_geo(xlim=c(-450,-260), ylim=c(0,Ntip(MCT_TRE2)+2), expand=FALSE,
            dat = list("epochs", "periods"), abbrv = list(TRUE, FALSE),
            skip = c("Pliocene", "Pleistocene", "Holocene","Quaternary"),
            pos = list("bottom", "bottom"),alpha = 1, height = unit(1, "line"),
            rot = 0, size = list(2.5,3), neg = TRUE) +
  scale_x_continuous(breaks=seq(-450,-260,20), labels=abs(seq(-450,-260,20)))+
  theme_tree2()+
  labs(title = "Moderate Selection")+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position=c(.9, .2),
        legend.title = element_text(size = 10, face = "bold"))
# Pl4R1b
Pl4R1b <- revts(Pl4R1b)
Pl4R1b


#Rates1 (Strong shift => 2SD)

Pl4R1c<-ggtree(MCT_TRE2, layout = "rectangular", ladderize=TRUE, right=TRUE, position = position_nudge(x = -offset),
              aes(color=rateTK02Brlens1_mean), size=2)%>% flip(63, 73)+
  geom_tiplab(size=3, linesize = 0.01, fontface = "italic", color="black", offset = -285)+
  #geom_tiplab(aes(x=branch, label=round(rateTK02Brlens1_mean, digits = 2)), color="black", vjust=-0.5, size=2, offset = -288)+
  scale_colour_steps2("Background Rate Threshold",
                      low ="blue", mid= "gray", high ="red",
                      breaks=c(RelLower2SD, RelUpper2SD),
                      labels =c("-2 SD", "+2 SD"),
                      limits = c(0,2),
                      midpoint = 1)+
  #geom_range(range='rateTK02Brlens_0.95HPD', color='red', alpha=.6, size=2) +
  #geom_nodelab(aes(x=branch, label=round(rateTK02Brlens1_mean, digits = 2)), color="black", vjust=-.5, hjust=+.5, size=2, nudge_x = -offset)+
  coord_geo(xlim=c(-450,-260), ylim=c(0,Ntip(MCT_TRE2)+2), expand=FALSE,
            dat = list("epochs", "periods"), abbrv = list(TRUE, FALSE),
            skip = c("Pliocene", "Pleistocene", "Holocene","Quaternary"),
            pos = list("bottom", "bottom"),alpha = 1, height = unit(1, "line"),
            rot = 0, size = list(2.5,3), neg = TRUE) +
  scale_x_continuous(breaks=seq(-450,-260,20), labels=abs(seq(-450,-260,20)))+
  theme_tree2()+
  labs(title = "Strong Selection")+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position=c(.9, .2),
        legend.title = element_text(size = 10, face = "bold"))
Pl4R1c
Pl4R1c <- revts(Pl4R1c)
Pl4R1c


par(mar = c(10, 10, 10, 10), )


A<-plot_grid(Pl4R1a, Pl4R1b, Pl4R1c, ncol=3,
             byrow = FALSE,
             labels = c('A', 'B', 'C'),
             label_size = 12)
A


#Rates1 (Selection gradient)
Pl4R1d<-ggtree(MCT_TRE2, layout = "rectangular", ladderize=TRUE, right=TRUE, position = position_nudge(x = -offset),
               aes(color=rateTK02Brlens1_mean), size=2)%>% flip(63, 73)+
  geom_tiplab(size=3, linesize = 0.01, fontface = "italic", color="black", offset = -285)+
  #geom_tiplab(aes(x=branch, label=round(rateTK02Brlens3_mean, digits = 2)), color="black", vjust=-0.5, size=2, offset = -288)+
  scale_colour_stepsn("Background Rate Threshold",
                      colours = c("blue", "cyan3", "cyan","gray", "orange", "red1" ,"red4"),
                      breaks=c(RelLower2SD,RelLower1SD,RelLower95CI,RelUpper95CI, RelUpper1SD, RelUpper2SD),
                      labels =c("-2 SD", "-1 SD", "Lower 95%CI","Upper 95%CI", "+1 SD", "+2 SD"),
                      limits = c(0,2),
                      guide = "coloursteps",
                      aesthetics = "colour",)+
  #geom_range(range='rateTK02Brlens_0.95HPD', color='red', alpha=.6, size=2) +
  #geom_nodelab(aes(x=branch, label=round(rateTK02Brlens3_mean, digits = 2)), color="black", vjust=-.5, hjust=+.5, size=2, nudge_x = -offset)+
  coord_geo(xlim=c(-450,-260), ylim=c(0,Ntip(MCT_TRE2)+2), expand=FALSE,
            dat = list("epochs", "periods"), abbrv = list(TRUE, FALSE),
            skip = c("Pliocene", "Pleistocene", "Holocene","Quaternary"),
            pos = list("bottom", "bottom"),alpha = 1, height = unit(1, "line"),
            rot = 0, size = list(2.5,3), neg = TRUE) +
  scale_x_continuous(breaks=seq(-450,-260,20), labels=abs(seq(-450,-260,20)))+
  theme_tree2()+
  labs(title = "Selection Strength Gradient")+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position=c(.9, .2),
        legend.title = element_text(size = 10, face = "bold"))
# Pl4R1d
Pl4R1d <- revts(Pl4R1d)
Pl4R1d



############## #Rates2, 3, etc.

#Rates2 (Sign shift = Outside 95%CI)

Pl4R2a<-ggtree(MCT_TRE2, layout = "rectangular", ladderize=TRUE, right=TRUE, position = position_nudge(x = -offset),
               aes(color=rateTK02Brlens3_mean), size=2)%>% flip(63, 73)+
  geom_tiplab(size=3, linesize = 0.01, fontface = "italic", color="black", offset = -285)+
  #geom_tiplab(aes(x=branch, label=round(rateTK02Brlens3_mean, digits = 2)), color="black", vjust=-0.5, size=2, offset = -288)+
  scale_colour_steps2("Background Rate Threshold",
                      low ="blue", mid= "gray", high ="red",
                      breaks=c(RelLower95CI, RelUpper95CI),
                      labels =c("Lower 95%CI","Upper 95%CI"),
                      limits = c(0,2),
                      midpoint = 1)+
  #geom_range(range='rateTK02Brlens_0.95HPD', color='red', alpha=.6, size=2) +
  #geom_nodelab(aes(x=branch, label=round(rateTK02Brlens3_mean, digits = 2)), color="black", vjust=-.5, hjust=+.5, size=2, nudge_x = -offset)+
  coord_geo(xlim=c(-450,-260), ylim=c(0,Ntip(MCT_TRE2)+2), expand=FALSE,
            dat = list("epochs", "periods"), abbrv = list(TRUE, FALSE),
            skip = c("Pliocene", "Pleistocene", "Holocene","Quaternary"),
            pos = list("bottom", "bottom"),alpha = 1, height = unit(1, "line"),
            rot = 0, size = list(2.5,3), neg = TRUE) +
  scale_x_continuous(breaks=seq(-450,-260,20), labels=abs(seq(-450,-260,20)))+
  theme_tree2()+
  labs(title = "Significant (weak) Selection")+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position=c(.9, .2),
        legend.title = element_text(size = 10, face = "bold"))
Pl4R2a
Pl4R2a <- revts(Pl4R2a)
Pl4R2a


#Rates2 (Robust Sign shift => 1SD)
Pl4R2b<-ggtree(MCT_TRE2, layout = "rectangular", ladderize=TRUE, right=TRUE, position = position_nudge(x = -offset),
               aes(color=rateTK02Brlens3_mean), size=2)%>% flip(63, 73)+
  geom_tiplab(size=3, linesize = 0.01, fontface = "italic", color="black", offset = -285)+
  #geom_tiplab(aes(x=branch, label=round(rateTK02Brlens3_mean, digits = 2)), color="black", vjust=-0.5, size=2, offset = -288)+
  scale_colour_steps2("Background Rate Threshold",
                      low ="blue", mid= "gray", high ="red",
                      breaks=c(RelLower1SD, RelUpper1SD),
                      labels =c( "-1 SD", "+1 SD"),
                      limits = c(0,2),
                      midpoint = 1)+
  #geom_range(range='rateTK02Brlens_0.95HPD', color='red', alpha=.6, size=2) +
  #geom_nodelab(aes(x=branch, label=round(rateTK02Brlens3_mean, digits = 2)), color="black", vjust=-.5, hjust=+.5, size=2, nudge_x = -offset)+
  coord_geo(xlim=c(-450,-260), ylim=c(0,Ntip(MCT_TRE2)+2), expand=FALSE,
            dat = list("epochs", "periods"), abbrv = list(TRUE, FALSE),
            skip = c("Pliocene", "Pleistocene", "Holocene","Quaternary"),
            pos = list("bottom", "bottom"),alpha = 1, height = unit(1, "line"),
            rot = 0, size = list(2.5,3), neg = TRUE) +
  scale_x_continuous(breaks=seq(-450,-260,20), labels=abs(seq(-450,-260,20)))+
  theme_tree2()+
  labs(title = "Moderate Selection")+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position=c(.9, .2),
        legend.title = element_text(size = 10, face = "bold"))
Pl4R2b
Pl4R2b <- revts(Pl4R2b)
Pl4R2b


#Rates2 (Strong shift => 2SD)

Pl4R2c<-ggtree(MCT_TRE2, layout = "rectangular", ladderize=TRUE, right=TRUE, position = position_nudge(x = -offset),
               aes(color=rateTK02Brlens3_mean), size=2)%>% flip(63, 73)+
  geom_tiplab(size=3, linesize = 0.01, fontface = "italic", color="black", offset = -285)+
  #geom_tiplab(aes(x=branch, label=round(rateTK02Brlens3_mean, digits = 2)), color="black", vjust=-0.5, size=2, offset = -288)+
  scale_colour_steps2("Background Rate Threshold",
                      low ="blue", mid= "gray", high ="red",
                      breaks=c(RelLower2SD, RelUpper2SD),
                      labels =c("-2 SD", "+2 SD"),
                      limits = c(0,2),
                      midpoint = 1)+
  #geom_range(range='rateTK02Brlens_0.95HPD', color='red', alpha=.6, size=2) +
  #geom_nodelab(aes(x=branch, label=round(rateTK02Brlens3_mean, digits = 2)), color="black", vjust=-.5, hjust=+.5, size=2, nudge_x = -offset)+
  coord_geo(xlim=c(-450,-260), ylim=c(0,Ntip(MCT_TRE2)+2), expand=FALSE,
            dat = list("epochs", "periods"), abbrv = list(TRUE, FALSE),
            skip = c("Pliocene", "Pleistocene", "Holocene","Quaternary"),
            pos = list("bottom", "bottom"),alpha = 1, height = unit(1, "line"),
            rot = 0, size = list(2.5,3), neg = TRUE) +
  scale_x_continuous(breaks=seq(-450,-260,20), labels=abs(seq(-450,-260,20)))+
  theme_tree2()+
  labs(title = "Strong Selection")+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position=c(.9, .2),
        legend.title = element_text(size = 10, face = "bold"))
Pl4R2c
Pl4R2c <- revts(Pl4R2c)
Pl4R2c


par(mar = c(10, 10, 10, 10), )


B<-plot_grid(Pl4R2a, Pl4R2b, Pl4R2c, ncol=3,
             byrow = FALSE,
             labels = c('D', 'E', 'F'),
             label_size = 12)
B


#Rates2 (Selection gradient)
Pl4R2d<-ggtree(MCT_TRE2, layout = "rectangular", ladderize=TRUE, right=TRUE, position = position_nudge(x = -offset),
               aes(color=rateTK02Brlens3_mean), size=2)%>% flip(63, 73)+
  geom_tiplab(size=3, linesize = 0.01, fontface = "italic", color="black", offset = -285)+
  #geom_tiplab(aes(x=branch, label=round(rateTK02Brlens3_mean, digits = 2)), color="black", vjust=-0.5, size=2, offset = -288)+
  scale_colour_stepsn("Background Rate Threshold",
                      colours = c("blue", "cyan3", "cyan","gray", "orange", "red1" ,"red4"),
                      breaks=c(RelLower2SD,RelLower1SD,RelLower95CI,RelUpper95CI, RelUpper1SD, RelUpper2SD),
                      labels =c("-2 SD", "-1 SD", "Lower 95%CI","Upper 95%CI", "+1 SD", "+2 SD"),
                      limits = c(0,2),
                      guide = "coloursteps",
                      aesthetics = "colour",)+
  #geom_range(range='rateTK02Brlens_0.95HPD', color='red', alpha=.6, size=2) +
  #geom_nodelab(aes(x=branch, label=round(rateTK02Brlens3_mean, digits = 2)), color="black", vjust=-.5, hjust=+.5, size=2, nudge_x = -offset)+
  coord_geo(xlim=c(-450,-260), ylim=c(0,Ntip(MCT_TRE2)+2), expand=FALSE,
            dat = list("epochs", "periods"), abbrv = list(TRUE, FALSE),
            skip = c("Pliocene", "Pleistocene", "Holocene","Quaternary"),
            pos = list("bottom", "bottom"),alpha = 1, height = unit(1, "line"),
            rot = 0, size = list(2.5,3), neg = TRUE) +
  scale_x_continuous(breaks=seq(-450,-260,20), labels=abs(seq(-450,-260,20)))+
  theme_tree2()+
  labs(title = "Selection Strength Gradient")+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position=c(.9, .2),
        legend.title = element_text(size = 10, face = "bold"))
Pl4R2d
Pl4R2d <- revts(Pl4R2d)
Pl4R2d


