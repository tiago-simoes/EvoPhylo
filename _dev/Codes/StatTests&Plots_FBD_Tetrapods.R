library(ggplot2)
library(cowplot)
library(plotly)
library(ggrepel)
library(plyr)
library(lattice)
library(Rmisc)
library(deeptime)
library(reshape2)
library(FSA)
library(car)
library(tidyr)
library(rstatix)


# Import combined log file of both tip+TipNode analyses 
AllRunsrMelted_MC <- read.csv("AllRuns_COMB_Melted.csv")
names(AllRunsrMelted_MC)

#Downsample posterior sample to 1k (necessary for some stats tests)
Down_AllRunsrMelted_MC<- AllRunsrMelted_MC[seq(1, nrow(AllRunsrMelted_MC), 270),]

write.csv(Down_AllRunsrMelted_MC, file= "Down1k_1p_Tip_Par_Melted.csv")

#Order groups by Time_bin

AllRunsrMelted_MC$Time_bin = factor(AllRunsrMelted_MC$Time_bin,
                                   levels=unique(AllRunsrMelted_MC$Time_bin))

summary(AllRunsrMelted_MC)



####### Density plot 
#clockrate tip vs tip+node

Pl1a<- ggplot(AllRunsrMelted_MC, aes(clockrate, fill = analysis)) +
  geom_density(position="dodge",  alpha = 0.6)+
  scale_x_continuous(limits=c(0,0.075))+
  theme(legend.position="top")
Pl1a




########## Tests for: FBD parameters #################
  
#Summary stats FBD

Sum_net_speciation<-Summarize(clockrate ~ Time_bin,
                           data=AllRunsrMelted_MC,
                           digits=3)
parameter <- c("net_speciation")
Sum_net_speciation <- cbind(parameter, Sum_net_speciation)

Sum_relative_extinction<-Summarize(relative_extinction ~ Time_bin,
                              data=AllRunsrMelted_MC,
                              digits=3)
parameter <- c("relative_extinction")
Sum_relative_extinction <- cbind(parameter, Sum_relative_extinction)

Sum_relative_fossilization<-Summarize(relative_fossilization ~ Time_bin,
                              data=AllRunsrMelted_MC,
                              digits=3)
parameter <- c("relative_fossilization")
Sum_relative_fossilization <- cbind(parameter, Sum_relative_fossilization)

Sum_FBD <- rbind (Sum_net_speciation,
                  Sum_relative_extinction,
                  Sum_relative_fossilization)

write.csv(Sum_FBD, file="Sum_FBD_MC.csv")


#Remove NA by col
AllRunsrMelted_MC<-drop_na(AllRunsrMelted_MC, net_speciation)

### Normality of distributions by group


####### Density plot for each group 
#net_speciation
Pl1a<- ggplot(AllRunsrMelted_MC, aes(net_speciation, fill = Time_bin)) +
  geom_density(position="dodge",  alpha = 0.6)+
  scale_x_continuous(limits=c(0,0.25))+
  theme(legend.position="top")
Pl1a

Pl1b<-histogram(~ net_speciation | Time_bin,
          data=AllRunsrMelted_MC,
          layout=c(1,3))
Pl1b

par(mar = c(10, 10, 10, 10))
A1<-plot_grid(Pl1a,Pl1b,ncol=2)  
A1

#relative_extinction
Pl1c<- ggplot(AllRunsrMelted_MC, aes(relative_extinction, fill = Time_bin)) +
  geom_density(position="dodge", alpha = 0.6)+
  scale_x_continuous(limits=c(0.5,1))+
  theme(legend.position="top")
Pl1c

Pl1d<-histogram(~ relative_extinction | Time_bin,
                data=AllRunsrMelted_MC,
                layout=c(1,3))
Pl1d

par(mar = c(10, 10, 10, 10))
A2<-plot_grid(Pl1c,Pl1d,ncol=2)  
A2

#relative_fossilization

Pl1e<- ggplot(AllRunsrMelted_MC, aes(relative_fossilization, fill = Time_bin)) +
  geom_density(position="dodge", alpha = 0.6)+
  scale_x_continuous(limits=c(0,0.15))+
  theme(legend.position="top")
Pl1e

Pl1f<-histogram(~ relative_fossilization | Time_bin,
                data=AllRunsrMelted_MC,
                layout=c(1,3))
Pl1f


par(mar = c(10, 10, 10, 10))
B1<-plot_grid(Pl1a,Pl1c,Pl1e, ncol=3)  
B1

par(mar = c(10, 10, 10, 10))
B2<-plot_grid(Pl1a,Pl1c,Pl1e, ncol=3)  
B2

par(mar = c(10, 10, 10, 10))
B12<-plot_grid(B1,B2, nrow=2)  
B12



############### Stats for final analysis (TipNode_1partition)

######Import log files
AllRuns <- read.table("1p_TipNode_MC_10k.p", header = TRUE)
names(AllRuns)

#Melt DATA 


AllRunsrMelted1<-melt(AllRuns, measure.vars = c(8,9,10),
                      variable.name = "Time_bin", value.name = "net_speciation")
AllRunsrMelted1$Time_bin<- sub("net_speciation_", "", AllRunsrMelted1$Time_bin)

AllRunsrMelted2<-melt(AllRuns, measure.vars = c(12,13,14),
                      variable.name = "Time_bin", value.name = "relative_extinction")
AllRunsrMelted2$Time_bin<- sub("relative_extinction_", "", AllRunsrMelted2$Time_bin)

AllRunsrMelted3<-melt(AllRuns, measure.vars = c(16,17,18),
                      variable.name = "Time_bin", value.name = "relative_fossilization")
AllRunsrMelted3$Time_bin<- sub("relative_fossilization_", "", AllRunsrMelted3$Time_bin)


AllRunsrMelted <- cbind (AllRunsrMelted1, 
                         relative_extinction= AllRunsrMelted2$relative_extinction, 
                         relative_fossilization = AllRunsrMelted3$relative_fossilization)

write.csv(AllRunsrMelted, file = "1p_TipNode_MC_Melted.csv")


#ViolinPlots - FOCAL RANGE 
#Use AllRunsrMelted object 

#ViolinPlots
Pl1a<- ggplot(AllRunsrMelted, aes(x=Time_bin, y=net_speciation)) + 
  geom_violin(alpha=0.8, colour="blue", fill = "cyan3", trim=TRUE, draw_quantiles = 0.5, size=0.8, )+
  ylim(0, 0.25)+
  stat_summary(fun=mean, geom="point", shape=16, size=2, color = "black")+
  labs(title = "Net speciation", x="Time bin", y = "Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="bottom")
Pl1a

Pl1b<- ggplot(AllRunsrMelted, aes(x=Time_bin, y=relative_extinction)) + 
  geom_violin(alpha=0.5, colour="firebrick", fill = "red", trim=TRUE, draw_quantiles = 0.5, size=0.8)+
  ylim(0.5, 1)+
  stat_summary(fun=mean, geom="point", shape=16, size=2, color = "black")+
  labs(title = "Relative extinction", x="Time bin", y = "Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="bottom")
Pl1b

Pl1c<- ggplot(AllRunsrMelted, aes(x=Time_bin, y=relative_fossilization)) + 
  geom_violin(alpha=0.6, colour = "green4", fill = "green", trim=TRUE, draw_quantiles = 0.5, size=0.8)+
  ylim(0, 0.15)+
  stat_summary(fun=mean, geom="point", shape=16, size=2, color = "black")+
  labs(title = "Relative fossilization", x="Time bin", y = "Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="bottom")
Pl1c

par(mar = c(10, 10, 10, 10))
B1<-plot_grid(Pl1a,Pl1b,Pl1c, ncol=3)  
B1

par(mar = c(10, 10, 10, 10))
B2<-plot_grid(Pl1a,Pl1b,Pl1c, ncol=3)  
B2

par(mar = c(10, 10, 10, 10))
B12<-plot_grid(B1,B2, nrow=2)  
B12




#ViolinPlots - FULL RANGE

Pl1a<- ggplot(AllRunsrMelted, aes(x=Time_bin, y=net_speciation)) + 
  geom_violin(alpha=0.8, colour="blue", fill = "cyan3", trim=FALSE, draw_quantiles = 0.5, size=0.8, )+
  ylim(0, 1)+
  stat_summary(fun=mean, geom="point", shape=16, size=2, color = "black")+
  #stat_summary(fun=mean_sdl, geom="pointrange", color="black")+
  #stat_summary(fun=median, geom="point", shape = 17, size=2, color="purple")+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
  #geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(title = "Net speciation", x="Time slices", y = "Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="bottom")
Pl1a

Pl1b<- ggplot(AllRunsrMelted, aes(x=Time_bin, y=relative_extinction)) + 
  geom_violin(alpha=0.5, colour="firebrick", fill = "red", trim=FALSE, draw_quantiles = 0.5, size=0.8)+
  ylim(0, 1)+
  stat_summary(fun=mean, geom="point", shape=16, size=2, color = "black")+
  #stat_summary(fun=mean_sdl, geom="pointrange", color="black")+
  #stat_summary(fun=median, geom="point", shape = 17, size=2, color="purple")+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
  #geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(title = "Relative extinction", x="Time slices", y = "Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="bottom")
Pl1b

Pl1c<- ggplot(AllRunsrMelted, aes(x=Time_bin, y=relative_fossilization)) + 
  geom_violin(alpha=0.6, colour = "green4", fill = "green", trim=FALSE, draw_quantiles = 0.5, size=0.8)+
  ylim(0, 1)+
  stat_summary(fun=mean, geom="point", shape=16, size=2, color = "black")+
  #stat_summary(fun=mean_sdl, geom="pointrange", color="black")+
  #stat_summary(fun=median, geom="point", shape = 17, size=2, color="purple")+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
  #geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(title = "Relative fossilization", x="Time slices", y = "Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="bottom")
Pl1c

par(mar = c(10, 10, 10, 10))
B1<-plot_grid(Pl1a,Pl1b,Pl1c, ncol=3)  
B1

par(mar = c(10, 10, 10, 10))
A1<-plot_grid(Pl1a,Pl1b,Pl1c, ncol=3)  
A1







#Shapiro-Wilk normality test

#Use downsampled file Down_AllRunsrMelted_MC

#Combined
shapiro.test(Down_AllRunsrMelted_MC$net_speciation)
shapiro.test(Down_AllRunsrMelted_MC$relative_extinction)
shapiro.test(Down_AllRunsrMelted_MC$relative_fossilization)

#By group/levels
tapply(Down_AllRunsrMelted_MC$net_speciation, Down_AllRunsrMelted_MC$Time_bin, shapiro.test)
tapply(Down_AllRunsrMelted_MC$relative_extinction, Down_AllRunsrMelted_MC$Time_bin, shapiro.test)
tapply(Down_AllRunsrMelted_MC$relative_fossilization, Down_AllRunsrMelted_MC$Time_bin, shapiro.test)


### Homogeneity of variance (homoscedasticity)-MOST IMPORTANT
#Plot of residuals vs. fitted values (Normality of residuals and homoscedasticity)
#Define a linear model

model = lm(net_speciation ~ Time_bin,
           data = Down_AllRunsrMelted_MC)
model
summary(model)
model$residuals
plot(model)

#Add residuals to object
Down_AllRunsrMelted_MC$Residuals <- model$residuals

#Bartlett's test for homogeneity of variance
bartlett.test(Residuals ~ Time_bin, data=Down_AllRunsrMelted_MC)

#Fligner-Killeen test for homogeneity of variance (robust to departures in normality of the data)

fligner.test(Residuals ~ Time_bin, data=Down_AllRunsrMelted_MC)


### Normality of residuals ()

#Plot1
ggplot(Down_AllRunsrMelted_MC, aes(Residuals, fill = Time_bin)) +
  geom_density(position="dodge",
               alpha = 0.6)

histogram(~ Residuals | Time_bin,
          data=Down_AllRunsrMelted_MC,
          layout=c(1,3))

#Shapiro-Wilk normality test
#Residuals
shapiro.test(residuals(model))

#By group/levels
tapply(Down_AllRunsrMelted_MC$Residuals, Down_AllRunsrMelted_MC$Time_bin, shapiro.test)
#
#




#### Tests across time bins (within analysis)

### Pairwise t-test

PWST_NetSpec <- AllRunsrMelted_MC %>% 
  group_by(analysis) %>%
  pairwise_t_test(net_speciation ~ Time_bin, 
                  paired = FALSE, p.adjust.method = "fdr")

PWST_NetSpec


PWST_RelExt <- AllRunsrMelted_MC %>% 
  group_by(analysis) %>%
  pairwise_t_test(relative_extinction ~ Time_bin, 
                  paired = FALSE, p.adjust.method = "fdr")
PWST_RelExt

PWST_RelFoss <- AllRunsrMelted_MC %>% 
  group_by(analysis) %>%
  pairwise_t_test(relative_fossilization ~ Time_bin, 
                  paired = FALSE, p.adjust.method = "fdr")
PWST_RelFoss

PWST_TimeBins_Table <- rbind (PWST_NetSpec, 
                        PWST_RelExt, 
                        PWST_RelFoss)
PWST_TimeBins_Table

write.csv(PWST_MC_Table, file="PWST_AllAna_Timebin.csv")


### Pairwise Wilcoxon tests (Mann-Whitney)

PWWI_NetSpec <- AllRunsrMelted_MC %>% 
  group_by(analysis) %>%
  pairwise_wilcox_test(net_speciation ~ Time_bin, 
                       paired = FALSE, p.adjust.method = "fdr")

PWWI_NetSpec


PWWI_RelExt <- AllRunsrMelted_MC %>% 
  group_by(analysis) %>%
  pairwise_wilcox_test(relative_extinction ~ Time_bin, 
                       paired = FALSE, p.adjust.method = "fdr")
PWWI_RelExt

PWWI_RelFoss <- AllRunsrMelted_MC %>% 
  group_by(analysis) %>%
  pairwise_wilcox_test(relative_fossilization ~ Time_bin, 
                       paired = FALSE, p.adjust.method = "fdr")
PWWI_RelFoss

PWWI_AllAna_Timebin <- rbind (PWWI_NetSpec, 
                        PWWI_RelExt, 
                        PWWI_RelFoss)
PWWI_AllAna_Timebin

write.csv(PWWI_AllAna_Timebin, file="PWWI_AllAna_Timebin.csv")







####################### Tests across analyses (same time bins)

RunsTip <- read.table("1p_Tip_10k.p", header = TRUE)
RunsMC <- read.table("1p_TipNode_MC_10k.p", header = TRUE)

analysis <- c("Tip")
RunsTip2 <- cbind(analysis, RunsTip)

analysis <- c("TipNode")
RunsMC2  <- cbind(analysis, RunsMC)

AllRuns_COMB <- rbind(RunsTip2, RunsMC2)
names(AllRuns_COMB)

write.csv(AllRuns_COMB, file = "AllRuns_COMB.csv")




#Summary stats (FBD by analysis and time bin)

summary(AllRunsrMelted_MC)

Sum_net_speciation<-Summarize(net_speciation ~ analysis*Time_bin,
                               data=AllRunsrMelted_MC,
                               digits=3)
Sum_net_speciation

Sum_relative_extinction<-Summarize(relative_extinction ~ analysis*Time_bin,
                               data=AllRunsrMelted_MC,
                               digits=3)
Sum_relative_extinction

Sum_relative_fossilization<-Summarize(relative_fossilization ~ analysis*Time_bin,
                               data=AllRunsrMelted_MC,
                               digits=3)
Sum_relative_fossilization

SumStats_FBD_AllAna <- rbind (net_speciation = Sum_net_speciation, 
                        relative_extinction = Sum_relative_extinction, 
                        relative_fossilization = Sum_relative_fossilization)
SumStats_FBD_AllAna

write.csv(SumStats_FBD_AllAna, file="SumStats_FBD_AllAna.csv")
          


### Pairwise t-test

PWST_NetSpec <- AllRunsrMelted_MC %>% 
  group_by(Time_bin) %>%
  pairwise_t_test(net_speciation ~ analysis, 
                  paired = FALSE, p.adjust.method = "fdr")

PWST_NetSpec


PWST_RelExt <- AllRunsrMelted_MC %>% 
  group_by(Time_bin) %>%
  pairwise_t_test(relative_extinction ~ analysis, 
                  paired = FALSE, p.adjust.method = "fdr")
PWST_RelExt

PWST_RelFoss <- AllRunsrMelted_MC %>% 
  group_by(Time_bin) %>%
  pairwise_t_test(relative_fossilization ~ analysis, 
                  paired = FALSE, p.adjust.method = "fdr")
PWST_RelFoss

PWST_MC_Table <- rbind (PWST_NetSpec, 
                        PWST_RelExt, 
                        PWST_RelFoss)
PWST_MC_Table

write.csv(PWST_MC_Table, file="PWST_AllAna_Table.csv")


### Pairwise Wilcoxon tests (Mann-Whitney)

PWWI_NetSpec <- AllRunsrMelted_MC %>% 
  group_by(Time_bin) %>%
  pairwise_wilcox_test(net_speciation ~ analysis, 
                  paired = FALSE, p.adjust.method = "fdr")

PWWI_NetSpec


PWWI_RelExt <- AllRunsrMelted_MC %>% 
  group_by(Time_bin) %>%
  pairwise_wilcox_test(relative_extinction ~ analysis, 
                  paired = FALSE, p.adjust.method = "fdr")
PWWI_RelExt

PWWI_RelFoss <- AllRunsrMelted_MC %>% 
  group_by(Time_bin) %>%
  pairwise_wilcox_test(relative_fossilization ~ analysis, 
                  paired = FALSE, p.adjust.method = "fdr")
PWWI_RelFoss

PWWI_MC_Table <- rbind (PWWI_NetSpec, 
                        PWWI_RelExt, 
                        PWWI_RelFoss)
PWWI_MC_Table

write.csv(PWWI_MC_Table, file="PWWI_AllAna_Table.csv")

