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



##################### TIP+Node MultiCons Plots ##  MrBayes MCT tree #########################

# For plotting trees without dummy extant, use the following only (not remove.MrB.extant) - SINGLE & MULTI CLOCK

MCT_TRE2<- read.mrbayes("3p_TK02_FT_SFBD4l_TipNode_MultiTopCons_AllCom.t")
names(MCT_TRE2@data)

#Drop extant "dummy" tip
MCT_TRE2 <- drop.tip(MCT_TRE2, "Dummyextant")

MCT_TRE2@data$`prob+-sd` <- as.factor(MCT_TRE2@data$`prob+-sd`)
MCT_TRE2@data <- MCT_TRE2@data %>% mutate_if(is.character,as.numeric)

offset = min(MCT_TRE2@data$age_median)
offset

center <- MCT_TRE2@data$age_median


######################### Mr Bayes RATE STATS #######################

### Median Rates by clade (Summary stats)

nodes <- MCT_TRE2@data$node
rates1 <- MCT_TRE2@data$rateTK02Brlens1_median
rates2 <- MCT_TRE2@data$rateTK02Brlens2_median
rates3 <- MCT_TRE2@data$rateTK02Brlens3_median

RateTable <- data.frame(nodes, rates1, rates2, rates3)

write.csv(RateTable, file="RateTable.csv")

#Import new RateTable with customized clade names included -> named here RateTable_Medians_Clades.csv (example provided)

CladeRate1 <- read.csv("RateTable_Medians_Clades.csv")

#Summary stats (medians)

Sum_CladeRate1<-Summarize(rates1 ~ clade,
                          data=CladeRate1,
                          digits=3)
Sum_CladeRate1

Sum_CladeRate2<-Summarize(rates2 ~ clade,
                          data=CladeRate1,
                          digits=3)
Sum_CladeRate2

Sum_CladeRate3<-Summarize(rates3 ~ clade,
                          data=CladeRate1,
                          digits=3)
Sum_CladeRate3

AllCladeRates1 <- cbind (Sum_CladeRate1, Sum_CladeRate2, Sum_CladeRate3)

write.csv(AllCladeRates, file="Sum_AllCladeRates_Medians.csv")

### Mean Rates by clade (Summary stats)

nodesM <- MCT_TRE2@data$node
ratesM1 <- MCT_TRE2@data$rateTK02Brlens1_mean
ratesM2 <- MCT_TRE2@data$rateTK02Brlens2_mean
ratesM3 <- MCT_TRE2@data$rateTK02Brlens3_mean

RateTable2 <- data.frame(nodesM, ratesM1, ratesM2, ratesM3)

write.csv(RateTable2, file="RateTable2.csv")

#Import new RateTable2 with customized clade names included -> named here RateTable_Means_Clades.csv (example provided)

CladeRate2 <- read.csv("RateTable_Means_Clades.csv")

#Summary stats (means)
Sum_CladeRateM1<-Summarize(rates1 ~ clade,
                          data=CladeRate2,
                          digits=3)
Sum_CladeRateM1

Sum_CladeRateM2<-Summarize(rates2 ~ clade,
                          data=CladeRate2,
                          digits=3)
Sum_CladeRateM2

Sum_CladeRateM3<-Summarize(rates3 ~ clade,
                          data=CladeRate2,
                          digits=3)
Sum_CladeRateM3

AllCladeRates2 <- cbind (Sum_CladeRateM1, Sum_CladeRateM2, Sum_CladeRateM3)

write.csv(AllCladeRates2, file="Sum_AllCladeRates_Means.csv")


################ Rates by morpho clock: medians
# Summary stats   


Sum_RateByClock1<-Summarize(rates1,
                            data=CladeRate1,
                            digits=3)
Sum_RateByClock1

Sum_RateByClock2<-Summarize(rates2,
                            data=RateByClock,
                            digits=3)
Sum_RateByClock2

Sum_RateByClock3<-Summarize(rates3,
                            data=RateByClock,
                            digits=3)
Sum_RateByClock3

AllRateByClocks <- cbind (Sum_RateByClock1, Sum_RateByClock2, Sum_RateByClock3)
AllRateByClocks <- as.data.frame(AllRateByClocks)
names(AllRateByClocks)

write.csv(AllRateByClocks, file="Sum_AllRateByClocks_medians.csv")


################ Rates by morpho clock: means

# Summary stats   


Sum_RateByClockM1<-Summarize(rates1,
                            data=CladeRate2,
                            digits=3)
Sum_RateByClockM1

Sum_RateByClockM2<-Summarize(rates2,
                            data=CladeRate2,
                            digits=3)
Sum_RateByClockM2

Sum_RateByClockM3<-Summarize(rates3,
                            data=CladeRate2,
                            digits=3)
Sum_RateByClockM3

AllRateByClocks2 <- cbind (Sum_RateByClockM1, Sum_RateByClockM2, Sum_RateByClockM3)
AllRateByClocks2 <- as.data.frame(AllRateByClocks2)
names(AllRateByClocks2)

write.csv(AllRateByClocks2, file="Sum_AllRateByClocks_means.csv")



# Clock and clade rate distribution


#Melt DATA-medians

CladeRate1_Melted<-melt(CladeRate1, measure.vars = c(3,4,5),
                        variable.name = "clock", value.name = "rate")

write.csv(CladeRate1_Melted, file = "RateByClockrMelted_medians.csv")

#Melt DATA-Means

CladeRate2_Melted<-melt(CladeRate2, measure.vars = c(3,4,5),
                        variable.name = "clock", value.name = "rate")

write.csv(CladeRate2_Melted, file = "RateByClockrMelted_means.csv")


#Hist (by clock)
Pl1a<- ggplot(CladeRate1_Melted, aes(x=rate, fill = clock, color = clock))+
  geom_density(position="dodge",  alpha = 0.6)+
  scale_x_continuous()+
  #geom_jitter(aes(y = nodes), height = 0.01)+
  #scale_x_continuous()+
  #scale_y_continuous()+
  #geom_text_repel(aes(label=nodes))+
  theme(legend.position="bottom")
Pl1a

#Hist (by clade)
Pl1b<- ggplot(CladeRate1_Melted, aes(x=rate, fill = clade, color = clade))+
  geom_density(alpha = 0.3)+
  scale_x_continuous()+
  #geom_jitter(aes(y = nodes), height = 0.01)+
  #scale_x_continuous()+
  #scale_y_continuous(trans = "sqrt")+
  #geom_text_repel(aes(label=nodes))+
  theme(legend.position="top")
Pl1b



# Clade rate by clock distribution (stacked Histo)

Pl2a<- ggplot(CladeRate1, aes(x=rates1, fill = clade, color = clade))+
  geom_density(position = "stack", alpha = 1, weight = 2)+
  scale_x_continuous()+
  #geom_jitter(aes(y = nodes), height = 0.01)+
  #scale_x_continuous()+
  #scale_y_continuous(trans = "sqrt")+
  #geom_text_repel(aes(label=nodes))+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  #theme_classic()+
  theme(legend.position="none")
Pl2a

Pl2b<- ggplot(CladeRate1, aes(x=rates2, fill = clade, color = clade))+
  geom_density(position = "stack", alpha = 1, weight = 2)+
  scale_x_continuous()+
  #geom_jitter(aes(y = nodes), height = 0.01)+
  #scale_x_continuous()+
  #scale_y_continuous(trans = "sqrt")+
  #geom_text_repel(aes(label=nodes))+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  #theme_classic()+
  theme(legend.position="bottom")
Pl2b

Pl2c<- ggplot(CladeRate1, aes(x=rates3, fill = clade, color = clade))+
  geom_density(position = "stack", alpha = 1, weight = 2)+
  scale_x_continuous()+
  #geom_jitter(aes(y = nodes), height = 0.01)+
  #scale_x_continuous()+
  #scale_y_continuous(trans = "sqrt")+
  #geom_text_repel(aes(label=nodes))+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  #theme_classic()+
  theme(legend.position="none")

Pl2c

#Combine plots
par(mar = c(7, 7, 7, 7))
A2<-plot_grid(Pl2a, Pl2c, Pl2b, ncol = 1)  
A2


# Clock rates linear regressions

Lm2a<- lm(data=CladeRate1, rates1 ~ rates2)
summary(Lm2a)


R2a<-ggplot(CladeRate1, aes(y = rates1, x=rates2))+
  geom_point()+
  geom_smooth(method ="lm", se=TRUE)+
  scale_x_continuous()+
  scale_y_continuous()+
  labs (x="Skull1 rates", y="Postcranial rates")+
  theme_classic()
R2a

Lm2b<- lm(data=CladeRate1, rates1 ~ rates3)
summary(Lm2b)

R2b<-ggplot(CladeRate1, aes(y = rates1, x=rates3))+
  geom_point()+
  geom_smooth(method ="lm", se=TRUE)+
  scale_x_continuous()+
  scale_y_continuous()+
  labs (x="Skull1 rates", y="Skull2 rates")+
  theme_classic()
R2b

Lm2c<- lm(data=CladeRate1, rates2~rates3)
summary(Lm2c)


R2c<-ggplot(CladeRate1, aes(y = rates2, x=rates3))+
  geom_point()+
  geom_smooth(method ="lm", se=TRUE)+
  scale_x_continuous()+
  scale_y_continuous()+
  labs (x="Postcranium rates", y="Skull2 rates")+
  theme_classic()
R2c

#Combine plots
par(mar = c(7, 7, 7, 7))
A1<-plot_grid(R2a, R2b, R2c, ncol = 3)  
A1
