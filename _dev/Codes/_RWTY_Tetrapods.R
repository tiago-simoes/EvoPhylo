library(ggplot2)
library(rwty)
library(reshape2)
library(cowplot)


rwty.processors <<- 8

#Import all .p and .t files from appropriate directory 
Tetra <- load.multi("D:\\Programas\\Rcodes\\Phylogenies\\PhyloParameters\\TetrapodsTipNode_MultiCons", format = "mb")



#### GENERAL WRAPPER
Tetra.RWTY <- analyze.rwty(Tetra, burnin=25, fill.color = 'net_speciation')


# to see which plots you have
names(Tetra.RWTY)

Tetra.RWTY$LnL.trace
Tetra.RWTY$LnPr.trace
Tetra.RWTY$prop_ancfossil.trace
Tetra.RWTY$net_speciation.trace
Tetra.RWTY$igrvar.trace
Tetra.RWTY$clockrate.trace
Tetra.RWTY$topology.trace.plot



#SPECIFICS 

names(Tetra$`1p_LN-Sym_TK02_FT_SFBD4l_TipNode_MultiTopCons.nex.run1.t`)

makeplot.param(Tetra$`1p_LN-Sym_TK02_FT_SFBD4l_TipNode_MultiTopCons.nex.run1.t`, burnin = 25, 
               parameter = "LnL", facet = TRUE, free_y = TRUE)

Tip_Run1_LnLClock <-makeplot.pairs(Tetra$`1p_LN-Sym_TK02_FT_SFBD4l_TipNode_MultiTopCons.nex.run1.t`, burnin = 50, 
            strip = c(1,3,4,5,6), params = c("LnL", "tk02var", "clockrate"))
Tip_Run1_LnLClock 

Tip_Run1_FBD1 <-makeplot.pairs(Tetra$`1p_LN-Sym_TK02_FT_SFBD4l_TipNode_MultiTopCons.nex.run1.t`, burnin = 50, 
           strip = c(1,3,4,5,6), params = c("net_speciation_1", "relative_extinction_1", "relative_fossilization_1"))
Tip_Run1_FBD1

Tip_Run1_FBD2 <-makeplot.pairs(Tetra$`1p_LN-Sym_TK02_FT_SFBD4l_TipNode_MultiTopCons.nex.run1.t`, burnin = 50, 
          strip = c(1,3,4,5,6), params = c("net_speciation_2", "relative_extinction_2", "relative_fossilization_2"))
Tip_Run1_FBD2

Tip_Run1_FBD3 <-makeplot.pairs(Tetra$`1p_LN-Sym_TK02_FT_SFBD4l_TipNode_MultiTopCons.nex.run1.t`, burnin = 50, 
           strip = c(1,3,4,5,6), params = c("net_speciation_3", "relative_extinction_3", "relative_fossilization_3"))
Tip_Run1_FBD3


