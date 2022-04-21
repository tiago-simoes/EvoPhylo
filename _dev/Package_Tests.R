library(devtools)

devtools::install_github("tiago-simoes/EvoPhylo", lib="C:/Program Files/R/R-4.1.0/library")

library(EvoPhylo)

# TESTS
devtools::check_man("E:/Git/EvoPhylo")


#Build package manual
build_manual(pkg = "E:/Git/EvoPhylo", path = "E:/Git/EvoPhylo")

