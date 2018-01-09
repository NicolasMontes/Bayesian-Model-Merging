rm(list=ls())
gc()
library("gtools")
library("HMM")


load("RData/strings.RData")


source(file = "../functions/m0.R")
hmm=m0(f)


source(file = "../functions/bmm.R")

#Profiling
Rprof("profiling.out")
(hmm2=bmm(hmm,f,10))
Rprof()
summaryRprof("profiling.out")




