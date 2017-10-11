# install.packages(c("ks", "picante"))
# install.packages("TESS")
library(TESS)
library(ks)
library(picante)
setwd("/home/angela/Downloads/treecomp-master/R")
source("sim.R")


t = rand_umt(10)

#If you get the ELF??? error
# GO to the cpp directory as in "treecomp-master/cpp"
# Delete the Rfunc.so file (or rename if you need a backup)
# Run "make"
# Ideally the makefile and the sim.r would checkif 64/32 bit version and load the appropriate library.
# see: http://www.linuxquestions.org/questions/programming-9/how-can-make-makefile-detect-64-bit-os-679513/ 
# and: https://stackoverflow.com/questions/18091614/how-can-i-know-if-r-is-running-on-64-bits-versus-32

sim1=sim(tree=t, sigma=1) #BM
sim1$tval
sim2=sim(tree=t, sigma=1, a=1) #comp model
sim2$tval
sd(sim1$tval)
sd(sim2$tval)
Kcalc(x=sim1$tval, phy=t)
Kcalc(x=sim2$tval, phy=t)

t=rand_umt(50)
d1=sim(t)$tval
d2=sim(t,a=2)$tval

# performs simulations on the tree with varying parameter values
param_stats(tree=t, file="mysims.txt", reps=1000,max_sigma = 5, max_a = 5)
# find the 100 simulations that most closely resemble the observed data, bsaed on summary stats: mean and SD of the distances between neighbouring trait values
lrt(tree=t, data=d1, file="mysims.txt", posteriorSize=100,max_sigma=5, max_a=5) #data with Brownian
lrt(tree=t, data=d2, file="mysims.txt", posteriorSize=100,max_sigma=5, max_a=5) #data with comp model

# find the 100 simulations that most closely resemble the observed data, bsaed on summary stats: 
  # mean and SD of the distances between neighbouring trait values ANS Blomberg's K
param_stats(tree=t, file="mysims2.txt", reps=1000, max_sigma=5,max_a=5, use_K=TRUE)
lrt(tree=t, data=d1, file="mysims2.txt", posteriorSize=100, max_sigma=5, max_a=5, use_K=TRUE)

## should use 1ml simulations and increase the #of accepted simulations
  # 1mil simulations --> 17h


## multivar?
lrt(tree=t, data=as.matrix(cbind(d1,d2)), file="mysims.txt", posteriorSize=100,max_sigma=5, max_a=5) 

read.table("mysims2.txt", sep=",")-> mysims2
hist(mysims2$V4)
plot(mysims2$V4)
