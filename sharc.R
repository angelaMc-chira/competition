rm(list=ls()) #clear

# packages
install.packages("ks")
install.packages("picante")
install.packages("TESS")

# source the functions
source("scripts/R/sim.R")
source("scripts/R/sim2.R")


# read the input
load("data/pc_phylist.RDS") 
target_tree<- pcphy_list$Anatinae #vary as I apply it to more groups
file<- "Anatinae" #vary as I apply it to more groups 
load("data/sig.RDS") ; load("data/atry.RDS")


# run the model
for(jj in 1: 10) #trial first 10
  {d = sim(tree=target_tree, a=atry[jj], sigma=sig[jj], ntraits=8)$tval 
   stats = summary_stats(data=d, ntraits=8, use_K = T, tree=target_tree) 
   write(c(sig[jj], atry[jj], stats), file=file, append=TRUE, sep=",")
  }

