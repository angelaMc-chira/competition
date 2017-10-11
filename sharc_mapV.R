rm(list=ls()) #clear


# source the functions - use relative path ...?
setwd("/home/angela/Downloads/treecomp-master/R")
source("sim.R")
source("sim2.R")
setwd("/home/angela/Downloads/Iceberg_docs")


# read the input
load("pc_phylist.RDS") 
target_tree<- pcphy_list$Anatinae #vary as I apply it to more groups
file<- "Anatinae" #vary as I apply it to more groups 
load("sig.RDS") ; load("atry.RDS")
#i from 0 to 99
sig<- sig[i*1000+1:(i+1)*1000]
atry<- atry[i*1000+1:(i+1)*1000]
file<- paste("Anatinae", i, sep="")

# run the model
for(jj in 1: length(sig)) 
  {d = sim(tree=target_tree, a=atry[jj], sigma=sig[jj], ntraits=8)$tval 
   stats = summary_stats(data=d, ntraits=8, use_K = T, tree=target_tree) 
   write(c(sig[jj], atry[jj], stats), file=file, append=TRUE, sep=",")
  }

