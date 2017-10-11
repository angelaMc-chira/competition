
# ICEBERG 32/64?
## How to run via command line
#  Rscript 35_IcebergFxn.R -i 1

# errors because of packages 
# ERROR: dependency ‘rgl’ is not available for package ‘ks’
# * removing ‘/home/angela/Downloads/Iceberg_docs/libs/ks’

# also : iceberg 32 or 64 bits operating system


## How do I expect the bash script to look like (call it script.sh)
# 
## #!/bin/bash
##
## -t 1-100000
##
# echo "Task id is $SGE_TASK_ID"
# Rscript 35_IcebergFxn.R -i "$SGE_TASK_ID"

## How do I expect the iceberg job to be submitted
# qsub script.sh


## 35_IcebergFxn.R:

# install packages needed
.libPaths("libs/")
install.packages("TESS", repos = "http://cran.us.r-project.org")
library(TESS, lib.loc="libs/")
# Dependency for ks
install.packages("ks", repos = "http://cran.us.r-project.org", dependencies=TRUE)
library(ks, lib.loc="libs/")
install.packages("picante", repos = "http://cran.us.r-project.org")
library(picante, lib.loc="libs/")
install.packages("optparse", repos = "http://cran.us.r-project.org")
library(optparse, lib.loc="libs/")
option_list = list(make_option(c("-i", "--i"), type="character", default=NULL,
                               help="The i", metavar="character"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser);

i<-opt$i
print(i)


# source the functions
home_path<- getwd()
setwd("treecomp-master/R")
source("sim.R")
source("sim2.R")
setwd(home_path)


# read the input
load("pc_phylist.RDS") # load("pc_dflist.RDS")

target_tree<- pcphy_list$Anserinae  # parameter that needs to change with different groups ...
#sig = runif(1e5, 0, 0.05) ; save(sig, file="sig.RDS")
#atry = runif(1e5, 0, 8) ; save(atry, file="atry.RDS")

load("sig.RDS") ; load("atry.RDS")
file<- "Anserinae" # parameter that needs to change with different groups ...


# run the model, uses i as a variable that changes
if (i==1) dir.create(paste("dir",file,"/", sep="")) #creates it in the working dir
d = sim(tree=target_tree, a=atry[i], sigma=sig[i], ntraits=8)$tval 
stats = summary_stats(data=d, ntraits=8, use_K = T, tree=target_tree) 
write(c(sig[i], atry[i], stats), file=paste(paste("dir",file,"/", sep=""), file,"_PS", i , sep=""), append=F, sep=",")

