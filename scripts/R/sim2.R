
### ALL IN ONE GO, PARALLEL
param_stats_par = function(tree, file='param_stats.out', reps=1e3, max_sigma=8, max_a=8, ncores=30, save_sim=F,...)
{
  sig = runif(reps, 0, max_sigma)
  atry = runif(reps, 0, max_a)
  
  dir.create(paste("dir",file,"/", sep="")) #creates it in the working dir
  
  rez<- mclapply(1:reps, function(i){
    d = sim(tree=tree, a=atry[i], sigma=sig[i], ...)$tval
    if(save_sim==T) { d1 = rbind( c( sig[i],atry[i], rep(-1, (length(d[1,])-2))),d) #add sigma and comp_a
                     write.table(d1, file=paste(paste("dir",file,"/", sep=""),file, "_", i, sep=""), sep=",", row.names = F, col.names = F)
                     rm(d1)
    }     
    
    stats = summary_stats(data=d, tree=tree,...) 
    write(c(sig[i], atry[i], stats), file=paste(paste("dir",file,"/", sep=""), file,"_PS", i , sep=""), append=F, sep=",")
  }, mc.cores=ncores, mc.preschedule = F)
  
  # add a for that reads each of those lines, appends them to final_txt, and deletes the file with 1 line
  for(i in 1: reps)
  {  as.numeric(read.table (paste(paste("dir",file,"/", sep=""), file,"_PS", i , sep=""), sep=",")) ->l
    write(l, file=file, append=TRUE, sep=",")
    rm(l)
    file.remove(paste(paste("dir",file,"/", sep=""), file,"_PS", i , sep=""))
  }}
#param_stats_par(tree=pcphy_list$Anserinae, file="test_small", reps=10, max_sigma=8, max_a=8, ntraits=8, cores=10, save_sim=T) 
#lrt_acc(tree=pcphy_list$Anserinae, data=data_use, file="test_small", posteriorSize=5,max_sigma=8, max_a=8, ntraits=8) ->acc
##-------------------------------------------------------------------------------------------------------------------------------------------------------


### SIM FIRST
param_stats_simonly = function(tree, file='param_stats.out', reps=1e3, max_sigma=8, max_a=8,...) 
{ sig = runif(reps, 0, max_sigma)
  atry = runif(reps, 0, max_a)

  dir.create(paste("dir",file,"/", sep=""))
  
for(i in 1:reps)
  {d = sim(tree=tree, a=atry[i], sigma=sig[i], ...)$tval
   d = rbind( c( sig[i],atry[i], rep(-1, (length(d[1,])-2))),d) #add sigma and comp_a
   write.table(d, file=paste(paste("dir",file,"/", sep=""),file, "_", i, sep=""), sep=",", row.names = F, col.names = F)
}}
#param_stats_simonly(tree=phy_list$ACCIPITRIFORMES, file="test_acc4", reps=2, max_sigma=8, max_a=8, ntraits=8) #or add ntraits to both


param_stats_simonly_parr = function(tree, file='param_stats.out', reps=1e3, max_sigma=8, max_a=8, ncores=30,...)
{
  sig = runif(reps, 0, max_sigma)
  atry = runif(reps, 0, max_a)
  dir.create(paste("dir",file,"/", sep=""))
  
  rez<- mclapply(1:reps, function(i){
    #d = sim(tree=tree, a=atry[i], sigma=sig[i], ...)$tval
    d = sim(tree=tree, a=atry[i], sigma=sig[i], ...)$tval
    d = rbind( c(sig[i],atry[i], rep(-1, (length(d[1,])-2))),d) #add sigma and comp_a
    write.table(d, file=paste(paste("dir",file,"/", sep=""),file, "_", i, sep=""), sep=",", row.names = F, col.names = F)
  }, mc.cores=ncores, mc.preschedule = F)
}
#param_stats_simonly_parr(tree=phy_list$ACCIPITRIFORMES, file="test_acc3", reps=30, max_sigma=8, max_a=8, ntraits=8, ncores=10) #or add ntraits to both
##-------------------------------------------------------------------------------------------------------------------------------------------------------


# SUM STATS
param_stats_sstat_only = function(file,tree,...)
{# for each dataset
  for(i in 1: length(list.files(paste("dir",file, sep=""))))
  {d = as.matrix(read.table(paste(paste("dir",file,"/", sep=""),file, "_", i, sep=""), sep=",")) 
  stats = summary_stats(data=d[2:length(d[,1]),], tree=tree,...) # from 2: because line 1 is the one with sigma, alpha and -1
  write(c(d[1,1], d[1,2], stats), file=file, append=TRUE, sep=",")
  }} 
#param_stats_sstat_only(file="test_acc3", tree=phy_list$ACCIPITRIFORMES, ntraits=8)


param_stats_sstat_only_parr = function(file,tree, ncores=30,...)
{#dir.create(paste("/home/angela/Downloads/treecomp-master/R/","dir",file,"/", sep="")) #exists
  length(list.files(paste("dir",file, sep=""))) -> reps
  rez<- mclapply(1:reps, function(i){
    d = as.matrix(read.table(paste(paste("dir",file,"/", sep=""),file, "_", i, sep=""), sep=",")) 
    stats = summary_stats(data=d[2:length(d[,1]),], tree=tree, ...) ## from 2: because line 1 is the one with sigma, alpha and -1
    write(c(d[1,1], d[1,2], stats), file=paste(paste("dir",file,"/", sep=""), "PS",file,"_", i , sep=""), append=F, sep=",") #save it with PS as the simle file might be the sim file !!
  }, mc.cores=ncores, mc.preschedule = F)
  
  # add a for that reads each of those lines, appends them to final_txt, and deletes the file with 1 line
  for(i in 1: reps)
  {  as.numeric(read.table (paste(paste("dir",file,"/", sep=""), "PS",file,"_", i , sep=""), sep=",")) ->l
    write(l, file=file, append=TRUE, sep=",")
    rm(l)
    file.remove(paste(paste("dir",file,"/", sep=""), "PS",file,"_", i , sep=""))
  }
}  
#param_stats_sstat_only_parr(file="test_acc3", tree=phy_list$ACCIPITRIFORMES, ntraits=8, ncores=20)
##-------------------------------------------------------------------------------------------------------------------------------------------------------


# LRT WITH OUTPUT ACCEPTED
lrt_acc = function(file, posteriorSize=500, use_K=FALSE, max_sigma=8, max_a=8, file_name, print_file=T, ...)
{
  # Read simulations into R
  #x   = read.csv(file)
  
  x = read.table(file, sep=",")
  sig = x[,1]
  atry = x[,2]
  sstat = x[, -(1:2)]
  
  # Get summary stats for the true data, and distance to sims
  tstat = summary_stats(use_K=use_K, ...)
  diff = sdiff(sstat, tstat) 
  
  if(print_file==T)
    {# save param stats from closest to end
    H1_save     = order(diff)
    Usig_save   = sig[H1_save]
    Uatry_save = atry[H1_save]
    Usstat_save = sstat[H1_save,]
    H1_save = matrix(ncol=2, nrow=length(Usig_save))
    H1_save[,1] = Usig_save
    H1_save[,2] = Uatry_save
    write.table(cbind(H1_save, as.matrix(Usstat_save)) , file=file_name, row.names = F, col.names=F)
   }
  
  
  # Get simulations from nth closest to closest 
  H1_post     = order(diff)[1:posteriorSize]
  
  Usig        = sig[H1_post]
  Uatry       = atry[H1_post]
  Usstat = sstat[H1_post,]
  
  H1_post = matrix(ncol=2, nrow=length(Usig))
  H1_post[,1] = Usig
  H1_post[,2] = Uatry
  
  #write.table(cbind(H1_post, as.matrix(Usstat)) , file=paste("accepted",posteriorSize,file, sep=""), row.names = F, col.names=F)
  
  k   = kde(H1_post, xmin=c(0, 0), xmax=c(max_sigma,max_a))
  k0  = kde(H1_post, xmin=c(0, 0), xmax=c(max_sigma,0))       # sigma to max, a to 0.
  
  # Use kernel smoothing to estimate likelihood maxima with and without competition.
  k_max_index     = which(k$estimate == max(k$estimate), arr.ind = TRUE)
  H1_lik          = k$estimate[k_max_index[1], k_max_index[2]]
  H1_est          = c(unlist(k$eval.points)[k_max_index[1]], unlist(k$eval.points)[length(k$estimate[,1]) + k_max_index[2]])
  
  k0_max_index    = which(k0$estimate == max(k0$estimate), arr.ind = TRUE)
  H0_lik          = k0$estimate[k0_max_index[1], k0_max_index[2]]
  H0_est          = c(unlist(k0$eval.points)[k0_max_index[1]], unlist(k0$eval.points)[length(k0$estimate[,1]) + k0_max_index[2]])
  
  LRT             = -2 * log( H0_lik / H1_lik )
  
  data.frame(H0_est, H0_lik, H1_est, H1_lik, LRT)
}
#lrt_acc(tree=pcphy_list$Anserinae, data=data_use, file="1MIL_Anserinae_0.001.txt", posteriorSize=100,max_sigma=8, max_a=8, ntraits=8, use_K=F, print_file=T,file_name="test_allsortedAnserinae")
#lrt(tree=pcphy_list$Anserinae, data=data_use, file="1MIL_Anserinae_0.001.txt", posteriorSize=100,max_sigma=8, max_a=8, ntraits=8, use_K=F)
