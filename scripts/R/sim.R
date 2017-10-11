# Functions to simulate density-dependent phylogenetic data
# magnusclarke@gmail.com
# modified 2016

require('ks')
require('picante')
source('scripts/R/tree.R')
if(.Platform$pkgType == "mac.binary")   dyn.load("scripts/cpp/Rfunc_mac.so")
if(.Platform$pkgType == "source")       dyn.load("scripts/cpp/Rfunc.so")

#--------------------------------------------------------------------------------------#
#--- Convert time-matrices of when species become sympatric/allopatric to vectors. ----#
#--- If matrix is NA then all elements are set to constant number n. ------------------#
#--------------------------------------------------------------------------------------#
vectortime = function(s, t, n)
{
    if(all(is.na(s)))
    {
        s = replicate(n^2, t)
    } else if(class(s)=="matrix") {
        s = as.vector(s)
    }
    s
}
    
#--------------------------------------------------------------------------------------#
#----- Reorder a dataset to match the order of tips in an ape-format tree -------------#
#--------------------------------------------------------------------------------------#
reorder_data = function(t, d, n)
{
    rd = matrix(ncol=n, nrow=length(t$data_order ) )
    for (i in 1:n) 
    {
        rd[,i] = d[seq(i, length(d), by=n)]
        rd[,i] = rd[,i][t$data_order]
    }
    rd
}

#--------------------------------------------------------------------------------------#
#--- Get a dataset simulated under BM + competition, for a given tree. ----------------#
#--- Returns trait values for tips in order corresponding to ape tree tips. -----------#
#--------------------------------------------------------------------------------------#
sim = function(tree, dt=0.001, sigma=1, a=0, ntraits=1, symp=NA, allo=NA, lim=NA, ...)
{
    tree = checktree(tree)

    num_tips = length(tree$data_order)
    splitting_nodes = tree$splitting_nodes - 1  # R counts from 1; c counts from 0.
    times = tree$times
    tval = rep(0, num_tips*ntraits) 

    symp = vectortime(symp, 0, num_tips)
    allo = vectortime(allo, 99e9, num_tips )

    if(is.na(lim))  lim = 9e99

    result = .C ("pathsim", ntip=as.integer(num_tips), dt=as.double(dt), 
                rate = as.double(sigma^2), a=as.double(a), r_intervals=as.double(times), 
                splitters=as.integer(splitting_nodes), tval = as.double(tval),
                ntraits=as.integer(ntraits), symp=as.double(symp), allo=as.double(allo),
                lim=as.numeric(lim)
                )

    result$tval = reorder_data(t=tree, d=result$tval, n=ntraits)
    result$symp = NULL
    result$allo = NULL

    result
}

#--------------------------------------------------------------------------------------#
#--- Generate a matrix of the times at which lineages come into sympatry if -----------#
#--- each new lineage starts interacting after a delay. -------------------------------#
#--------------------------------------------------------------------------------------#
symp_matrix = function(tree, delay=0)
{
    t = checktree(tree)

    ntip = length(t$data_order)
    s = matrix(nrow=ntip, ncol=ntip, 0)
    if(class(tree)=="phylo")
    {
        rownames(s) = colnames(s) = tree$tip.label
    }

    # find ages (time from root) of tip lineages
    age = 1:ntip
    age[1]=age[2]=0
    for(i in 1:ntip)
    {
        age[i+2] = sum(t$times[1:i])
    }

    # delay measured in terms of mean time between speciation events.
    time_between_splits = mean(abs(diff(t$times)))
    delay = delay * time_between_splits

    # apply starttimes and delay to symp-matrix s
    for (i in 1:ntip) 
    {
        for (j in 1:ntip)
        {
            s[i,j] = max(c(age[i], age[j])) + delay
        }
    }
    return(s)
}

#--------------------------------------------------------------------------------------#
#---------- Generate vcv-matrix from simulated trees. ---------------------------------#
#--------------------------------------------------------------------------------------#
as_vcv  = function(reps=1e4, ...)
{
    args = as.list(substitute(list(...)))[-1L]
    x = t( replicate(reps, do.call(sim, args)$tval, simplify='matrix') )
    cov(x)
}

#--------------------------------------------------------------------------------------#
#-------- Get summary statistics for a given tree and dataset -------------------------#
#--------------------------------------------------------------------------------------#
summary_stats = function(tree, data, use_K=FALSE, ...)
{
    difs = as.matrix(dist(data))            # Euclidian distance
    difs[which(difs==0)] = NA               # Ignore matrix diagonal
    gap = apply(difs, 1, min, na.rm=T)  

    data = as.matrix(data)

    ntraits = length(data[1,])

    # Summary statistics: mean and sd of gaps between neighbours. Plus Blomberg's K optionally.
    k = c()
    if(use_K)
    {
        k = 1:ntraits
        for (i in 1:ntraits) 
        {
            k[i] = tryCatch(Kcalc(data[,i], tree, F), error=function(err){return(1)})
        }
        k = mean(k)
    }
    c(mean(gap), sd(gap), k)
}

#--------------------------------------------------------------------------------------#
#-------- Generate file with many simulation parameters and summary statistics --------#
#--------------------------------------------------------------------------------------#
param_stats = function(tree, file='param_stats.out', reps=1e3, max_sigma=8, max_a=8, ...)
{
    sig = runif(reps, 0, max_sigma)
    atry = runif(reps, 0, max_a)

    for(i in 1:reps)
    {
        d = sim(tree=tree, a=atry[i], sigma=sig[i], ...)$tval
        stats = summary_stats(data=d, ...) 
        write(c(sig[i], atry[i], stats), file=file, append=TRUE, sep=",")
    }
}


#--------------------------------------------------------------------------------------#
#-- Get sum of squared distances between summary stats. s1 has many rows, s2 only 1. --#
#--------------------------------------------------------------------------------------#
sdiff = function(s1, s2)
{
    colSums(abs(t(s1) - s2)^2)  
}

#--------------------------------------------------------------------------------------#
#---------- Likelihood ratio: BM versus competition -----------------------------------#
#--------------------------------------------------------------------------------------#
lrt = function(file, posteriorSize=500, use_K=FALSE, max_sigma=8, max_a=8, ...)
{
    # Read simulations into R
    x   = read.csv(file)
    sig = x[,1]
    atry = x[,2]
    sstat = x[, -(1:2)]

    # Get summary stats for the true data, and distance to sims
    tstat = summary_stats(use_K=use_K, ...)
    diff = sdiff(sstat, tstat)

    # Get simulations from nth closest to closest 
    H1_post     = order(diff)[1:posteriorSize]

    Usig        = sig[H1_post]
    Uatry       = atry[H1_post]
    H1_post     = matrix(ncol=2, nrow=length(Usig))
    H1_post[,1] = Usig
    H1_post[,2] = Uatry

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

