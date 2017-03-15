

library(ape)
library(PBD)
library(doParallel)
library(stats4)


# load(file = "analysis/data/Botero14.RData")
load(file = "Botero14.RData")

# load reduced set of Jetz trees (only 10)
#jetz = read.tree("analysis/data/Jetz_etal/EricsonStage2_0001_1000/AllBirdsEricson1_reduced.tre")
jetz = read.tree("AllBirdsEricson1_reduced.tre")
# check which sp are present in the table and drops the ones that are not
# miss = lapply(jetz, function(x) which(!x$tip.label %in% birds$SPECIES))
# jj = mapply(function(p, m) drop.tip(p, m), p = jetz, m = miss, SIMPLIFY = FALSE)
# sapply(jj, Ntip)




# run the analysis
cores = 5
rep = 5e+4
chains = 10
files = paste0("Jetz_Botero_first_try_phy_", 1:length(jetz))
#source("analysis/R/pbd_Bayes.R")
source("pbd_Bayes.R")
prior_mu1 = function(mu1){
  dexp(mu1, rate = 1/0.5, log = TRUE)
}
prior_la1 = function(la1){
  dlnorm(la1, meanlog = 0, sdlog = 1, log = TRUE)
}
prior_mu2 = function(mu2){
  dexp(mu2, rate = 1/0.5, log = TRUE)
}
# select the species that occur in that given region
sp = lapply(b1names, function(xx) birds$SPECIES[birds$LAT.RANGE %in% xx])
for(i in 1:length(jetz)){
  # get the list of species that are NOT in the list
  miss = lapply(sp, function(x) which(!jetz[[i]]$tip.label %in% x))
  # drop the tips, keeping only what occurs in the region
  jj = lapply(miss, drop.tip, phy = jetz[[i]])
  # get the branch length for the each region
  branches = lapply(jj, branching.times)
  
  for(k in 1:ncol(prior4b1)){
    # create the prior with the estimated parameters, for the given number of subspecies
    prior_b = function(b){
      #dexp(b, rate = 1/sim.pars[1], log = TRUE)
      dlnorm(b, meanlog = as.numeric(prior4b1[1, k]), sdlog = as.numeric(prior4b1[2, k]), log = TRUE)
    }
    
    
    
    mclapply(1:chains, FUN = function(xxx){
      pbd_Bayes(brts = branches[[k]],
                initparsopt = rep(xxx, 4),
                prior_b = prior_b, prior_mu1 = prior_mu1, 
                prior_la1 = prior_la1, prior_mu2 = prior_mu2,
                step = 0.7, rep = rep, file = paste0(files[i], "_region_", paste(b1names[[k]], collapse = "_"), ".", xxx))
    }, mc.cores = cores)
  }
}


