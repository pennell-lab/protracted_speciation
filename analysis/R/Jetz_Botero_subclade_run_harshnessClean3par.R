
# Subclade clean harshness 3 parameters analysis

library(parallel)
library(geiger)

# LOCAL
# setwd("~/Desktop/gitHub/protracted_sp/analysis/R/")
# input_data = "../data/"
# input_code = ""
# ouput_data = "../data/"

# PARALLEL
input_data = "input/"
input_code = ""
ouput_data = "output/"


source(paste0(input_code, "pbd_Bayes3par.R"))
# source(paste0(input_code, "auxiliary_functions.R")) # it was here for "unlist_branches"

# branches, branches.unity, branchesHarsh, branchesHarsh.unity, branchesHarshClean, branchesHarshClean.unity
load(file = paste0(input_data, "Jetz_Botero14_30min_0.3mixed_0multilat_branches.RData"))
# these_branches = unlist_branches(branchesHarsh.unity) # no longer necessary......
# data.frame with info for harshness
propMixedHarsh = read.csv(file = paste0(input_data, "Jetz_Botero14_propMixedHarsh.csv"))
# exponential priors for b1
load(file = paste0(input_data, "Jetz_Botero14_subclades_priorb1expHarshClean.RData"))
#curve(priorb1expHarsh[[7]][[19]](x), 0, 10)



# run the analysis
steps = 5e+6
chains = 2
limits = c(10, 10, 20, 20) # upper bound for parameters; order = "b", "mu1", "la1", "mu2"
rates = structure(c(0.1,0.4,0.1), names = c("la1", "mu1", "mu2")) 
prior_mu1 = function(mu1){
  dexp(mu1, rate = rates["mu1"], log = TRUE)
}
prior_la1 = function(la1){
  dexp(la1, rate = rates["la1"], log = TRUE)
  #dlnorm(la1, meanlog = 1, sdlog = 1, log = TRUE)
}
prior_mu2 = function(mu2){
  dexp(mu2, rate = rates["mu2"], log = TRUE)
}


lapply(unique(propMixedHarsh$phy), function(i){
  # selects only this phylogeny
  this_data = propMixedHarsh[propMixedHarsh$phy == i, ]
  # creates file names
  files = paste0(ouput_data, "Jetz_Botero_subclade_phy_", i,
                 "_harshnessClean3par_", this_data$which.max,
                 "_est_", this_data$estimate)
  # gets the indices for the "Stable" and "Variable" subclades (ie, excludes "Regular")
  ind = which(this_data$which.max != "Regular")
  ind = ind[-which(ind == 111)]
  mclapply(ind, FUN = function(xxx){
    for(chn in 1:chains){
      pbd_Bayes3par(brts = branchesHarshClean.unity[[i]][[xxx]],
                    initparsopt = rep(chn, 4),
                    prior_b = priorb1expHarshClean[[i]][[xxx]],
                    prior_mu1 = prior_mu1, 
                    prior_la1 = prior_la1,
                    prior_mu2 = prior_mu2,
                    step = c(0.5, 0.5, 1.5, 0.5), 
                    upper = limits, 
                    rep = steps, 
                    file = paste(files[xxx], chn, sep = "_chn_"))
    }
  }, mc.cores = 10)
})
