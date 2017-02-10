




library(PBD)
library(TreeSim)
library(doParallel)

#load("analysis/data/pbd_sim_best_scenario.RData")
load("pbd_sim_best_scenario.RData")

parameters = parameters[-3, ]

rep = 5e+4
chains = 5
files = paste0("Bayes_test_par", c(1,2,4), "_phy")
#source("analysis/R/pbd_Bayes.R")
source("pbd_Bayes.R")
for(i in 1:3){
  prior_b = function(b){
    #dexp(b, rate = 1/parameters[i, 1], log = TRUE)
    dlnorm(b, meanlog = log(parameters[i, 1]), sdlog = 0.1, log = TRUE)
  }
  prior_mu1 = function(mu1){
    #dexp(mu1, rate = 1/parameters[i, 4], log = TRUE)
    dlnorm(mu1, meanlog = log(parameters[i, 4]), sdlog = 0.1, log = TRUE)
  }
  prior_la1 = function(la1){
    dlnorm(la1, meanlog = log(parameters[i, 2]), sdlog = 0.1, log = TRUE)
  }
  prior_mu2 = function(mu2){
    #dexp(mu2, rate = 1/parameters[i, 5], log = TRUE)
    dlnorm(mu2, meanlog = log(parameters[i, 5]), sdlog = 0.1, log = TRUE)
  }
  
  
  
  filog = seq(from = (((i-1)*10)+1), to = (((i-1)*10))+10)
  #for(k in filog){
  mclapply(filog, FUN = function(k){
    foo = Sys.Date()
    cat(paste0(foo, ": start phylo", k, "\t"))
    for(w in 1:chains){
      pbd_Bayes(brts = branches[[k]],
                initparsopt = parameters[i, c(1,4,2,5)],
                prior_b = prior_b, prior_mu1 = prior_mu1,
                prior_la1 = prior_la1, prior_mu2 = prior_mu2,
                step = 0.7, rep = rep, file = paste0(files[i], k, ".", w))
    }
    time = Sys.Date() - foo
    cat(paste0("finish:", time, "\n"))
  }, mc.cores = 2)
  cat(paste0("DONE: parameters ", i, "\n\n\n"))
}







## Maximum likelihood
parameters = matrix(data = c(c(0.3,1,0,0.1,0.1),
                             c(0.3,1,0,0.2,0.1),
                             #c(0.5,1,0,0.1,0.7),
                             c(0.3,0.5,0,0.1,0.1)),
                    ncol = 5, byrow = TRUE)
results = mclapply(1:3, FUN = function(x){
  range = (1:10)+(10*(x-1))
  ml = lapply(branches[range], FUN = pbd_ML, 
              initparsopt = parameters[x, c(1, 4, 2, 5)])
}, mc.cores = 3)



save(results, file = "pbd_best_ML.RData")









