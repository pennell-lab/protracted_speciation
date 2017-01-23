




library(PBD)
library(TreeSim)
library(doParallel)

#load("analysis/data/pbd_sim.RData")
load("pbd_sim.RData")
sim.pars = apply(parameters, MARGIN = 2, mean)


# library(msm) for "dtnorm" - truncated normal
prior_b = function(b){
  dexp(b, rate = 1/sim.pars[1], log = TRUE)
}
prior_mu1 = function(mu1){
  dexp(mu1, rate = 1/sim.pars[4], log = TRUE)
}
prior_la1 = function(la1){
  dlnorm(la1, meanlog = log(sim.pars[2]), sdlog = log(15), log = TRUE)
}
prior_mu2 = function(mu2){
  dexp(mu2, rate = 1/sim.pars[5], log = TRUE)
}


rep = 5e+4
chains = 5
files = paste0("Bayes_test_phy", 1:6)
times = vector(mode = "numeric")
#source("analysis/R/pbd_Bayes.R")
source("pbd_Bayes.R")
for(i in 1:6){
  for(k in 1:chains){
    foo = Sys.time()
    cat(paste0(foo, "\n"))
    pbd_Bayes(brts = branches[[i]],
              initparsopt = rep(k, 4),
              prior_b = prior_b, prior_mu1 = prior_mu1, 
              prior_la1 = prior_la1, prior_mu2 = prior_mu2,
              step = 0.7, rep = rep, file = paste0(files[i], ".", k))
    times[i] = Sys.time() - foo
    cat(paste0(times[i], "\n"))
  }
}