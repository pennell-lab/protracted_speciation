

library(ape)
library(PBD)
library(doParallel)

#output = "~/Desktop/gitHub/protracted_sp/analysis/R/"
output = ""
out_sim = "simulations/"
out_post = "posterior/"

source("pbd_Bayes.R")
load(file = paste0(output, out_sim, "simulation_Cook2.RData"))

rep = 5e+4
cores = 16
burnin = ceiling(rep*0.3)
samp = 100
chains = 2
convergence = 1.1
sel = seq(burnin, rep, by = samp) # gets the rows that will be sampled
limits = c(10, 10, 20, 20) # upper bound for parameters; order = "b", "mu1", "la1", "mu2"
rates = structure(c(0.1,0.1,0.4,0.1), names = c("b", "la1", "mu1", "mu2")) # OLD = structure(c(0.1,0.1,0.4,0.1), names = c("b", "la1", "mu1", "mu2"))

#quant_b = function(b){plnorm(b, meanlog = 0, sdlog = 1)}
quant_b = function(b){pexp(b, rate = rates["b"])}
#quant_la1 = function(la1){plnorm(la1, meanlog = 1, sdlog = 1)}
quant_la1 = function(la1){pexp(la1, rate = rates["la1"])}
quant_mu1 = function(mu1){pexp(mu1, rate = rates["mu1"])}
quant_mu2 = function(mu2){pexp(mu2, rate = rates["mu2"])}
get_true_quantile = function(values){
  quantiles = c(quant_b, quant_la1, quant_mu1, quant_mu2)
  out = mapply(function(f,x) f(x), x = values, f = quantiles)
  out
}
get_posterior_quantile = function(values, posterior){
  nr = nrow(posterior)
  out = sapply(1:4, function(i){sum(posterior[ , i] > values[i])})/nr
  out
}
prior_b = function(b){
  dexp(b, rate = rates["b"], log = TRUE)
  #dlnorm(b, meanlog = 0, sdlog = 1, log = TRUE)
}
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

bbb = which(apply(true.par, 1, function(x) !any(x[c(1,3,2,4)] > limits)))
files = apply(true.par[bbb, ], 1, FUN = function(x) paste0(output, out_post, "Bayes_simulation_study_", paste0(round(x,3), collapse = "_")))
mclapply(1:length(bbb), FUN = function(j){
  for(i in 1:chains){
    pbd_Bayes(brts = branches[bbb[j], ],
              initparsopt = rep(i, 4),
              prior_b = prior_b, prior_mu1 = prior_mu1, 
              prior_la1 = prior_la1, prior_mu2 = prior_mu2,
              step = c(0.5, 0.5, 1.5, 1.5), 
              upper = limits, 
              rep = rep, 
              file = paste(files[[j]], i, sep = "."))
  }
}, mc.cores = cores)

# combine them into one single data.frame to store all the info
final = mclapply(1:length(files), FUN = function(j){
  res = lapply(1:chain, function(i) {
      read.csv(file = paste0(paste(files[[j]], i, sep = "."), "_PBD_Bayes.txt"), sep = "\t", row.names = NULL, as.is = TRUE)
    })
  my.mcmc = lapply(lapply(res, "[", 4:7), mcmc, start = burnin, thin = samp)
  my.mcmc = mcmc.list(my.mcmc)
  test = gelman.diag(my.mcmc)$mpsrf
  rm(my.mcmc)
  
  sam = lapply(X = res, FUN = function(x) return(x[sel, c("b", "la1", "mu1", "mu2")]))
  samples = do.call(rbind, sam)
  parameters = true.par[bbb[j], ]
  out = matrix(c(get_posterior_quantile(parameters, samples), get_true_quantile(parameters)), ncol = 2)
  rm(samples)
  return(list(mat = out, conv = (test <= convergence)))
}, mc.cores = cores)
save(final, "final.RData")

convergence = sapply(final, "[[", 2)
cook = do.call(rbind, lapply(final, "[[", 1))
ml = mclapply(bbb, FUN = function(j){pbd_ML(brts = branches[j, ], initparsopt = c(1,0.1,1,0.1), exteq = 0, btorph = 0, soc = 2, verbose = FALSE)}, mc.cores = cores)
ml = do.call(rbind, ml)
save(convergence, cook, ml, "final.RData")


cat("DONE!!!")
