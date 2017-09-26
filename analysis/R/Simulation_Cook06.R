



library(ape)
library(PBD)
library(Rcpp)
library(RcppArmadillo)
library(doParallel)
library(coda)



setwd("~/Desktop/gitHub/protracted_sp/analysis/R/")
output = "~/Desktop/gitHub/protracted_sp/analysis/data/"

source("opt_pbd_sim_cpp.R")


limits = c(20, 10, 20, 20) # upper bound for parameters; order = "b", "mu1", "la1", "mu2"
cores = 4
tries = 100
ntaxa = 256 # 2^(5:10) = c(32, 64, 128, 256, 512, 1024)
Nphy = 1000
samp_priors = function(){
  while(1){
    #b = rlnorm(1, meanlog = 0, sdlog = 1)
    b = rexp(1, rate = 0.2)
    mu1 = rexp(1, rate = 0.5)
    #la1 = rlnorm(1, meanlog = 1, sdlog = 1)
    la1 = rexp(1, rate = 0.2)
    mu2 = rexp(1, rate = 0.2)
    out = c(b, la1, 0, mu1, mu2)
    if(all(b <= limits[1], mu1 <= limits[2], la1 <= limits[3], mu2 <= limits[4])){
      return(out) 
    }
  }
}
sim = function(x){
  out1 = try(opt_pbd_sim_cpp(pars = x, taxa = ntaxa, ntry = tries), silent = TRUE)
  if(class(out1) != "try-error"){
    this_age = out1$L[1, 3]
    phy = try(get_phylo(out1$L, this_age), silent = TRUE)
    return(phy)
  }
  return(out1)
}
simulations = mclapply(1:Nphy, FUN = function(x){
  while(1){
    parameters = samp_priors()
    phy = sim(parameters)
    if(class(phy) == "phylo"){
      phy = drop.fossil(phy)
      if(Ntip(phy) == ntaxa){
        b = branching.times(phy)
        return(list(branches = b/max(b), true.par = parameters[-3]))
      }
    }
  }
}, mc.cores = cores)
branches = do.call(rbind, lapply(simulations, "[[", 1))
true.par = do.call(rbind, lapply(simulations, "[[", 2))
colnames(true.par) = c("b", "la1", "mu1", "mu2")
rm(simulations)
save(branches, true.par, file = paste0(output, "simulation_Cook.RData"))









source("pbd_Bayes.R")
rep = 5e+4
burnin = ceiling(rep*0.3)
samp = 100
chains = 2
convergence = 1.1
sel = seq(burnin, rep, by = samp) # gets the rows that will be sampled

#quant_b = function(b){plnorm(b, meanlog = 0, sdlog = 1)}
quant_b = function(b){pexp(b, rate = 0.2)}
#quant_la1 = function(la1){plnorm(la1, meanlog = 1, sdlog = 1)}
quant_la1 = function(la1){pexp(la1, rate = 0.2)}
quant_mu1 = function(mu1){pexp(mu1, rate = 0.5)}
quant_mu2 = function(mu2){pexp(mu2, rate = 0.2)}
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
  dexp(b, rate = 0.2, log = TRUE)
  #dlnorm(b, meanlog = 0, sdlog = 1, log = TRUE)
}
prior_mu1 = function(mu1){
  dexp(mu1, rate = 0.5, log = TRUE)
}
prior_la1 = function(la1){
  dexp(la1, rate = 0.2, log = TRUE)
  #dlnorm(la1, meanlog = 1, sdlog = 1, log = TRUE)
}
prior_mu2 = function(mu2){
  dexp(mu2, rate = 0.2, log = TRUE)
}

files = apply(true.par, 1, FUN = function(x) paste0(output, "Bayes_simulation_study_", paste0(round(x,3), collapse = "_")))
mclapply(1:nrow(branches), FUN = function(j){
  lapply(1:chains, FUN = function(i){
    pbd_Bayes(brts = branches[j, ],
              initparsopt = rep(i, 4),
              prior_b = prior_b, prior_mu1 = prior_mu1, 
              prior_la1 = prior_la1, prior_mu2 = prior_mu2,
              step = c(1.5, 1.5, 2.5, 2.5), 
              upper = limits,
              rep = rep,
              file = paste(files[[j]], i, sep = "."))
  })
}, mc.cores = cores)

cook = mclapply(1:nrow(branches), FUN = function(j){
  res = lapply(1:chain, function(i) {
    read.csv(file = paste0(paste(files[[j]], i, sep = "."), "_PBD_Bayes.txt"), sep = "\t", row.names = NULL, as.is = TRUE)
  })
  # combine them into one single data.frame to store all the info
  my.mcmc = lapply(lapply(res, "[", 4:7), mcmc, start = burnin, thin = samp)
  my.mcmc = mcmc.list(my.mcmc)
  test = gelman.diag(my.mcmc)$mpsrf
  rm(my.mcmc)
  if(test <= convergence){
    samples = mclapply(X = res, FUN = function(x) return(x[sel, c("b", "la1", "mu1", "mu2")]), mc.cores = cores)
    samples = do.call(rbind, samples)
    parameters = true.par[j, ]
    out = matrix(c(get_posterior_quantile(parameters, samples), get_true_quantile(parameters)), ncol = 2)
    rm(samples)
  } else{
    out = matrix(ncol = 2, nrow = 4)
  }
  return(out)
}, mc.cores = cores)

cook.quantiles = data.frame(variable = rep(c("b", "la1", "mu1", "mu2"), nrow(branches)), posterior = unlist(lapply(cook, function(x)as.numeric(x[,1]))), true = unlist(lapply(cook, function(x)as.numeric(x[,2]))))

write.csv(x = cook.quantiles, file = paste0(output, "simulation_Cook.csv"))



