

library(ape)
library(PBD)
library(doParallel)

#output = "~/Desktop/gitHub/protracted_sp/analysis/R/"
output = ""
out_sim = "simulations/"
out_post = "posterior/"


rep = 5e+4
cores = 15
burnin = ceiling(rep*0.3)
samp = 100
chains = 2
convergence = 1.1
sel = seq(burnin, rep, by = samp) # gets the rows that will be sampled
limits = c(20, 10, 20, 20) # upper bound for parameters; order = "b", "mu1", "la1", "mu2"


# source("opt_pbd_sim_cpp.R")
# tries = 1000
# ntaxa = 256 # 2^(5:10) = c(32, 64, 128, 256, 512, 1024)
# nphy = 500
# samp_priors = function(){
#   while(1){
#     #b = rlnorm(1, meanlog = 0, sdlog = 1)
#     b = rexp(1, rate = 0.1)
#     mu1 = rexp(1, rate = 0.4)
#     #la1 = rlnorm(1, meanlog = 1, sdlog = 1)
#     la1 = rexp(1, rate = 0.1)
#     #mu2 = rexp(1, rate = 0.1)
#     mu2 = mu1
#     #out = c(b, mu1, la1, mu2)
#     out = c(b, la1, 0, mu1, mu2)
#     if(b < limits[1] & mu1 < limits[2] & la1 < limits[3] & mu2 < limits[4]){
#       return(out)
#     }
#   }
# }
# simulations = mclapply(1:Nphy, FUN = function(k){
#   while(1){
#     par = samp_priors()
#     sim = try(opt_pbd_sim_cpp(pars = par, taxa = ntaxa, ntry = tries), silent = TRUE)
#     if(class(sim) == "list"){
#       this_age = sim$L[1, 3]
#       p = try(get_phylo(sim$L, this_age), silent = TRUE)
#       if(class(p) == "phylo"){
#         ppp = drop.fossil(p)
#         if(Ntip(ppp) == ntaxa){
#           br = branching.times(ppp)
#           return(list(brnchs = br/max(br), true.par = par[-3], phylogeny = p))
#         }
#       }
#     }
#   }
# }, mc.cores = cores)
# save(simulations, file = paste0(output, out_sim, "simulation3.RData"))
# true.par = do.call(rbind, lapply(simulations, "[[", 2))
# colnames(true.par) = c("b", "la1", "mu1", "mu2")
# phylogenies = lapply(simulations, "[[", 3)
# branches = do.call(rbind, lapply(simulations, "[[", 1))
# save(branches, true.par, phylogenies, file = paste0(output, out_sim, "simulation3.RData"))












load(file = paste0(output, out_sim, "simulation3.RData"))
source("auxiliary_BPBD.R")
###########################################
#   Function to fit the protracted BD model
# using Bayesian statistics developed based
# on code from the 'PBD' package
##########################################

pbd_Bayes = function(brts, # branching times
                      initparsopt = c(0.2, 0.1, 1), # initial parameters
                      prior_b,   # logLik function for 'b' - split
                      prior_mu1, # logLik function for 'mu1' - ext good
                      prior_la1, # logLik function for 'la1' - sp completion
                      prior_mu2 = prior_mu1, # logLik function for 'mu2' - ext incipient
                      step = 0.5, # standard deviation of the normal distribution used to generate new proposed values; can be a function; can also be passed different sd/functions for each variable (i.e., b, mu1, la1, mu2)
                      sampler = NULL, # controls how should the parameters be updated; 'NULL' the number of parameters updated is sampled from a binomial distribution; integer between 1:4 updates integer parameters at a time; function to determine how many parameters should be updated each time
                      upper = c(20, 10, 20, 20), # upper limit for the parameters proposal; can be NULL (ie, NO limit) or a vector of 4 in order = "b", "mu1", "la1", "mu2".
                      rep = 1e+5,
                      file = NULL,
                      ...) {
  
  # making functions
  generate_proposal = make_generate_proposal(step, upper)
  prior_logLik = make_prior_logLik(prior_b, prior_mu1, prior_la1, prior_mu2)
  logLik_fun = opt_loglik(brts = brts, ...)
  sampler_fun = function(){
    out = sample(c("b", "mu1", "la1"), size = rbinom(1, 2, .15)+1)
    
    return(out)
  }
  get_ratio_proposal = make_ratio_proposal(upper)
  
  if(class(file) == "character"){
    # initializing the output
    logLik1 = logLik_fun(pars1 = initparsopt)
    prior1 = prior_logLik(initparsopt)
    post1 = logLik1 + prior1
    output = start.output(outputName = file, logLik = logLik1,
                          prior = prior1, posterior = post1,
                          parameters = initparsopt)
    
    new.pars = setNames(initparsopt, c("b", "mu1", "la1", "mu2"))
    old.pars = new.pars
    old.Lik = c(logLik1, prior1, post1)
    for(i in 1:rep){
      par = sampler_fun()
      proposal = generate_proposal(var = new.pars, par = par)
      ratio.proposal = get_ratio_proposal(new.pars, proposal, par)
      if("mu1" %in% par){
        par = c(par, "mu2")
        proposal = structure(c(proposal, proposal[which(par=="mu1")]), names = par)
      }
      new.pars[par] = proposal
      
      new.logLik = logLik_fun(pars1 = new.pars)
      new.prior = prior_logLik(new.pars)
      new.post = new.logLik + new.prior
      
      ratio.logLik = new.logLik - old.Lik[1]
      ratio.prior = new.prior - old.Lik[2]
      ratio = ratio.logLik + ratio.prior + ratio.proposal
      
      accept = exp(ratio) > runif(1)
      if(is.na(accept) | is.nan(accept)){
        converged = FALSE
        cat(old.Lik, old.pars, accept, paste(par, collapse = "_"), converged,
            file = output, sep = "\t")
        new.pars = old.pars
      } else{
        converged = TRUE
        if(accept){
          cat(new.logLik, new.prior, new.post, new.pars,
              accept, paste(par, collapse = "_"), converged,
              file = output, sep = "\t")
          old.Lik = c(new.logLik, new.prior, new.post)
          old.pars = new.pars
        } else{
          cat(old.Lik, old.pars, accept, paste(par, collapse = "_"), converged,
              file = output, sep = "\t")
          new.pars = old.pars
        }
      }
      
      cat("\n", file = output)
    }
    
    
    close(output)
    
    return(invisible())
  }
}
#quant_b = function(b){plnorm(b, meanlog = 0, sdlog = 1)}
quant_b = function(b){pexp(b, rate = 0.1)}
#quant_la1 = function(la1){plnorm(la1, meanlog = 1, sdlog = 1)}
quant_la1 = function(la1){pexp(la1, rate = 0.1)}
quant_mu1 = function(mu1){pexp(mu1, rate = 0.4)}
quant_mu2 = function(mu2){pexp(mu2, rate = 0.1)}
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
  dexp(b, rate = 0.1, log = TRUE)
  #dlnorm(b, meanlog = 0, sdlog = 1, log = TRUE)
}
prior_mu1 = function(mu1){
  dexp(mu1, rate = 0.4, log = TRUE)
}
prior_la1 = function(la1){
  dexp(la1, rate = 0.1, log = TRUE)
  #dlnorm(la1, meanlog = 1, sdlog = 1, log = TRUE)
}
prior_mu2 = function(mu2){
  #dexp(mu2, rate = 0.1, log = TRUE)
  0
}

files = paste0(output, out_post, "Bayes_simulation3par_phy_", 1:nrow(branches))
mclapply(1:nrow(branches), FUN = function(j){
  for(i in 1:chains){
    pbd_Bayes(brts = branches[j, ],
              initparsopt = rep(i, 4),
              prior_b = prior_b, prior_mu1 = prior_mu1,
              prior_la1 = prior_la1, prior_mu2 = prior_mu2,
              step = c(0.5, 0.5, 1.5, 1.5),
              upper = limits,
              rep = rep,
              file = paste(files[[j]], i, sep = "_chain"))
  }
}, mc.cores = cores)

# combine them into one single data.frame to store all the info
final = mclapply(1:length(files), FUN = function(j){
  res = lapply(1:chain, function(i) {
    read.csv(file = paste0(paste(files[[j]], i, sep = "_chain"), "_PBD_Bayes.txt"), sep = "\t", row.names = NULL, as.is = TRUE)
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
save(final, "final3.RData")

convergence = sapply(final, "[[", 2)
cook = do.call(rbind, lapply(final, "[[", 1))
ml = mclapply(bbb, FUN = function(j){
  pbd_ML(brts = branches[j, ], initparsopt = c(1,0.1,1), exteq = 1, btorph = 0, soc = 2, verbose = FALSE)
}, mc.cores = cores)
ml = do.call(rbind, ml)
save(convergence, cook, ml, "final3.RData")

cat("DONE!!!")
