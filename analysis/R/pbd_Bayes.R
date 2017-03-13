###########################################
#   Function to fit the protracted sp model
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
                     rep = 1e+5,
                     file = NULL,
                     ...) {
  
  if(length(initparsopt) == 3){
    initparsopt[4] = initparsopt[2]
  }
  
  # making functions
  generate_proposal = make_generate_proposal(step)
  prior_logLik = make_prior_logLik(prior_b, prior_mu1, prior_la1, prior_mu2)
  logLik_fun = opt_loglik(brts = brts, ...)
  sampler_fun = make_sampler(sampler)
  
  
  if(is.null(file)){
    # initializing the output
    out = data.frame("logLik" = rep(999, rep+1),
                     "prior" = rep(999, rep+1),
                     "posterior" = rep(999, rep+1),
                     "b" = rep(-999, rep+1),
                     "mu1" = rep(-999, rep+1),
                     "la1" = rep(-999, rep+1),
                     "mu2" = rep(-999, rep+1),
                     "accepted" = logical(rep+1),
                     "proposed.par" = character(rep+1),
                     "converged" = logical(rep+1),
                     stringsAsFactors = FALSE
    )
    logLik1 = logLik_fun(pars1 = initparsopt)
    prior1 = prior_logLik(initparsopt)
    post1 = logLik1 + prior1
    out[1, 1:7] = c(logLik1, prior1, post1, initparsopt)
    
    
    new.pars = setNames(initparsopt, c("b", "mu1", "la1", "mu2"))
    for(i in 1:rep){
      par = sampler_fun()
      proposal = generate_proposal(var = new.pars, par = par)
      new.pars[par] = proposal
      
      new.logLik = logLik_fun(pars1 = new.pars)
      new.prior = prior_logLik(new.pars)
      new.post = new.logLik + new.prior
      
      ratio.logLik = new.logLik - out[i, "logLik"]
      ratio.prior = new.prior - out[i, "prior"]
      ratio = ratio.logLik + ratio.prior
      
      accept = exp(ratio) > runif(1)
      out[i+1, "accepted"] = accept
      out[i+1, "proposed.par"] = par
      if(is.na(accept) | is.nan(accept)){
        out[i+1, "converged"] = FALSE
        out[i+1, 1:7] = out[i, 1:7]
        new.pars = setNames(out[i, 4:7], c("b", "mu1", "la1", "mu2"))
      } else{
        out[i+1, "converged"] = TRUE
        if(accept){
          out[i+1, 1:7] = c(new.logLik, new.prior, new.post, new.pars)
        } else{
          out[i+1, 1:7] = out[i, 1:7]
          new.pars = setNames(out[i, 4:7], c("b", "mu1", "la1", "mu2"))
        }
      }
    }
    
    return(out)
  }
  
  
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
      cat(paste0(i, "\n"))
      par = sampler_fun()
      proposal = generate_proposal(var = new.pars, par = par)
      new.pars[par] = proposal
      
      new.logLik = logLik_fun(pars1 = new.pars)
      new.prior = prior_logLik(new.pars)
      new.post = new.logLik + new.prior
      
      ratio.logLik = new.logLik - old.Lik[1]
      ratio.prior = new.prior - old.Lik[2]
      ratio = ratio.logLik + ratio.prior
      
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
  }
}


source("analysis/R/auxiliary_functions.R")
source("analysis/R/opt_logLik.R")



