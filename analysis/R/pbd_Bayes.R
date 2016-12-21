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
                     rep = 1e+5, step = 1,
                     ...) {
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
  if(length(initparsopt) == 3){
    initparsopt[4] = initparsopt[2]
  }
  logLik1 = pbd_loglik(pars1 = initparsopt, brts = brts)
  prior1 = prior_b(initparsopt[1]) + prior_mu1(initparsopt[2]) + 
    prior_la1(initparsopt[3]) + prior_mu2(initparsopt[4])
  post1 = logLik1 + prior1
  out[1, 1:7] = c(logLik1, prior1, post1, initparsopt)
  
  
  new.pars = setNames(initparsopt, c("b", "mu1", "la1", "mu2"))
  for(i in 1:rep){
    par = sample(c("b", "mu1", "la1", "mu2"), size = 1)
    proposal = generate_proposal(step = step, var = out[i, par])
    new.pars[par] = proposal
    
    new.logLik = pbd_loglik(pars1 = new.pars, brts = brts)
    new.prior = prior_b(as.numeric(new.pars[1])) + 
      prior_mu1(as.numeric(new.pars[2])) + 
      prior_la1(as.numeric(new.pars[3])) + 
      prior_mu2(as.numeric(new.pars[4]))
    new.post = new.logLik + new.prior
    
    diff.logLik = out[i, "logLik"] - new.logLik
    diff.prior = out[i, "prior"] - new.prior
    diff = diff.logLik + diff.prior
    
    accept = exp(diff) > runif(1)
    out[i+1, "accepted"] = accept
    out[i+1, "proposed.par"] = par
    if(is.na(accept) | is.nan(accept)){
      out[i+1, "converged"] = FALSE
      out[i+1, 1:7] = out[i, 1:7]
      new.pars = out[i, 4:7]
    } else{
      out[i+1, "converged"] = TRUE
      if(accept){
        out[i+1, 1:7] = c(new.logLik, new.prior, new.post, new.pars)
      } else{
        out[i+1, 1:7] = out[i, 1:7]
        new.pars = out[i, 4:7]
      }
    }
  }
  
  return(out)
}

generate_proposal = function(step, var){
  new.var = runif(1, min = max(var - (step/2), 0), max = var + (step/2) )
  return(new.var)
}