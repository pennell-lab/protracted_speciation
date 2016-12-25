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
                     step = 1, # step size for new proposed values; can be a function; can also be passed different values/functions for each variable (i.e., b, mu1, la1, mu2) 
                     rep = 1e+5,
                     ...) {
  
  if(length(initparsopt) == 3){
    initparsopt[4] = initparsopt[2]
  }
  
  # making functions
  generate_proposal = make_generate_proposal(step)
  prior_logLik = make_prior_logLik(prior_b, prior_mu1, prior_la1, prior_mu2)
  logLik_fun = opt_loglik(brts = brts, ...)
  
  
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
    par = sample(c("b", "mu1", "la1", "mu2"), size = 1)
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

make_generate_proposal = function(step){
  if(is.numeric(step)){
    if(length(step) == 1){
      s = step/2
      fun = function(var, par){
        x = as.numeric(var[par])
        new.var = runif(1, min = max(x - s, 0), max = x + s )
        return(new.var)
      }
    } else{
      step2 = setNames(step/2, c("b", "mu1", "la1", "mu2"))
      fun = function(var, par){
        x = as.numeric(var[par])
        s = as.numeric(step2[par])
        new.var = runif(1, min = max(x - s, 0), max = x + s )
        return(new.var)
      }
    }
  } else if(class(step) == "function"){
    fun = function(var, par){
      x = as.numeric(var[par])
      new.var = step(x)
      return(new.var)
    }
  } else if(class(step) == "list"){
    names(step) = c("b", "mu1", "la1", "mu2")
    fun = function(var, par){
      x = as.numeric(var[par])
      new.var = step[[par]](x)
      return(new.var)
    }
  } else{
    stop("'step' must be one of the three: a vector with the size of the step; a function; or a list of functions.")
  }
  
  return(fun)
}


make_prior_logLik = function(prior_b, prior_mu1, prior_la1, prior_mu2){
  logLik = function(vector){
    b = prior_b(as.numeric(vector[1]))
    mu1 = prior_mu1(as.numeric(vector[2]))
    la1 = prior_la1(as.numeric(vector[3]))
    mu2 = prior_mu2(as.numeric(vector[4]))
    
    return(b + mu1 + la1 + mu2)
  }
  
  return(logLik)
}




source("analysis/R/opt_logLik.R")



