






opt_loglik = function(pars1, pars1f = c(function(t, pars) {pars[1]},
                                        function(t, pars) {pars[2]},
                                        function(t, pars) {pars[3]},
                                        function(t, pars) {pars[4]}),
                      pars2 = c(1, 0, 2, 0, "lsoda", 0, 0), 
                      brts, missnumspec = 0){
  brts = sort(abs(brts))
  abstol = 1e-16
  reltol = 1e-10
  methode = as.character(pars2[5])
  cond = as.numeric(pars2[1])
  btorph = as.numeric(pars2[2])==0
  soc = as.numeric(pars2[3])
  m = missnumspec
  probs = c(1, 1, 0, 0)
  lbrts = length(brts)
  
  
  if(cond == 0){
    cond_fun = function(y){
      soc * (log(y[lbrts + 1, 2]) + log(1 - y[lbrts + 1, 3]))
    }
  } else{
    cond_fun = function(y){
      soc * log(y[(lbrts + 1), 2])
    }
  }
  
  if (cond == 2) {
    if (as.numeric(pars2[6]) == 0 & as.numeric(pars2[7]) == 0) {
      n_l_fun = function(S) return(S + m)
      n_u_fun = function(S) return(S + m)
    } else{
      n_l_fun = function(S) return(as.numeric(pars2[6]))
      n_u_fun = function(S) return(as.numeric(pars2[7]))
    }
    
    if (m > 0) {
      if (soc == 1) {
        y2 = as.numeric(c(1 - y[2:(lbrts + 1), 2]))
      }
      if (soc == 2) {
        y2 = as.numeric(c(1 - y[2:(lbrts + 1), 2],
                          1 - y[lbrts + 1, 2]))
      }
      x = rep(0, m + 1)
      x[1] = 1
      for (j in 1:S) {
        x = conv(x, (1:(m + 1)) * (y2[j]^(0:m)))[1:(m + 1)]
      }
      
      fun = function (pars1)
      {
        pars1 = c(pars1f, pars1)
        b = pars1[[1]](brts, as.numeric(pars1[5:length(pars1)]))
        S = lbrts + (soc - 1)
        y = deSolve::ode(probs, c(0, brts), PBD:::pbd_loglik_rhs, c(pars1),
                         rtol = reltol,
                         atol = abstol, method = methode)
        if (dim(y)[1] < lbrts + 1) {
          return(-Inf)
        }
        
        n_u = n_u_fun(S)
        n_l = n_l_fun(S)
        if (as.numeric(pars2[7]) == Inf) {
          n_u = n_l - 1
          n_l = soc
        }
        if (n_u < n_l) {
          logcond = 1
        } else {
          one = 1
          if (n_l == 1 & n_u == 1) {
            one = 0
          }
          logcond = ((soc - 1) * log((n_l:n_u) - one) +
                       soc * log(y[(lbrts + 1), 2]) +
                       ((n_l:n_u) - soc) * log(1 - y[(lbrts + 1), 2]))
        }
        if (as.numeric(pars2[7]) == Inf) {
          logcond = 1 - logcond
        }
        
        loglik = btorph * lgamma(S) +
          cond_fun(y) +
          sum(log(b) + log(y[2:lbrts, 2]) + log(1 - y[2:lbrts, 3])) -
          logcond + 
          lgamma(S + 1) + lgamma(m + 1) -
          lgamma(S + m + 1) + log(x[m + 1])
        
        
        return(as.numeric(loglik))
      }
      
    } else{ # m==0
      fun = function (pars1)
      {
        pars1 = c(pars1f, pars1)
        b = pars1[[1]](brts, as.numeric(pars1[5:length(pars1)]))
        S = lbrts + (soc - 1)
        y = deSolve::ode(probs, c(0, brts), PBD:::pbd_loglik_rhs, c(pars1),
                         rtol = reltol,
                         atol = abstol, method = methode)
        if (dim(y)[1] < lbrts + 1) {
          return(-Inf)
        }
        
        n_u = n_u_fun(S)
        n_l = n_l_fun(S)
        if (as.numeric(pars2[7]) == Inf) {
          n_u = n_l - 1
          n_l = soc
        }
        if (n_u < n_l) {
          logcond = 1
        } else {
          one = 1
          if (n_l == 1 & n_u == 1) {
            one = 0
          }
          logcond = ((soc - 1) * log((n_l:n_u) - one) +
                       soc * log(y[(lbrts + 1), 2]) +
                       ((n_l:n_u) - soc) * log(1 - y[(lbrts + 1), 2]))
        }
        if (as.numeric(pars2[7]) == Inf) {
          logcond = 1 - logcond
        }
        loglik = btorph * lgamma(S) +
          cond_fun(y) +
          sum(log(b) + log(y[2:lbrts, 2]) + log(1 - y[2:lbrts, 3])) -
          logcond
        
        
        return(as.numeric(loglik))
      }
    }
  } else if (m > 0) {
    if (soc == 1) {
      y2 = as.numeric(c(1 - y[2:(lbrts + 1), 2]))
    }
    if (soc == 2) {
      y2 = as.numeric(c(1 - y[2:(lbrts + 1), 2],
                        1 - y[lbrts + 1, 2]))
    }
    x = rep(0, m + 1)
    x[1] = 1
    for (j in 1:S) {
      x = conv(x, (1:(m + 1)) * (y2[j]^(0:m)))[1:(m + 1)]
    }
    
    fun = function (pars1)
    {
      pars1 = c(pars1f, pars1)
      b = pars1[[1]](brts, as.numeric(pars1[5:length(pars1)]))
      S = lbrts + (soc - 1)
      y = deSolve::ode(probs, c(0, brts), PBD:::pbd_loglik_rhs, c(pars1),
                       rtol = reltol,
                       atol = abstol, method = methode)
      if (dim(y)[1] < lbrts + 1) {
        return(-Inf)
      }
      loglik = btorph * lgamma(S) +
        cond_fun(y) +
        sum(log(b) + log(y[2:lbrts, 2]) + log(1 - y[2:lbrts, 3])) +
        lgamma(S + 1) + lgamma(m + 1) -
        lgamma(S + m + 1) + log(x[m + 1])
      
      
      return(as.numeric(loglik))
    }
    
  } else{
    fun = function (pars1)
    {
      pars1 = c(pars1f, pars1)
      b = pars1[[1]](brts, as.numeric(pars1[5:length(pars1)]))
      S = lbrts + (soc - 1)
      y = deSolve::ode(probs, c(0, brts), PBD:::pbd_loglik_rhs, c(pars1),
                       rtol = reltol,
                       atol = abstol, method = methode)
      if (dim(y)[1] < lbrts + 1) {
        return(-Inf)
      }
      loglik = btorph * lgamma(S) +
        cond_fun(y) +
        sum(log(b) + log(y[2:lbrts, 2]) + log(1 - y[2:lbrts, 3]))
      
      
      return(as.numeric(loglik))
    }
    
  }
  
  
  return(fun)
}






start.output = function(outputName, logLik, prior, posterior, parameters){
  output = file(paste(outputName, "_PBD_Bayes.txt", sep = ""), "w")
  
  b = parameters[1]
  mu1 = parameters[2]
  la1 = parameters[3]
  mu2 = parameters[4]
  
  cat("logLik", "prior", "posterior", "b", "mu1", "la1", "mu2", "accepted", "proposed.par", "converged", file = output, sep = "\t")
  cat("\n", file = output)
  cat(logLik, prior, posterior, b, mu1, la1, mu2, rep("NA", 3), file = output, sep = "\t")
  cat("\n", file = output)
  
  return(output)
}

make_generate_proposal = function(step, upp = NULL){
  if(!is.null(upp)){
    if(is.null(names(upp))){
      upp = structure(upp, names = c("b", "mu1", "la1", "mu2"))
    }
  }
  if(class(step) %in% c("numeric", "list", "function")){
    UseMethod(generic = "make_generate_proposal", object = step)
  } else{
    stop("'step' must be one of the three: a vector with the standard deviation of the step; 
         a function; or a list of 4 functions.")
  }
}
make_generate_proposal.numeric = function(step, upp = NULL){
  if(is.null(upp)){
    if(length(step) == 1){
      fun = function(var, par){
        x = as.numeric(var[par])
        n = length(par)
        r = rnorm(n, mean = 0, sd = step)
        new.var = x + ifelse(r < -x, -r, r)
        return(new.var)
      }
    } else{
      if(is.null(names(step))){
        names(step) = c("b", "mu1", "la1", "mu2")
      }
      fun = function(var, par){
        x = as.numeric(var[par])
        s = as.numeric(step[par])
        r = sapply(s, FUN = function(yyy) rnorm(1, mean = 0, sd = yyy))
        new.var = x + ifelse(r < -x, -r, r)
        return(new.var)
      }
    }
  } else{
    zeros = rep(0, 4)
    if(length(step) == 1){
      fun = function(var, par){
        while(1){
          x = as.numeric(var[par])
          n = length(par)
          r = rnorm(n, mean = 0, sd = step)
          new.var = x + ifelse(r < -x, -r, r)
          if(all(new.var > zeros & new.var <  upp)){
            return(new.var)
          }
        }
      }
    } else{
      if(is.null(names(step))){
        names(step) = c("b", "mu1", "la1", "mu2")
      }
      fun = function(var, par){
        while(1){
          x = as.numeric(var[par])
          s = as.numeric(step[par])
          r = sapply(s, FUN = function(yyy) rnorm(1, mean = 0, sd = yyy))
          new.var = x + ifelse(r < -x, -r, r)
          if(all(new.var > zeros & new.var <  upp)){
            return(new.var)
          }
        }
      }
    }
  }
  return(fun)
} 
make_generate_proposal.function = function(step, ...){
  # eg, step = function(x) x+runif(length(x), -1, 1)
  fun = function(var, par){
    x = as.numeric(var[par])
    new.var = step(x)
    return(new.var)
  }
  return(fun)
}
make_generate_proposal.list = function(step, ...){
  # eg, step = list(function(x) x+runif(1, -1, 1),
  #                 function(x) x+runif(1, -1, 1),
  #                 function(x) x+runif(1, -1, 1),
  #                 function(x) x+rexp(1, 1))
  if(is.null(names(step))){
    names(step) = c("b", "mu1", "la1", "mu2")
  }
  fun = function(var, par){
    x = as.numeric(var[par])
    new.var = mapply(function(p, xxx) step[[p]](xxx), p = par, xxx = x)
    return(new.var)
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

make_sampler = function(sampler){
  if(is.null(sampler)){
    fun = function(){
      out = sample(c("b", "mu1", "la1", "mu2"), size = rbinom(1, 3, .15)+1)
      
      return(out)
    }
  } else if(is.numeric(sampler) | is.integer(sampler)){
    if(sampler > 4 | sampler < 1){
      stop("'sampler' must be between 1 and 4")
    } else{
      sampler = as.integer(sampler)
      fun = function(){
        out = sample(c("b", "mu1", "la1", "mu2"), size = sampler)
        
        return(out)
      }
    }
  } else if(class(sampler) == "function"){
    fun = function(){
      out = sampler(c("b", "mu1", "la1", "mu2"))
      
      return(out)
    }
  } else{
    stop("'sampler' must be one of the three: NULL; an integer between 1 and 4; or a function.")
  }
  
  return(fun)
}




make_ratio_proposal = function(upp){
  UseMethod("make_ratio_proposal", upp)
}
make_ratio_proposal.null = function(upp){
  fun = function(old, new, id){
    old = old[id]
    diff = abs(old - new)
    old2new = max(1, sum(old < diff) * 2)
    new2old = max(1, sum(new < diff) * 2)
    return(old2new / new2old)
  }
  return(fun)
} 
make_ratio_proposal.numeric = function(upp){
  if(is.null(names(upp))){
    upp = structure(upp, names = c("b", "mu1", "la1", "mu2"))
  }
  fun = function(old, new, id){
    old = old[id]
    lim = upp[id]
    diff = abs(old - new)
    old2new = max(1, sum(old < diff | diff > (lim - old)) * 2)
    new2old = max(1, sum(new < diff | diff > (lim - new)) * 2)
    return(old2new / new2old)
  }
  return(fun)
}



