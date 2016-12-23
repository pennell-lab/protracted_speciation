
opt_loglik = function(pars1, pars1f = c(function(t, pars) {pars[1]},
                                        function(t, pars) {pars[2]},
                                        function(t, pars) {pars[3]},
                                        function(t, pars) {pars[4]}),
                      pars2 = c(1, 1, 2, 1, "lsoda", 0, 0), 
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
        y = deSolve::ode(probs, c(0, brts), PBD::pbd_loglik_rhs, c(pars1),
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
      
    } else{ # m<=0
      fun = function (pars1)
      {
        pars1 = c(pars1f, pars1)
        b = pars1[[1]](brts, as.numeric(pars1[5:length(pars1)]))
        S = lbrts + (soc - 1)
        y = deSolve::ode(probs, c(0, brts), PBD::pbd_loglik_rhs, c(pars1),
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
      y = deSolve::ode(probs, c(0, brts), PBD::pbd_loglik_rhs, c(pars1),
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
      y = deSolve::ode(probs, c(0, brts), PBD::pbd_loglik_rhs, c(pars1),
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
