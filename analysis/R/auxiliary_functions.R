





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

make_generate_proposal = function(step){
  if(is.numeric(step)){
    if(length(step) == 1){
      fun = function(var, par){
        x = as.numeric(var[par])
        r = rnorm(1, mean = 0, sd = step)
        new.var = x + ifelse(r < -x, -r, r)
        return(new.var)
      }
    } else{
      step2 = setNames(step, c("b", "mu1", "la1", "mu2"))
      fun = function(var, par){
        x = as.numeric(var[par])
        s = as.numeric(step2[par])
        r = rnorm(1, mean = 0, sd = s)
        new.var = x + ifelse(r < -x, -r, r)
        return(new.var)
      }
    }
  } else if(class(step) == "function"){
    # eg, step = function(x) x+runif(length(x), -1, 1)
    fun = function(var, par){
      x = as.numeric(var[par])
      new.var = step(x)
      return(new.var)
    }
  } else if(class(step) == "list"){
    # eg, step = list(function(x) x+runif(1, -1, 1),
    #                 function(x) x+runif(1, -1, 1),
    #                 function(x) x+runif(1, -1, 1),
    #                 function(x) x+rexp(1, 1))
    names(step) = c("b", "mu1", "la1", "mu2")
    fun = function(var, par){
      x = as.numeric(var[par])
      new.var = mapply(function(p, xxx) step[[p]](xxx), p = par, xxx = x)
      return(new.var)
    }
  } else{
    stop("'step' must be one of the three: a vector with the standard deviation of the step; 
         a function; or a list of 4 functions.")
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







