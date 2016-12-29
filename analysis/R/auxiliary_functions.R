





start.output <- function(outputName, logLik, prior, posterior, parameters){
  output <- file(paste(outputName, "_PBD_Bayes.txt", sep = ""), "w")
  
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










