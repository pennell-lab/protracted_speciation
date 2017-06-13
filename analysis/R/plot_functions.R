
require(reshape2)
require(ggplot2)
require(ape)
require(gridExtra)

check_MCMC = function(mat, mode = c("all", "acceptance"), Ncore = 2){
  require(dplyr)
  require(parallel)
  if(mode == "all"){
    out = data.frame(t(sapply(mat, function(x){
      ll = nrow(x)
      foo = x %>% summarize(prior = sum(is.na(prior))/ll, 
                            logLik = sum(is.na(logLik))/ll,
                            posterior = sum(is.na(posterior))/ll,
                            converged = sum(converged == FALSE, na.rm = TRUE)/ll,
                            accep.rate = sum(accepted, na.rm = T)/ll)
      as.numeric(foo)
    })))
    colnames(out) = c("prior", "logLik","posterior","converged","accep.rate")
    
  }else if(grepl(mode, "acceptance")){
    accepted = lapply(mat, function(x) x$accepted)
    chsn = mclapply(mat, function(x){
      foo = sapply(x$proposed.par, strsplit, split = "_")
      foo = do.call(rbind, lapply(foo, function(x) as.numeric(c("b"%in%x,"mu1"%in%x,"la1"%in%x,"mu2"%in%x))))
      foo
    }, mc.cores = Ncore)
    
    out = data.frame(t(mcmapply(function(m, b){
      foo = apply(m, MARGIN = 2, function(x) sum(x*b)/sum(x))
      foo2 = sum(b)/length(b)
      c(foo2, foo)
    }, m = chsn, b = accepted, mc.cores = Ncore)))
    colnames(out) = c("all","b","mu1","la1","mu2")
    
  } else{
    stop("'mode' must be either 'all' - to check all - or 'acceptance' - to check the acceptance ratio per parameter.")
  }
  
  return(out)
}

get_param = function(output = "",
                     out_post = "",
                     initial.patt = "",
                     end.patt = "_PBD_Bayes.txt",
                     patt = ".txt",
                     sep.csv = "\t"){
  files = list.files(path = paste0(output, out_post), pattern = patt)
  list.names = gsub(end.patt, "", gsub(initial.patt, "", files))
  par.factor = gsub("_phy", "", list.names)
  parameters = .get_these_par(par.factor, files)
  
  return(parameters)
}
.get_these_par = function(these_par, files = NULL, full = TRUE){
  param = do.call(rbind, strsplit(these_par, split = "_"))
  if(full){
    param = cbind(param, sapply(strsplit(param[ , 6], split = "\\."), "[[", 2))  # adds a cloumn with the chain number
    param[ , 6] = as.integer(param[ , 6])
  }
  param = apply(param, 2, as.numeric)
  if(full){
    parameters = data.frame("phylogeny" = param[ , 6],
                            "chain" = param[ , 7],
                            "r1" = param[ , 1] - param[ , 4],
                            "r2" = (1 / param[ , 2]) - param[ , 5],
                            "b" = param[ , 1],
                            "la1" = param[ , 2],
                            "b2" = param[ , 3],
                            "mu1" = param[ , 4],
                            "mu2" = param[ , 5],
                            "control" = apply(param[ , 1:6], 1, paste0, collapse = "_"),
                            "file" = files,
                            stringsAsFactors = FALSE
    )
    parameters = parameters[order(parameters$b,parameters$la1,parameters$mu1,parameters$mu2, parameters$phylogeny, parameters$chain), ]
    parameters$control = as.factor(parameters$control)
  } else{
    parameters = data.frame("b" = param[ , 1],
                            "la1" = param[ , 2],
                            "mu1" = param[ , 4],
                            "mu2" = param[ , 5],
                            stringsAsFactors = FALSE
    )
    parameters = parameters[order(parameters$b,parameters$la1,parameters$mu1,parameters$mu2), ]
  }
  return(parameters)
}

read.BPBD = function(samp = 100, 
                     burnin = 0.3, 
                     chains = NULL,
                     folder = "",
                     par = NULL,
                     Ncore = NULL,
                     ...){
  require(coda)
  require(doParallel)
  
  if(is.null(chains)){
    stop("number of 'chains' must be provided!")
  }
  if(is.null(par)){
    par = get_param(...)
    .f = function(output, out_post, ...) paste0(output, out_post)
    folder = .f(...)
  }
  if(is.null(Ncore)){
    Ncore = max(detectCores()-1, 1)
  }
  
  Nloop = length(levels(par$control))
  my.mcmc = my.par = rep(list(NA), Nloop)
  
  
  out = lapply(X = 1:Nloop, FUN = function(k){
    ind = which(par$control == levels(par$control)[k])
    if(length(ind) == chains){
      # files to be used
      thisfiles = par$file[ind]
      # read the files
      raw = mclapply(paste0(folder, thisfiles), read.csv, sep = "\t", row.names = NULL, as.is = TRUE, mc.cores = Ncore)
      reads = mclapply(raw, function(x) return(x[-1, ]), mc.cores = Ncore)
      # gets the number of steps for each independent run
      nrep = sapply(reads, nrow)
      # check nrow of reads and drops the 'wrong' ones
      if(any(nrep != nrep[1])){
        warning("The following files were incomplete and were not read:\n", paste(par$file[ind], collapse = "\n"))
      } else{
        # gets the rows that will be sampled
        burn = ifelse(burnin > 1, burnin, ceiling(burnin*nrep[1]))
        sel = seq(burnin, nrep[1], by = samp)
        # combine them into one single data.frame to store all the info
        pp = mclapply(reads, function(x) return(x[sel, ]), mc.cores = Ncore)
        my.par = do.call(rbind, pp)
        
        # transforms all runs into "mcmc" objects
        # then combines the runs by region into "mcmc.list" objects
        aux = mclapply(lapply(reads, "[", 4:7), mcmc, start = burnin, thin = samp, mc.cores = Ncore)
        my.mcmc = mcmc.list(aux)
        # to name the simulated parameters
        names(my.mcmc) = paste(levels(par$control)[k], 1:length(ind), sep = ".")
        
        return(list(my.mcmc, my.par))
      }
    } else{
      warning("The following files do not have the correct number of chains:\n", paste(par$file[ind], collapse = "\n"))
    }
    return(list(NA, NA))
  })
  
  my.mcmc = lapply(out, "[[", 1)
  my.par = lapply(out, "[[", 2)
  
  names(my.par) = levels(par$control)
  
  if(  any(sapply(my.par, function(x) sum(is.na(x))))  ){
    my.mcmc = my.mcmc[!sapply(my.mcmc, function(x) sum(is.na(x)))]
    my.par = my.par[!sapply(my.par, function(x) sum(is.na(x)))]
  }
  
  return(list(mcmc = my.mcmc, samples = my.par))
}


plot.prior.post = function(samples, priors, par, ...){
  if(is.null(par)){
    par = as.numeric(strsplit(names(samples)[1], split = "_")[[1]])[c(1,4,2,5)]
  }
  if(is.character(priors)){
    if(priors %in% c("lnorm", "lognormal", "lognorm")){
      prior_b = function(b){
        dlnorm(b, meanlog = par[1], sdlog = 1, log = TRUE)
      }
      prior_mu1 = function(mu1){
        dlnorm(mu1, meanlog = par[4], sdlog = 1, log = TRUE)
      }
      prior_la1 = function(la1){
        dlnorm(la1, meanlog = par[2], sdlog = 1, log = TRUE)
      }
      prior_mu2 = function(mu2){
        dlnorm(mu2, meanlog = par[5], sdlog = 1, log = TRUE)
      }
    }
    if(priors %in% c("exp", "exponencial")){
      prior_b = function(b){
        dexp(b, rate = 1/par[1], log = TRUE)
      }
      prior_mu1 = function(mu1){
        dexp(mu1, rate = 1/max(0.01, par[4]), log = TRUE)
      }
      prior_la1 = function(la1){
        dlnorm(la1, meanlog = par[2], sdlog = 1, log = TRUE)
      }
      prior_mu2 = function(mu2){
        dexp(mu2, rate = 1/max(0.01, par[5]), log = TRUE)
      }
    }
    priors = try(c(prior_b, prior_mu1, prior_la1, prior_mu2))
    if(class(priors) == "try-error"){
      stop("'priors' must be either 4 funcitions or a charcter string choosing between 'lognormal' or 'exponencial'.")
    }
  }
  samples = do.call(rbind, samples)
  medians = apply(samples, 2, median)
  ggplot(melt(samples), aes(value)) + geom_density() + geom_vline(xintercept = medians) + facet_grid(. ~ variable)
  plot(density(this_vector), main = colnames(samples[[1]])[i], ...)
  abline(v = median(this_vector))
  abline(v = par[i], col = "green")
  curve(exp(priors[[i]](x)), add = TRUE, col = "red")
  legend("topright", legend = c("prior", "posterior"),
         lty = c(1, 1), col = c("red", "black"),
         cex = 0.7)
}




#############################################################
#   Function to inspect for possible correlations between
# parameters of a Bayesian model. Better than the function
# coda::pairs.
#
#   Developed by Florian Hartig
#############################################################
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="blue4", ...)
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, method = "spearman"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex * r)
}

betterPairs <- function(YourData){
  require(IDPmisc)
  return(pairs(YourData, lower.panel=function(...) {par(new=TRUE);ipanel.smooth(...)}, diag.panel=panel.hist, upper.panel=panel.cor))
}


simulatedXposterior = function(samples, 
                               parameters = NULL,
                               statistics = "mode",
                               shapes0 = 1,
                               colour0 = "black",
                               size0 = 1,
                               shapes1 = 19,
                               colour1 = "black",
                               size1 = 1.5,
                               rates = FALSE,
                               plot = TRUE,
                               save = NULL,
                               ...){
  require(ggplot2)
  require(reshape2)
  require(gridExtra)
  
  if(any(sapply(samples, ncol) < 4)){
    stop("'samples' must have at least 4 columns; one for each parameter.")
  }
  if(any(sapply(samples, ncol) > 4)){
    warning("'samples' has more than 4 columns; I will try to extract the columns for each parameter.")
    samples = try(lapply(samples, function(x) x[c("b", "mu1", "la1", "mu2")]))
    if(class(samples) == "try-error"){
      stop("I could not extract the columns for the 4 parameters. Check name of the columns!")
    }
  }
  if(is.null(parameters)){
    warning("No parameters given, I'll try to extract it from the names of 'samples'.")
    parameters = .get_these_par(these_par = names(samples), full = FALSE)
  }
  
  
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  stat.fun = switch (statistics,
                     "mode" = Mode,
                     "median" = median,
                     "mean" = mean,
                     "variance" = var,
                     "var" = var,
                     "sd" = sd
  )
  
  if(rates){
    rates = data.frame(t(sapply(samples, function(x) structure(c(stat.fun(x$b - x$mu2), stat.fun(x$la1 - x$mu1)), names = c("samp.r1","samp.r2")) )))
    samples = data.frame(t(sapply(samples, function(x) apply(x, MARGIN = 2, FUN = stat.fun))))
    samples = cbind(samp = samples, rates, parameters[c("r1", "r2", "b", "la1", "mu1", "mu2")])
    melt0 = melt(samples, measure.vars = c("samp.b", "samp.mu1", "samp.la1", "samp.mu2", "samp.r1", "samp.r2"))
    
  } else{
    samples = data.frame(t(sapply(samples, function(x) apply(x, MARGIN = 2, FUN = stat.fun))))
    samples = cbind(samp = samples, parameters[c("b", "la1", "mu1", "mu2")])
    melt0 = melt(samples, measure.vars = c("samp.b", "samp.mu1", "samp.la1", "samp.mu2"))
  }
  
  # B
  df.b = get_fixed_rows(melt0, "b")
  dfPb = get_fixed_points(unique(df.b$b), "b", names = paste0("samp.",c("b", "la1", "mu1", "mu2")))
  b = ggplot(df.b, aes(b, value), ...) + 
    geom_point(data = dfPb, aes(x=x, y=y), shape = shapes1, color = colour1, size = size1) +
    geom_point(shape = shapes0, colour = colour0, size = size0) +
    facet_grid(. ~ variable) #geom_jitter(width = .3, height = 0)
  # LA1
  df.la1 = get_fixed_rows(melt0, "la1")
  dfPla1 = get_fixed_points(unique(df.la1$la1), "la1", names = paste0("samp.",c("b", "la1", "mu1", "mu2")))
  la1 = ggplot(df.la1, aes(la1, value), ...) + 
    geom_point(data = dfPla1, aes(x=x, y=y), shape = shapes1, color = colour1, size = size1) +
    geom_point(shape = shapes0, colour = colour0, size = size0) +
    facet_grid(. ~ variable) #geom_jitter(width = .3, height = 0)
  # MU1
  df.mu1 = get_fixed_rows(melt0, "mu1")
  dfPmu1 = get_fixed_points(unique(df.mu1$mu1), "mu1", names = paste0("samp.",c("b", "la1", "mu1", "mu2")))
  mu1 = ggplot(df.mu1, aes(mu1, value), ...) + 
    geom_point(data = dfPmu1, aes(x=x, y=y), shape = shapes1, color = colour1, size = size1) +
    geom_point(shape = shapes0, colour = colour0, size = size0) +
    facet_grid(. ~ variable) #geom_jitter(width = .3, height = 0)
  # MU2
  df.mu2 = get_fixed_rows(melt0, "mu2")
  dfPmu2 = get_fixed_points(unique(df.mu2$mu2), "mu2", names = paste0("samp.",c("b", "la1", "mu1", "mu2")))
  mu2 = ggplot(df.mu2, aes(mu2, value), ...) + 
    geom_point(data = dfPmu2, aes(x=x, y=y), shape = shapes1, color = colour1, size = size1) +
    geom_point(shape = shapes0, colour = colour0, size = size0) +
    facet_grid(. ~ variable) #geom_jitter(width = .3, height = 0)
  
  
  
  if(rates){
    out = grid.arrange(b, la1, mu1, mu2, r1, r2)
  } else{
    out = grid.arrange(b, la1, mu1, mu2)
  }
  
  if(plot){
    out
  }
  if(!is.null(save)){
    ggsave(filename = save, plot = out)
  }
  
  return(out)
}

get_fixed_rows = function(df, id, fixed = c(0.1, 0.05, 0.01, 0.05), names = c("b", "la1", "mu1", "mu2")){
  names(fixed) = names
  if(class(id) == "character"){
    id = switch (id,
                 "b" = 1,
                 "la1" = 2,
                 "mu1" = 3,
                 "mu2" = 4
    )
  }
  fixed = fixed[-id]
  bool = apply(df[names(fixed)], 1, function(x) sum(x==fixed)==3)
  return(df[bool, ])
}
get_fixed_points = function(varies, id, fixed = c(0.1, 0.05, 0.01, 0.05), names = c("b", "la1", "mu1", "mu2")){
  names(fixed) = names
  ll = length(varies)
  if(class(id) == "character"){
    id = switch (id,
                 "b" = 1,
                 "la1" = 2,
                 "mu1" = 3,
                 "mu2" = 4
    )
  }
  out = data.frame(x = rep(varies, 4), y = c(varies, rep(fixed[-id], each = ll)), var = rep(c(names[id], names[-id]), each = ll))
  return(out)
}


