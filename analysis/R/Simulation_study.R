

library(ape)
library(PBD)
library(Rcpp)
library(RcppArmadillo)
library(doParallel)

#output = "~/Desktop/gitHub/protracted_sp/analysis/R/"
output = ""
out_sim = "simulations/"
out_post = "posterior/"

# source("opt_pbd_sim_cpp.R")
# theta0 = list("b1" = c(0.1, 0.3, 0.5, 0.7, 1, 3, 5, 7),
#               "la1" = c(10, 7, 5, 3, 1, 0.5, 0.3, 0.1, 0.05, 0),
#               "b2" = c(0),
#               "mu1" = c(0, 0.01, 0.05, 0.1, 0.2, 0.4, 1, 2),
#               "mu2" = c(0, 0.01, 0.05, 0.1, 0.2, 0.7, 1, 2, 3))
# 
# 
# expand = function(x, default = 1, rm.neg = FALSE){
#   xl = length(x)
#   ll = sapply(x, length)
#   if(length(default) != xl){
#     default = sapply(x, "[[", default)
#   }
#   out = lapply(1:xl, function(y){
#     foo = data.frame(matrix(rep(default, each = ll[y]), ncol = xl))
#     foo[ , y] = x[[y]]
#     foo
#   })
#   out = do.call(rbind, out)
#   out = out[!duplicated(out), ]
#   if(!is.null(names(x))){
#     colnames(out) = names(x)
#   }
#   if(rm.neg){
#     if(is.null(colnames(out))){
#       warning("The columns have NO name!\nI'll assume they're organized as follow: b1 la1 b2 mu1 mu2.")
#       colnames(out) = c("b1", "mu1", "b2", "la1", "mu2")
#     } 
#     if(any(! c("b1", "mu1", "b2", "la1", "mu2") %in% colnames(out))){
#       warning("The names of the columns are incorrect!\nI'll return the full data.frame instead.")
#     } else{
#       bool = (out$b1 - out$mu2) < 0 & (1/out$la1 - out$mu1) < 0
#       out = out[!bool, ]
#     }
#   }
#   
#   return(out)
# }
# parameters = expand(theta0, default = c(1,3,0,0.5,2))
# nrow(parameters)
# 
# tries = 1000
# ntaxa = 128 # 2^(5:10) = c(32, 64, 128, 256, 512, 1024)
# nphy = 30
# sim = mclapply(X = 1:nrow(parameters), FUN = function(s){
#   cat(paste0(Sys.time(), "\t parameters: ", s, "\n"))
#   cat(paste0(parameters[s,]))
#   out1 = lapply(1:nphy, FUN = function(x){
#     try(opt_pbd_sim_cpp(pars = as.matrix(parameters[s, ]), taxa = ntaxa, ntry = tries), silent = TRUE)
#   })
#   out2 = lapply(X = out1, FUN = function(y){
#     if(class(y) != "try-error"){
#       this_age = y$L[1, 3]
#       p = try(get_phylo(y$L, this_age), silent = TRUE)
#       if(class(p) == "phylo"){
#         p = drop.fossil(p)
#         if(Ntip(p) == ntaxa){
#           return(p)
#         }
#       }
#     }
#     return(NULL)
#   })
#   bool = sapply(out2, FUN = function(x){!is.null(x)})
#   names(out2)[bool] = paste0(parameters[s,], collapse = "_")
#   
#   cat("\n\n")
#   save(out1, out2, file = paste0(output, out_sim, "pbd_sim_cpp", paste0(parameters[s,], collapse = "_"), ".RData"))
#   return(list(simulation = out1, phylogeny = out2))
# }, mc.cores = 3)
# save(sim, file = paste0(output, out_sim, "pbd_sim_cpp.RData"))
# 
# # get the phylogenies
# clean = lapply(lapply(sim, "[[", 2), function(x){
#   nn = names(x)
#   if(all(is.na(nn))){
#     return(NULL)
#   } else{
#     return(x[!is.na(nn)])
#   }
# })
# trees = unlist(clean, recursive=FALSE)
# phylos = trees[sapply(trees, function(x) ifelse(is.null(x), FALSE, class(x) == "phylo"))] # filter the parameters for the ones that have 10 phylogenies
# branches = lapply(phylos, branching.times) # get the branches ages
# save(phylos, branches, file = paste0(output, out_sim, "phylogenies.RData"))




load(file = paste0(output, out_sim, "phylogenies.RData"))
breaks = 1:298
cores = 16
rep = 5e+3
chains = 2
namphy = names(phylos)
repl = duplicated(namphy)
for(i in 1:length(namphy)){
  if(repl[i]){
    namphy[i] = paste0(namphy[i], "_phy_", count)
    count = count + 1
  } else{
    namphy[i] = paste0(namphy[i], "_phy_1")
    count = 2
  }
}
files = paste0("Bayes_simulation_study_", namphy)
source("pbd_Bayes.R")
prior_b = function(b){
  #dexp(b, rate = 1/sim.pars[1], log = TRUE)
  dlnorm(b, meanlog = 0, sdlog = 1, log = TRUE)
}
prior_mu1 = function(mu1){
  dexp(mu1, rate = 0.5, log = TRUE)
  #dlnorm(mu1, meanlog = 0, sdlog = 1, log = TRUE)
}
prior_la1 = function(la1){
  dlnorm(la1, meanlog = 1, sdlog = 1, log = TRUE)
}
prior_mu2 = function(mu2){
  dexp(mu2, rate = 0.2, log = TRUE)
  #dlnorm(mu2, meanlog = .25, sdlog = .5, log = TRUE)
}
mcmapply(FUN = function(phy, brchs, np, this_file){
  for(i in 1:chains){
    pbd_Bayes(brts = brchs/max(brchs),
              initparsopt = rep(i, 4),
              prior_b = prior_b, prior_mu1 = prior_mu1, 
              prior_la1 = prior_la1, prior_mu2 = prior_mu2,
              step = c(2.5, 2.5, 3, 3), 
              rep = rep, 
              file = paste0(output, out_post, this_file, ".", i))
  }
}, phy = phylos[breaks], brchs = branches[breaks], np = namphy[breaks], this_file = files[breaks],
SIMPLIFY = FALSE, mc.cores = cores)

curve(prior_b, 0, 10)
curve(prior_la1, 0, 10, add=T, lty=2)
curve(prior_mu2, 0, 10, add=T, lty=3)
curve(prior_mu1, 0, 10, add=T, lty=4)
plot.s(c(2.5, 2.5, 3, 3))
plot.s = function(vec, lim = 10){
  fff = function(x) dnorm(x, mean = 0, sd = vec[1])
  curve(fff, 0, lim)
  for(i in 2:3){
    fff = function(x) dnorm(x, mean = 0, sd = vec[i])
    curve(fff, 0, lim, add=T, lty=i)
  }
}

logLik_fun = opt_loglik(brts = branches[[100]])
ff = function(x) logLik_fun(pars1 = c(5, 0.5, x))
la1 = seq(0, 7, by = 0.1)
mu2 = seq(0, 5, by = 0.1)
z = matrix(nrow = length(la1), ncol = length(mu2))
for(i in 1:length(la1)){
  z[i, ] = sapply(mu2, function(y) ff(c(la1[i], y)))
}

contour(z)
