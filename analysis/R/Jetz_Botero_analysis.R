

library(ape)
library(PBD)
library(doParallel)

#output = "~/Desktop/gitHub/protracted_sp/analysis/R/"
out_post = "posterior/"


# load(file = "analysis/data/Botero14.RData")
load(file = "Botero14.RData")

# load reduced set of Jetz trees (only 10)
#jetz = read.tree("analysis/data/Jetz_etal/EricsonStage2_0001_1000/AllBirdsEricson1_reduced.tre")
# jetz = read.tree("AllBirdsEricson1_reduced.tre")
# check which sp are present in the table and drops the ones that are not
# there's 3300sp missing from 9993sp total
# miss = lapply(jetz, function(x) which(!x$tip.label %in% birds$SPECIES))
# jj = mapply(function(p, m) drop.tip(p, m), p = jetz, m = miss, SIMPLIFY = FALSE)
# sapply(jj, Ntip)
# 
# regions = sapply(b1names, paste0, collapse = "+")
# # select the species that occur in that given region
# sp = lapply(b1names, function(xx) birds$SPECIES[birds$LAT.RANGE %in% xx])
# branches = mclapply(1:length(jetz), function(i){
#   # get the list of species that are NOT in the list
#   miss = lapply(sp, function(x) which(!jetz[[i]]$tip.label %in% x))
#   # drop the tips, keeping only what occurs in the region
#   jj = lapply(miss, drop.tip, phy = jetz[[i]])
#   # get the branch length for the each region
#   out = lapply(jj, branching.times)
#   names(out) = regions
#   return(out)
# }, mc.cores = 12)
# branches = lapply(regions, function(x){
#   out = do.call(rbind, lapply(branches, "[[", x))
#   rownames(out) = paste0("phy", 1:10)
#   out
# })
# save(branches, file = "analysis/data/Jetz_Botero_branches.RData")


#load(file = "analysis/data/Jetz_Botero_branches.RData")
load(file = "Jetz_Botero_branches.RData")
# run the analysis
regions = sapply(b1names, paste0, collapse = "+")
rep = 5e+5
chains = 2
limits = c(10, 10, 20, 20) # upper bound for parameters; order = "b", "mu1", "la1", "mu2"
rates = structure(c(0.1,0.4,0.1), names = c("la1", "mu1", "mu2")) # OLD = structure(c(0.1,0.1,0.4,0.1), names = c("b", "la1", "mu1", "mu2"))
prior_mu1 = function(mu1){
  dexp(mu1, rate = rates["mu1"], log = TRUE)
}
prior_la1 = function(la1){
  dexp(la1, rate = rates["la1"], log = TRUE)
  #dlnorm(la1, meanlog = 1, sdlog = 1, log = TRUE)
}
prior_mu2 = function(mu2){
  dexp(mu2, rate = rates["mu2"], log = TRUE)
}

#source("analysis/R/pbd_Bayes.R")
source("pbd_Bayes.R")


lapply(1:length(branches), function(i){
  files = paste0("Jetz_Botero_", regions[i], "_phy", 1:nrow(branches[[i]]))
  # create the prior with the estimated parameters, for the given number of subspecies
  prior_b = function(b){
    dexp(b, rate = prior4b1.exp[regions[i]], log = TRUE)
    #dlnorm(b, meanlog = 0, sdlog = 1, log = TRUE)
  }
  
  mclapply(1:length(files), FUN = function(xxx){
    for(chn in 1:chains){
      pbd_Bayes(brts = branches[[i]][xxx, ],
                initparsopt = rep(chn, 4),
                prior_b = prior_b, prior_mu1 = prior_mu1, 
                prior_la1 = prior_la1, prior_mu2 = prior_mu2,
                step = c(0.5, 0.5, 1.5, 1.5), 
                upper = limits, 
                rep = rep, 
                file = paste0(out_post, paste(files[xxx], chn, sep = ".")))
    }
  }, mc.cores = 10)
})


## ML
# get_ML = function(brn, tries = 100, ...){
#   while(tries > 0){
#     foo = pbd_ML(brts = brn, ...)
#     if(foo$conv == 0){
#       break()
#     }
#     tries = tries-1
#     start = as.numeric(foo[1:4])
#   }
#   return(foo)
# }
# ml.estimates = lapply(1:length(branches), function(i){
#   out = mclapply(1:nrow(branches[[i]]), FUN = function(xxx){
#     get_ML(branches[[i]][xxx, ], initparsopt = rep(0.5, 4), exteq = 0, btorph = 0, soc = 2, verbose = FALSE)
#   }, mc.cores = 4)
#   do.call(rbind, out)
# })
# ml3estimates = lapply(1:length(branches), function(i){
#   out = mclapply(1:nrow(branches[[i]]), FUN = function(xxx){
#     get_ML(branches[[i]][xxx, ], initparsopt = rep(0.5, 3), exteq = 1, btorph = 0, soc = 2, verbose = FALSE)
#   }, mc.cores = 4)
#   do.call(rbind, out)
# })
# names(ml.estimates) = names(ml3estimates) = sapply(b1names, paste0, collapse = "+")
# save(ml.estimates, ml3estimates, file = "analysis/data/Jetz_Botero_ML.RData")




# my.birthdeath = function(brn) {
#   N <- length(brn)
#   x <- c(NA, brn)
#   dev <- function(a, r) {
#     if (r < 0 || a > 1) 
#       return(1e+100)
#     -2 * (lfactorial(N - 1) + (N - 2) * log(r) + r * sum(x[3:N]) + 
#             N * log(1 - a) - 2 * sum(log(exp(r * x[2:N]) - a)))
#   }
#   out <- nlm(function(p) dev(p[1], p[2]), c(0.1, 0.2), hessian = TRUE)
#   if (out$estimate[1] < 0) {
#     out <- nlm(function(p) dev(0, p), 0.2, hessian = TRUE)
#     para <- c(0, out$estimate)
#     inv.hessian <- try(solve(out$hessian))
#     se <- if (class(inv.hessian) == "try-error") 
#       NA
#     else sqrt(diag(inv.hessian))
#     se <- c(0, se)
#   }
#   else {
#     para <- out$estimate
#     inv.hessian <- try(solve(out$hessian))
#     se <- if (class(inv.hessian) == "try-error") 
#       c(NA, NA)
#     else sqrt(diag(inv.hessian))
#   }
#   Dev <- out$minimum
#   foo <- function(which, s) {
#     i <- 0.1
#     if (which == 1) {
#       p <- para[1] + s * i
#       bar <- function() dev(p, para[2])
#     }
#     else {
#       p <- para[2] + s * i
#       bar <- function() dev(para[1], p)
#     }
#     while (i > 1e-09) {
#       while (bar() < Dev + 3.84) p <- p + s * i
#       p <- p - s * i
#       i <- i/10
#     }
#     p
#   }
#   CI <- mapply(foo, c(1, 2, 1, 2), c(-1, -1, 1, 1))
#   dim(CI) <- c(2, 2)
#   names(para) <- names(se) <- rownames(CI) <- c("d/b", "b-d")
#   colnames(CI) <- c("lo", "up")
#   obj <- list(tree = deparse(substitute(phy)), N = N, dev = Dev, 
#               para = para, se = se, CI = CI)
#   class(obj) <- "birthdeath"
#   obj
# }
# bd.model = lapply(1:length(branches), function(i){
#   out = mclapply(1:nrow(branches[[i]]), FUN = function(xxx){
#     my.birthdeath(branches[[i]][xxx, ])
#   }, mc.cores = 3)
#   do.call(rbind, lapply(out, "[[", 4))
# })
# get_bd = function(div, sub){
#   birth = sub / (1 - div)
#   death = birth * div
#   return(structure(c(birth, death), names = c("birth", "death")))
# }
# bd.estimates = do.call(rbind, lapply(bd.model, function(y) t(apply(y, 1, function(k) get_bd(k[1], k[2])))))
# bd.estimates = data.frame(bd.estimates, regions = rep(sapply(b1names, paste0, collapse = "+"), each = 10))
# save(bd.estimates, file = "analysis/data/Jetz_Botero_BirthDeath.RData")
