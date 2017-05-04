

library(dplyr)
library(ape)
library(nlme)
library(doParallel)
library(parallel)

cores = 14L

# setwd("~/Desktop/gitHub/")
#path.data = "protracted_sp/analysis/data/"
path.data = "data/"
#path.code = "protracted_sp/SSP_DR/R/"
path.code = "data/"
#path.out = "protracted_sp/SSP_DR/output/"
path.out = "output/"

load("output/DR_age_estimate.RData")
dat = read.csv("output/dat.csv", as.is = TRUE)


backbone = c("Ericson", "Hackett")
# ericson = read.tree(paste0(path.data, "Jetz_etal/EricsonStage2_0001_1000/AllBirdsEricson1_100.tre"))
# hackett = read.tree(paste0(path.data, "Jetz_etal/HackettStage2_0001_1000/AllBirdsHackett1_100.tre"))
ericson = read.tree(paste0(path.data, "AllBirdsEricson1_100.tre"))
hackett = read.tree(paste0(path.data, "AllBirdsHackett1_100.tre"))
trees = list(ericson, hackett)
names(trees) = backbone

# estimate the PIC for dr and sp age
botero = read.csv(paste0(path.data, "Botero14.csv"), as.is = TRUE)
# select only birds (ie, drop mammals)
birds.raw = botero %>% 
  filter(TAXON == "birds") %>%
  arrange(SPECIES)
# select only the species with data
absent = trees[[1]][[1]]$tip.label[which(!trees[[1]][[1]]$tip.label %in% birds.raw$SPECIES)]
pruned = lapply(trees, function(x){
  # there's 1 sp present in Botero dataset that doesn't appear in the phylogeny "Pseudohirundo_griseopyga"
  mclapply(x, drop.tip, tip = c(absent, "Pseudohirundo_griseopyga"), mc.cores = cores)
})
save(pruned, file = "output/pruned_100_phy.RData")




gls.age.dr.log = mapply(function(aa, dd, pp){
  out1 = mcmapply(function(aaa, ddd, ppp){
    gls(log(ddd) ~ log(aaa), correlation = corPagel(1, ppp, fixed = FALSE))
  }, aaa = aa, ddd = dd, ppp = pp, SIMPLIFY = FALSE, mc.cores = cores)
  foo = lapply(out1, function(x) if(class(x) == "gls"){anova(x)})
  out2 = t(mapply(function(x, y){
    if(class(x) == "gls"){
      data.frame(intercept = coef(x)[1], data = coef(x)[2], lambda = attributes(x$apVar)$Pars[1], residual.se = x$sigma, logLik = x$logLik, p.inter = y$`p-value`[1], p.data = y$`p-value`[2])
    } else{
      data.frame(rep(NA, 7))
    }
    }, x = out1, y = foo))
  return(list(model = out1, summary = out2))
}, aa = age, dd = dr, pp = trees, SIMPLIFY = FALSE)
save(gls.age.dr.log, file = "output/age_DR_gls.RData")



togls = dat$ssp; togls = names(togls) = dat$sp

gls.ssp.dr.log = mapply(function(dd, pp){
  out1 = mcmapply(function(ddd, ppp){
    ddd = ddd[!names(ddd) %in% absent]
    gls(log(togls) ~ log(ddd), correlation = corPagel(1, ppp, fixed = FALSE))
  }, ddd = dd, ppp = pp, SIMPLIFY = FALSE, mc.cores = cores)
  foo = lapply(out1, function(x) if(class(x) == "gls"){anova(x)})
  out2 = t(mapply(function(x, y){
    if(class(x) == "gls"){
      data.frame(intercept = coef(x)[1], data = coef(x)[2], lambda = attributes(x$apVar)$Pars[1], residual.se = x$sigma, logLik = x$logLik, p.inter = y$`p-value`[1], p.data = y$`p-value`[2])
    } else{
      data.frame(rep(NA, 7))
    }
  }, x = out1, y = foo))
  return(list(model = out1, summary = out2))
}, dd = dr, pp = pruned, SIMPLIFY = FALSE)
save(gls.ssp.dr, file = "output/SSP_DR_gls.RData")

gls.ssp.age.log = mapply(function(aa, pp){
  out1 = mcmapply(function(aaa, ppp){
    ddd = ddd[!names(ddd) %in% absent]
    gls(log(togls) ~ log(aaa), correlation = corPagel(1, ppp, fixed = FALSE))
  }, aaa = aa, ppp = pp, SIMPLIFY = FALSE, mc.cores = cores)
  foo = lapply(out1, function(x) if(class(x) == "gls"){anova(x)})
  out2 = t(mapply(function(x, y){
    if(class(x) == "gls"){
      data.frame(intercept = coef(x)[1], data = coef(x)[2], lambda = attributes(x$apVar)$Pars[1], residual.se = x$sigma, logLik = x$logLik, p.inter = y$`p-value`[1], p.data = y$`p-value`[2])
    } else{
      data.frame(rep(NA, 7))
    }
  }, x = out1, y = foo))
  return(list(model = out1, summary = out2))
}, aa = age, pp = pruned, SIMPLIFY = FALSE)
save(gls.ssp.age, file = "output/SSP_age_gls.RData")
