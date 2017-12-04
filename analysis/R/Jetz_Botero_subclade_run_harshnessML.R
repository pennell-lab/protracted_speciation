
#   Subclade harshness/clean/3par analyses for pbd ML
# as well as tradinal spec- and extinction rates
# Possible effect of AGE


library(parallel)
library(PBD)


setwd("~/Desktop/gitHub/protracted_sp/")
input_data = "analysis/data/"
input_code = "analysis/R/"
ouput_data = "analysis/data/"



# branches, branches.unity, branchesHarsh, branchesHarsh.unity, branchesHarshClean, branchesHarshClean.unity
load(file = paste0(input_data, "Jetz_Botero14_30min_0.3mixed_0multilat_branches.RData"))

# data.frame with info for harshness
propMixedHarsh = read.csv(file = paste0(input_data, "Jetz_Botero14_propMixedHarsh.csv"))

# Maximum likelihood
################ ONLY phylogeny 1
# # run the analysis
# Ncore = 3
# # selects only this phylogeny
# this_data = propMixedHarsh[propMixedHarsh$phy == 1, ]
# # gets the indices for the "Stable" and "Variable" subclades (ie, excludes "Regular")
# ind = which(this_data$which.max != "Regular")
# 
# # ML for the harshness full phylogenies with 4 parameters
# harshML = mclapply(ind, FUN = function(xxx){
#   pbd_ML(brts = branchesHarsh.unity[[1]][[xxx]], # set of branching times of a phylogeny
#          initparsopt = rep(1, 4), # initial values of the parameters that must be optimized
#          exteq = 0, # incipient have the same (1) or different (0) ext rate as good species
#          btorph = 0, # likelihood is for the branching times (0) or the phylogeny (1)
#          soc = 2) # first element of the branching times is the stem (1) or the crown (2) age
# }, mc.cores = Ncore)
# harshML = do.call(rbind, harshML)
# harshML$type = "harshness"
# harshML$est = ind
# harshML$category = propMixedHarsh$which.max[ind]
# 
# # ML for the harshness pruned phylogenies with 4 parameters
# harshMLclean = mclapply(ind, FUN = function(xxx){
#   pbd_ML(brts = branchesHarshClean.unity[[1]][[xxx]], # set of branching times of a phylogeny
#          initparsopt = rep(1, 4), # initial values of the parameters that must be optimized
#          exteq = 0, # incipient have the same (1) or different (0) ext rate as good species
#          btorph = 0, # likelihood is for the branching times (0) or the phylogeny (1)
#          soc = 2) # first element of the branching times is the stem (1) or the crown (2) age
# }, mc.cores = Ncore)
# harshMLclean = do.call(rbind, harshMLclean)
# harshMLclean$type = "harshnessClean"
# harshMLclean$est = ind
# harshMLclean$category = propMixedHarsh$which.max[ind]
# 
# # ML for the harshness pruned phylogenies with 3 parameters
# harshMLclean3par = mclapply(ind, FUN = function(xxx){
#   pbd_ML(brts = branchesHarshClean.unity[[1]][[xxx]], # set of branching times of a phylogeny
#          initparsopt = rep(1, 3), # initial values of the parameters that must be optimized
#          exteq = 1, # incipient have the same (1) or different (0) ext rate as good species
#          btorph = 0, # likelihood is for the branching times (0) or the phylogeny (1)
#          soc = 2) # first element of the branching times is the stem (1) or the crown (2) age
# }, mc.cores = Ncore)
# harshMLclean3par = do.call(rbind, harshMLclean3par)
# harshMLclean3par$type = "harshClean3par"
# harshMLclean3par$est = ind
# harshMLclean3par$category = propMixedHarsh$which.max[ind]
# 
# save(harshML, harshMLclean, harshMLclean3par, file = "analysis/data/subclade_harshness/subclade_harshness_ML.RData")



load(file = "analysis/data/subclade_harshness/subclade_harshness_ML.RData")

library(dplyr)
library(ggplot2)
harshAll = rbind(harshML, harshMLclean, harshMLclean3par)
harshAllmelt = harshAll %>% 
  select(type, est, category, b, mu_1, lambda_1, mu_2) %>%
  reshape2::melt(., id.vars = 1:3)

ggplot(harshAllmelt, aes(value, fill = category)) + 
  geom_histogram(aes(y = ..density..)) +
  facet_grid(variable ~ type) +
  xlim(0, 20)





# # Traditional models: fits by maximum likelihood a birth-death model to the branching times computed from a phylogenetic tree using the method of Nee et al. (1994).
# ################ ONLY phylogeny 1
# source("analysis/R/auxiliary_functions.R")
# # run the analysis
# Ncore = 3
# ##### Subclades from the phylogenies; instead of branches
# # these_subclades, subclades, these_subcladesHarsh, subcladesHarsh, subcladesHarshClean
# load(file = "analysis/data/Jetz_Botero14_30min_0.3mixed_0multilat_subclades.RData")
# # selects only this phylogeny
# this_data = propMixedHarsh[propMixedHarsh$phy == 1, ]
# # gets the indices for the "Stable" and "Variable" subclades (ie, excludes "Regular")
# ind = which(this_data$which.max != "Regular")
# 
# # ML for the harshness full phylogenies 
# harshSpcExt = mclapply(ind, FUN = function(xxx){
#   foo = birthdeath(phy = subcladesHarsh[[1]][[xxx]])
#   get_SpecExt4birthdeath(foo$para)
# }, mc.cores = Ncore)
# harshSpcExt = do.call(rbind, harshSpcExt)
# harshSpcExt = data.frame(speciation = harshSpcExt[ , 1],
#                          extinction = harshSpcExt[ , 2],
#                          type = "harshness",
#                          est = ind,
#                          category = propMixedHarsh$which.max[ind])
# 
# # ML for the harshness pruned phylogenies
# harshSpcExtClean = mclapply(ind, FUN = function(xxx){
#   foo = birthdeath(phy = subcladesHarshClean[[1]][[xxx]])
#   get_SpecExt4birthdeath(foo$para)
# }, mc.cores = Ncore)
# harshSpcExtClean = do.call(rbind, harshSpcExtClean)
# harshSpcExtClean = data.frame(speciation = harshSpcExtClean[ , 1],
#                               extinction = harshSpcExtClean[ , 2],
#                               type = "harshnessClean",
#                               est = ind,
#                               category = propMixedHarsh$which.max[ind])
# 
# save(harshSpcExt, harshSpcExtClean, file = "analysis/data/subclade_harshness/subclade_harshness_SpecExt.RData")


load(file = "analysis/data/subclade_harshness/subclade_harshness_SpecExt.RData")

library(ggplot2)
harshAllSpecExt = rbind(harshSpcExt, harshSpcExtClean)
harshAllSpecExtmelt = reshape2::melt(harshAllSpecExt, id.vars = 3:5)

ggplot(harshAllSpecExtmelt, aes(value, fill = category)) +
  geom_histogram(aes(y = ..density..)) +
  facet_grid(variable ~ type) 











