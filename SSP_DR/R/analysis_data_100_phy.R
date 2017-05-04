

library(ape)
library(geiger)
library(phytools)
library(dplyr)
library(doParallel)

cores = 14

# setwd("~/Desktop/gitHub/")
#path.data = "protracted_sp/analysis/data/"
path.data = "data/"
#path.code = "protracted_sp/SSP_DR/R/"
path.code = "data/"
#path.out = "protracted_sp/SSP_DR/data/"
path.out = "output/"

backbone = c("Ericson", "Hackett")
# ericson = read.tree(paste0(path.data, "Jetz_etal/EricsonStage2_0001_1000/AllBirdsEricson1_100.tre"))
# hackett = read.tree(paste0(path.data, "Jetz_etal/HackettStage2_0001_1000/AllBirdsHackett1_100.tre"))
ericson = read.tree(paste0(path.data, "AllBirdsEricson1_100.tre"))
hackett = read.tree(paste0(path.data, "AllBirdsHackett1_100.tre"))
trees = list(ericson, hackett)
names(trees) = backbone

# estimate the DR and sp age
source(paste0(path.code, "auxiliary_functions.R"))
dr = lapply(trees, function(x) lapply(x, get_DR))
age = lapply(trees, function(x) lapply(x, get_sp_age))
names(dr) = names(age) = backbone
DR_age = mapply(function(dd, aa){
  out = data.frame(phy = rep(1:length(dd), each = 9993),
                   sp = unlist(lapply(dd, function(x) names(x))),
                   DR = unlist(dd),
                   age = unlist(aa))
  out = out[order(out$phy, out$sp), ]
}, dd = dr, aa = age, SIMPLIFY = FALSE)
save(dr, age, DR_age, file = paste0(path.out, "DR_age_estimate.RData"))



# estimate the PIC for dr and sp age
botero = read.csv(paste0(path.data, "Botero14.csv"), as.is = TRUE)
# select only birds (ie, drop mammals)
birds.raw = botero %>% 
  filter(TAXON == "birds") %>%
  arrange(SPECIES)
# build the dataset to be used in the analysis
present = names(dr[[1]][[1]]) %in% birds.raw$SPECIES
dat = data.frame("sp" = names(dr[[1]][[1]])[present],
                 #"DR" = dr[[1]][[1]][present],
                 "ssp" = numeric(sum(present)), 
                 "lat" = numeric(sum(present)),
                 "region" = numeric(sum(present)),
                 #"age" = age[[1]][[1]][present],
                 stringsAsFactors = FALSE)
dat = dat[order(dat$sp), ]
# there's 1 sp present in Botero dataset that doesn't appear in the phylogeny "Pseudohirundo_griseopyga"
birds = birds.raw[-which(birds.raw$SPECIES == "Pseudohirundo_griseopyga"), ]
dat$ssp = birds$SUBSPECIES
dat$lat = birds$CENTROID
dat$region = birds$LAT.RANGE
write.csv(dat, row.names = FALSE, file = paste0(path.out, "dat.csv"))
# select only the species with data
absent = trees[[1]][[1]]$tip.label[which(!names(dr[[1]][[1]]) %in% birds.raw$SPECIES)]
pruned = lapply(trees, function(x){
  # there's 1 sp present in Botero dataset that doesn't appear in the phylogeny "Pseudohirundo_griseopyga"
  lapply(x, drop.tip, tip = c(absent, "Pseudohirundo_griseopyga"))
})
dr.pic = mapply(function(dd, pp){
  mapply(function(ddd, ppp){ pic(x = ddd, phy = ppp)}, ddd = dd, ppp = pp, SIMPLIFY = FALSE)
}, dd = dr, pp = trees, SIMPLIFY = FALSE)
age.pic = mapply(function(aa, pp){
  mapply(function(aaa, ppp){ pic(x = aaa, phy = ppp)}, aaa = aa, ppp = pp, SIMPLIFY = FALSE)
}, aa = age, pp = trees, SIMPLIFY = FALSE)
# ssp2pic = dat$ssp; names(ssp2pic) = dat$sp
# ssp.pic = lapply(pruned, function(p) lapply(p, pic, x = ssp2pic))
my.pic = mapply(function(dd, aa, ss){
  out = data.frame(phy = rep(1:length(dd), each = 9993 - 1),
                   DR = unlist(dd),
                   age = unlist(aa))
}, dd = dr.pic, aa = age.pic, ss = ssp.pic, SIMPLIFY = FALSE)
save(my.pic, file = paste0(path.out, "pic.RData"))




# estimate PHYLOGENETIC SIGNAL
ssp2physig = dat$ssp; names(ssp2physig) = dat$sp
ssp.physig = lapply(X = pruned, FUN = function(pp){
  out = mclapply(pp, FUN = function(ppp){ phylosig(tree = ppp, x = ssp2physig, method = "lambda", test = TRUE) }, mc.cores = cores)
  do.call(rbind,out)
})
dr.physig = mapply(function(tt, pp){
  out = mapply(function(ttt, ppp) phylosig(tree = ppp, x = ttt, method = "lambda", test = TRUE), ttt = tt, ppp = pp, SIMPLIFY = FALSE)
  do.call(rbind, out)
}, tt = dr, pp = trees, SIMPLIFY = FALSE)
save(ssp.physig, dr.physig, file = paste0(path.out, "PhySig_results.RData"))



###### FAMILY level
load(file = paste0(path.code, "Reference_species_per_Family.RData"))
dat$Family = character(nrow(dat))
for(x in refBirdLife){
  ind = na.exclude(match(x$Scientific.name, dat$sp))
  if(!is.null(ind) & length(ind) > 0){
    dat$Family[ind] = as.character(x$Family.name[which(!is.na(ind))])
  }
}

####################################
# get the FAMILY phylogeny
nfam = table(Family)[-1]
phyFam = mclapply(pruned, function(pp){
  lapply(pp, function(ppp){
    mapply(function(nf, f){
      aux = sp[Family != f]
      aux = drop.tip(phy = ppp, tip = aux)
      return(aux)
    }, nf = nfam[nfam>2], f = names(nfam[nfam>2]), SIMPLIFY = FALSE)
  })
}, mc.cores = cores)
# Checking relation to Family-level Diversification rates
bdFam = lapply(phyFam, function(x) lapply(x, function(xx) mclapply(xx, birthdeath, mc.cores = cores)))
save(phyFam, bdFam, file = paste0(path.out, "Family_bdLik.RData"))
# estimate ancestral character to use as the phylogenetic mean for the FAMILY
sspFam = lapply(phyFam, function(p) lapply(p, function(pp) mclapply(pp, ace, x = ssp2physig, mc.cores = cores)))
save(sspFam, file = paste0(path.out, "Family_ace.RData"))

















