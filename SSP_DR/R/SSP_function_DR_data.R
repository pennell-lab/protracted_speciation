


# setwd("/Users/mauro/Desktop/gitHub")

library(parallel)
library(doParallel)
library(reshape2)
library(ape)
library(dplyr)
library(MASS)
library(phytools)
library(PBD)



# get 1 of Jetz1s trees
trees = read.tree("protracted_sp/analysis/data/Jetz_etal/EricsonStage2_0001_1000/AllBirdsEricson1_reduced.tre")
jetz = trees[[1]]


# estimate the DR
# source("protracted_sp/SSP_DR/R/auxiliary_functions.R")
# dr = get_DR(jetz)
# age = get_sp_age(jetz)
# save(dr, age, file = "protracted_sp/SSP_DR/DR_estimate.RData")
load("protracted_sp/SSP_DR/data/DR_estimate.RData")


botero = read.csv("protracted_sp/analysis/data/Botero14.csv", as.is = TRUE)
# select only birds (ie, drop mammals)
birds = botero %>% 
  filter(TAXON == "birds") %>%
  arrange(SPECIES)
summary(birds)



# select only the species with data
present = names(dr) %in% birds$SPECIES
dat = data.frame("sp" = names(dr)[present],
                 "DR" = dr[present],
                 "ssp" = numeric(sum(present)),
                 "lat" = numeric(sum(present)),
                 "region" = numeric(sum(present)),
                 "age" = age[present],
                 stringsAsFactors = FALSE)
dat = dat[sort(dat$sp), ]
# there's 1 sp present in Botero dataset that doesn't appear in the phylogeny "Pseudohirundo_griseopyga"
birds = birds[-which(birds$SPECIES == "Pseudohirundo_griseopyga"), ]
dat$ssp = birds$SUBSPECIES
dat$lat = birds$CENTROID
dat$region = birds$LAT.RANGE
dat$genus = sapply(dat$sp, function(x) strsplit(x, split = "_")[[1]][1])
load(file = "protracted_sp/SSP_DR/data/Reference_species_per_Family.RData")
dat$Order = dat$Family = dat$Genus = character(nrow(dat))
for(x in refBirdLife){
  ind0 = match(x$Scientific.name, dat$sp)
  ind = na.exclude(ind0)
  if(!is.null(ind) & length(ind) > 0){
    dat$Order[ind] = as.character(x$Order[which(!is.na(ind0))])
    dat$Family[ind] = as.character(x$Family.name[which(!is.na(ind0))])
    dat$Genus[ind] = as.character(x$Genus[which(!is.na(ind0))])
  }
}
write.csv(dat, row.names = FALSE, file = "protracted_sp/SSP_DR/data/DR_ssp_centroid.csv")
# library(rgbif)
# gg = dat %>% filter(Genus == "") %>% group_by(genus)
# missg = unique(gg$genus)
# gbif = lapply(missg, function(x) name_lookup(query=x, rank="genus", return="data"))
# gbifAves = lapply(gbif, function(x){x[x$class == "Aves", c(1,2,4:11,16,17,18,21:24,28,29)]})
# gbifRep = lapply(gbifAves, function(x) c(na.exclude(unique(x$family))))
# gbifFail = sapply(gbifRep, is.null)
# missing = missg[gbifFail]
# gbifWin = sapply(gbifRep, function(x)length(x)==1)
# found = missg[gbifWin]
# dat2 = dat
# for(i in 1:length(found)){
#   dat2$Family[which(dat2$Family == "" & dat2$genus == found[i])] = unlist(gbifRep[gbifWin])[i]
# }
# missg2 = apply(dat2[dat2$Family == "", ], MARGIN = 1, function(x) paste(strsplit(x[1], split = "_")[[1]], collapse = " "))
# gbif2 = list()
# #lapply(missg2, function(x) name_lookup(query=x, rank="species", return="data", config = httr::verbose()))
# for(i in 1:length(missg2)){ 
#   #478: Parus davidi  --> Paridae
#   #674: Vanellus malarbaricus  --> Charadriidae
#   gbif2[[i]] = name_lookup(query=missg2[i], rank="species", return="data", status = "ACCEPTED", config = httr::verbose())
# }
# bif = lapply(gbif2, function(x)if(is.null(x)){NULL}else{data.frame(x[ , c(1,2,4:11,16,17,18,21:24)])})
# fff = lapply(bif, function(x) if(is.null(x)){NA}else{na.exclude(unique(x$family))})
# gbifWin2 = sapply(fff, function(x) if(is.null(x) | length(x) > 1){FALSE}else{TRUE})
# found2 = sapply(missg2[gbifWin2], function(x) paste0(strsplit(x, " ")[[1]], collapse = "_"))
# for(i in 1:length(found2)){
#   dat2$Family[which(dat2$Family == "" & dat2$sp == found2[i])] = unlist(fff[gbifWin2])[i]
# }
# for(i in unique(dat2$Family)){
#   ind = which(dat2$Family == i)
#   if(any(dat2$Order[ind] != "")){
#     oo = unique(dat2$Order[ind])[unique(dat2$Order[ind]) != ""]
#     dat2$Order[ind] = gsub(" ", "", oo)
#   }
# }
# load("protracted_sp/SSP_DR/data/DR_age_estimate.RData")
# dr.age.200 = do.call(rbind, DR_age) %>% group_by(sp) %>%
#   summarize(dr = mean(log(DR)),
#             age = mean(age))
# for(i in 1:nrow(dr.age.200)){
#   ind = which(dat$sp == dr.age.200$sp[i])
#   dat$DR[ind] = dr.age.200$dr[i]
#   dat$age[ind] = dr.age.200$age[i]
# }

# load("protracted_sp/SSP_DR/data/DR_age_estimate.RData")
# str(DR_age$Ericson)
# #ddd = DR_age$Ericson %>% group_by(sp) %>% summarise(mean = mean(dr), log = mean(log(dr)))
# ddd = data.frame(mean = NA, log = NA)
# for(i in 1:9993){
#   id = levels(DR_age$Ericson$sp)[i]
#   rrr = c(DR_age$Ericson$DR[which(DR_age$Ericson$sp == id)], DR_age$Hackett$DR[which(DR_age$Hackett$sp == id)])
#   ddd[i, ] = c(mean(rrr), mean(log(rrr)))
# }
# dat$DR[match(levels(DR_age$Ericson$sp), dat$sp)] = ddd$log

### write.csv(dat2, "protracted_sp/SSP_DR/data/Updated_Botero_dat.csv")
dat=read.csv(file = "protracted_sp/SSP_DR/data/Complete_Botero.csv", row.names = 1, as.is = TRUE)


#### Now considering the phylogenetic dependency of the observations
# pruned = drop.tip(jetz, tip = c(jetz$tip.label[!present], "Pseudohirundo_griseopyga"))
# save(pruned, file = "protracted_sp/SSP_DR/data/pruned.RData")
load(file = "protracted_sp/SSP_DR/data/pruned.RData")
my.pic = data.frame("DR" = pic(dat$DR, pruned),
                    "ssp" = pic(dat$ssp, pruned))
summary(my.pic)
# estimate PHYLOGENETIC SIGNAL
# tophysig = list(dat$DR, dat$ssp)
# tophysig = lapply(tophysig, function(x){names(x) = dat$sp;x})
# physig = mclapply(X = tophysig, FUN = phylosig, tree = pruned, method = "lambda", mc.cores = 2)
load("protracted_sp/analysis/data/physig_results.RData")
physig



###### DR vs. SSP family level
datsummary = dat %>%
  group_by(Family) %>%
  summarise(Nsp = length(sp),
            Nssp = sum(ssp, na.rm = TRUE),
            ssp.mean = mean(ssp, na.rm = TRUE),
            ssp.median = median(ssp, na.rm = TRUE),
            ssp.sd = sd(ssp, na.rm = TRUE),
            dr.mean = mean(DR, na.rm = TRUE),
            dr.median = median(DR, na.rm = TRUE),
            dr.sd = sd(DR, na.rm = TRUE),
            age.mean = mean(age, na.rm = TRUE),
            age.median = median(age, na.rm = TRUE),
            age.sd = sd(age, na.rm = TRUE),
            centroid = mean(lat, na.rm = TRUE),
            centroid.median = median(lat, na.rm = TRUE),
            centroid.sd = sd(lat, na.rm = TRUE))
datsummary$ssp.sd[is.na(datsummary$ssp.sd)] = 0
datsummary$dr.sd[is.na(datsummary$dr.sd)] = 0
datsummary$age.sd[is.na(datsummary$age.sd)] = 0
datsummary = datsummary[order(datsummary$Nsp, decreasing = T), ]
datsummary
####################################
# Checking relation to Family-level Diversification rates
# nfam = table(dat$Family)
# tipmin = 6
# phyFam = mcmapply(FUN = function(nf, f){
#   if(nf > tipmin){
#     aux = dat$sp[dat$Family != f]
#     aux = drop.tip(phy = jetz, tip = aux)
#     return(aux)
#   } else {
#     return(NA)
#   }
# }, nf = nfam, f = names(nfam), SIMPLIFY = FALSE, mc.cores = 3)
# bdFam = lapply(phyFam, function(x) if(is.na(x)){NA}else{birthdeath(x)})
####################################
# Checking relation to Protracted BD rates
# brnchsFam = lapply(phyFam, function(x) if(is.na(x)){NA}else{branching.times(x)})
# PBDfam = lapply(brnchsFam, function(x) if(is.na(x)){NA}else{pbd_ML(x, initparsopt = c(0.2,0.1,1,0.1), exteq = 0, btorph = 0)})
####
# save(PBDfam, phyFam, bdFam, file = "protracted_sp/SSP_DR/data/Family.RData")
load(file = "protracted_sp/SSP_DR/data/Family.RData")
# the SPECIATION rate fit with ape::birthdeath, using the method of Nee etal 1994
datsummary$b = NA
# the DIVERSIFICATION rate fit with ape::birthdeath, using the method of Nee etal 1994
datsummary$bdLik = NA
# the SPECIATION rate fit with PBD::pbd_ML
datsummary$pbd.specia = NA
# the SPLIT rate fit with PBD::pbd_ML
datsummary$pbd.b = NA
# the TIME TO SPECIATION rate fit with PBD::pbd_ML
datsummary$pbd.la1 = NA
# the EXTINCTION rate fit with PBD::pbd_ML
datsummary$pbd.mu = NA
# the EXTINCTION2 rate fit with PBD::pbd_ML
datsummary$pbd.mu2 = NA
for(i in 1:length(bdFam)){
  ind = which(datsummary$Family == names(bdFam)[i])
  if(!is.na(bdFam[[i]]))
    if(all(c(PBDfam[[i]]$b, PBDfam[[i]]$lambda_1, PBDfam[[i]]$mu_1, PBDfam[[i]]$mu_2) < 100)){
      datsummary$b[ind] = bdFam[[i]][[4]][2] / (1 - bdFam[[i]][[4]][1])
      datsummary$bdLik[ind] = bdFam[[i]][[4]][2]
      datsummary$pbd.specia[ind] = (PBDfam[[i]]$b - PBDfam[[i]]$mu_2) / PBDfam[[i]]$lambda_1
      datsummary$pbd.b[ind] = PBDfam[[i]]$b
      datsummary$pbd.la1[ind] = PBDfam[[i]]$lambda_1
      datsummary$pbd.mu[ind] = PBDfam[[i]]$mu_1
      datsummary$pbd.mu2[ind] = PBDfam[[i]]$mu_2
    }
}
summary(datsummary)

####################################
# using the model from Phillimore 2010 to estimate SSP origination and extinction
source("protracted_sp/SSP_DR/R/auxiliary_functions.R")
library(stats4)
phill = lapply(datsummary$Family, function(xxx){
  aux = na.exclude(dat[which(dat$Family == xxx), ])
  out = NA
  if(nrow(aux) > 3){
    fff = phillimore_model(n = aux$ssp, t = aux$age)
    out = try(mle(fff,
                  start = list(rate = 0.1, epsilon = 0.2), 
                  method = "L-BFGS-B",
                  lower = c(1e-10, 1e-10), upper = c(Inf, 0.99999)))
  }
  return(out)
})
phill.res = do.call(rbind, lapply(phill, function(x) if(class(x) == "mle"){coef(x)}else(c(NA,NA))))
datsummary$ssp.rate = phill.res[ , 1]
datsummary$ssp.epsilon = phill.res[ , 2]
datsummary$ssp.lambda = phill.res[ , 1] / (1 - phill.res[ , 2])
summary(datsummary)
save(datsummary, file = "protracted_sp/SSP_DR/data/datsummary.RData")



####################################
# check if richness is nested
# SSP per species within GENUS
datgenus = dat %>%
  group_by(genus) %>%
  summarise(Nsp = length(sp),
            Nssp = sum(ssp, na.rm = TRUE),
            ssp.mean = mean(ssp, na.rm = TRUE),
            ssp.median = median(ssp, na.rm = TRUE),
            ssp.sd = sd(ssp, na.rm = TRUE),
            dr.mean = mean(DR, na.rm = TRUE),
            dr.median = median(DR, na.rm = TRUE),
            dr.sd = sd(DR, na.rm = TRUE),
            age.mean = mean(age, na.rm = TRUE),
            age.median = median(age, na.rm = TRUE),
            age.sd = sd(age, na.rm = TRUE),
            centroid = mean(lat, na.rm = TRUE),
            centroid.median = median(lat, na.rm = TRUE),
            centroid.sd = sd(lat, na.rm = TRUE))
datgenus = datgenus[order(datgenus$Nsp, decreasing = T), ]
datgenus
save(datgenus, file = "protracted_sp/SSP_DR/data/datgenus.RData")
# SP per GENERA within FAMILY
datfamily = dat %>%
  group_by(Family) %>%
  summarise(Ngen = length(unique(genus)),
            Nsp = length(sp),
            sp.mean = mean(table(genus)),
            sp.median = median(table(genus)),
            sp.sd = sd(table(genus)),
            ssp.mean = mean(ssp, na.rm = TRUE),
            ssp.median = median(ssp, na.rm = TRUE),
            ssp.sd = sd(ssp, na.rm = TRUE),
            dr.mean = mean(DR, na.rm = TRUE),
            dr.median = median(DR, na.rm = TRUE),
            dr.sd = sd(DR, na.rm = TRUE),
            age.mean = mean(age, na.rm = TRUE),
            age.median = median(age, na.rm = TRUE),
            age.sd = sd(age, na.rm = TRUE),
            centroid = mean(lat, na.rm = TRUE),
            centroid.median = median(lat, na.rm = TRUE),
            centroid.sd = sd(lat, na.rm = TRUE))
datfamily = datfamily[order(datfamily$Ngen, decreasing = T), ]
datfamily
save(datfamily, file = "protracted_sp/SSP_DR/data/datfamily.RData")





#########################
# counts the number of sp in each region per FAMILY
rrr = dat %>% group_by(Family) %>% 
  summarise(Trop = round(mean(region == "Tropical", na.rm = TRUE), 3)*100,
            Temp = round(mean(region == "Temperate", na.rm = TRUE), 3)*100,
            Mixed = round(mean(region == "Mixed", na.rm = TRUE), 3)*100,
            known = sum(region %in% c("Tropical", "Temperate", "Mixed")),
            N = length(region))
print.data.frame(rrr)





###########################
# COMPLETE table adding taxonomy for the Jetz' sp not present in Botero's dataset
#
#
# trees = read.tree("protracted_sp/analysis/data/Jetz_etal/EricsonStage2_0001_1000/AllBirdsEricson1_reduced.tre")
# jetz = trees[[1]]
# dat.new = read.csv("protracted_sp/SSP_DR/data/Updated_Botero_dat.csv", as.is = T)
# miss = jetz$tip.label[which(!jetz$tip.label %in% dat.new$sp)]
# m2gbif = gsub(pattern = "_", replacement = " ", x = miss)
# gbif3 = list()
# for(i in 3234:length(m2gbif)){gbif3[[i]] = name_lookup(query=m2gbif[i], rank="species", return="data", status = "ACCEPTED", config = httr::verbose())}
# # "Updated_Boreto_aux.txt" --> aa
# found3 = data.frame(sp = "", genus = "", Family = "", Order = "", stringsAsFactors = F)
# prob = list();aux=1;
# nn = c()
# for(i in 1:length(gbif3)){
#   if(i %in% aa$id){
#     iii = which(aa$id == i)
#     found3 = rbind(found3, c(aa$sp[iii], gsub(" +[a-z]*", "", aa$sp[iii]), aa$Family[iii], aa$Order[iii]))
#   } else{
#     fc = na.exclude(unique(gbif3[[i]]$family))
#     if(length(fc) == 1){
#       found3 = rbind(found3, c(m2gbif[i], gsub(" +[a-z]*", "", m2gbif[i]), fc, ""))
#     } else{
#       nn[aux] = m2gbif[i]
#       prob[[aux]] = gbif3[[i]];aux=aux+1
#     }
#   }
# }
# found3 = found3[-1, ]
# gg = gsub(" +[a-z]*", "", nn)
# cntr = c()
# for(i in 1:length(gg)){
#   if(gg[i] %in% found3$genus){
#     www = found3$Family[which(found3$genus == gg[i])]
#     found3 = rbind(found3, c(nn[i], gg[i], www[1], ""))
#     cat(i, "\t")
#     cntr = c(cntr, i)
#   } 
# }
# nn = nn[-cntr]
# gg = gg[-cntr]
# prob = prob[-cntr]
# # "Complete_Botero_aux.txt" --> bb
# for(i in 1:nrow(bb)){
#   ind = which(gg == bb$genus[i])
#   for(kkk in ind){
#     found3 = rbind(found3, c(nn[kkk], bb$genus[i], bb$family[i], bb$order[i]))
#   }
# }
# found2dat.new = dat.newa.frame(sp = found3$sp,
#                                DR = NA,
#                                ssp = NA,
#                                lat = NA,
#                                region = "NA",
#                                age = NA,
#                                genus = found3$genus,
#                                Genus = "",
#                                Family = found3$Family,
#                                Order = found3$Order,
#                                stringsAsFactors = FALSE)
# dat.new = rbind(dat.new, found2dat.new)
# 
# elton = read.csv("/Users/mauro/Dropbox/Masters_project/EltonTraits/BirdFuncDat.txt", sep = "\t", as.is = T)
# dat.new$eltonF = dat.new$eltonO = ""
# for(i in 1:nrow(elton)){
#   ind = which(dat.new$sp == gsub(" ","_",elton$Scientific[i]))
#   if(length(ind) == 1){
#     dat.new$eltonF[ind] = elton$BLFamilyLatin[i]
#     dat.new$eltonO[ind] = elton$IOCOrder[i]
#   }
# }
# summary(dat.new)
# 
# for(i in unique(dat.new$Family)){
#   ind = which(dat.new$Family == i)
#   if(any(dat.new$Order[ind] != "")){
#     oo = unique(dat.new$Order[ind])[unique(dat.new$Order[ind]) != ""]
#     dat.new$Order[ind] = gsub(" ", "", oo)
#   }
# }
# write.csv(dat.new, file = "protracted_sp/SSP_DR/data/Complete_Botero.csv")











#### AVIBASE taxonomic reference
# https://avibase.bsc-eoc.org/avibase.jsp?lang=EN&pg=families
# classification follows Howard and Moore for extant families
foo = read.csv(url("https://avibase.bsc-eoc.org/avibase.jsp?lang=EN&pg=families"))
ind.fam = grep(x = foo[,1], pattern = "<a href=avibase.jsp?pg=search&fam=", fixed = TRUE)
ind.ord = grep(x = foo[,1], pattern = "<P><b>", fixed = TRUE)[-1]
avi.fam = sapply(ind.fam, function(x){gsub(pattern = ".*lang=EN>", replacement = "", x = gsub(pattern = "</a> <br>", replacement = "", x = foo[x, 1]))})
avi.fam = sapply(avi.fam, function(xxx){gsub(pattern = "</a>  &#134;<br>", replacement = "", x = xxx)})
avi.order = sapply(ind.ord, function(x){gsub(pattern = ".*<P><b>", replacement = "", x = gsub(pattern = "</b><br>", replacement = "", x = foo[x, 1]))})
avibase = data.frame(family = avi.fam, order = "Passeriformes", stringsAsFactors = FALSE)
for(i in 2:length(ind.ord)){
  avibase[which(ind.fam > ind.ord[i-1] & ind.fam < ind.ord[i]), "order"] = avi.order[i-1]
}
unique(dat$Order)[which(!unique(dat$Order) %in% avibase$order)]
avibase[which(avibase$family=="CiconSuliformesiidae"),]
# ORDER.now       ORDER.Avibase         FAMILY                            Comment
# Ciconiiformes   Pelecaniformes        "Ciconiidae"
# Suliformes      Pelecaniformes        "Sulidae", "Fregatidae", "Phalacrocoracidae", "Anhingidae", "Plotopteridae"            Changed recently, making Pelecaniformes polyphyletic
# Cathartiformes    "Teratornithidae", "Cathartidae"
# Pterocliformes




