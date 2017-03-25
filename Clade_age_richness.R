


library(ape)
library(rgbif)
#library(xlsx) # to read Excel files into R


aves = read.csv("protracted_sp/analysis/data/BirdLife_Checklist_Version_90/BirdLife_Checklist_Version_9_reduced.csv", header = TRUE, as.is = TRUE)
aves$Order = as.factor(aves$Order)
aves$Family.name = as.factor(aves$Family.name)
head(aves)
str(aves)
nrow(aves)


ref = read.csv("protracted_sp/analysis/data/Clade_age_richness.csv", header = TRUE, as.is = TRUE)
head(ref)
clades = strsplit(ref$Clade, split = " \\+ ")
ll = sapply(clades, length)
clades = data.frame("taxa" = unlist(clades),
                    "contr" = rep(1:length(ll), ll))
gbif = lapply(clades$taxa, function(x) occ_search(scientificName = x))
gbifsp = lapply(gbif, function(x) if(is.null(x$data)){NA}else{unique(x$data$name)})
sp = list()
for(i in 1:length(ll)){
  ind = which(clades$contr == i)
  if(ll[i] == 1){
    sp[[i]] = gbifsp[[ind]]
  } else{
    sp[[i]] = na.exclude(unlist(gbifsp[ind]))
  }
}
bool = sapply(sp, function(x) class(x) == "character")
sp[bool] = lapply(sp[bool], gsub, pattern = " ", replacement = "_")

trees = read.tree("protracted_sp/analysis/data/Jetz_etal/EricsonStage2_0001_1000/AllBirdsEricson1_reduced.tre")
jetz = trees[[1]]

mrca = mapply(function(vec, contr){
  if(contr){
    aaa = getMRCA(phy = jetz, tip = vec)
    return(ifelse(is.null(aaa), NA, aaa))
  } else{
    return(NA)
  }}, vec = sp, contr = bool)
phy = lapply(na.exclude(mrca), extract.clade, phy = jetz)
ref$ntip = numeric(length(clades))
ref$ntip[bool] = sapply(phy, Ntip)
ref
