

# setwd("~/Desktop/gitHub/")


library(ape)
library(dplyr)
library(doParallel)


# taxonomic list from Birdlife
aves = read.csv("protracted_sp/analysis/data/BirdLife_Checklist_Version_90/BirdLife_Checklist_Version_9_reduced.csv", header = TRUE, as.is = TRUE)
aves$Order = as.factor(aves$Order)
aves$Family.name = as.factor(aves$Family.name)
aux = sapply(aves$Scientific.name, strsplit, split = " ")
aves$Scientific.name = sapply(aux, paste0, collapse = "_")
aves$Genus = as.factor(sapply(aux, "[", 1))
head(aves)
str(aves)
nrow(aves)

# get the reference SPECIES per family
refBirdLife = lapply(levels(aves$Family.name), function(x){
  out = aves %>% 
    filter(Family.name == x) %>% 
    filter(!duplicated(Scientific.name)) %>% 
    select(Order, Family.name, Family, Genus, Common.name, Scientific.name)
  out
})
save(refBirdLife, file = "SSP_DR/data/Reference_species_per_Family.RData")
sapply(refBirdLife, nrow) # there's less because I removed duplicated entries
as.numeric(table((aves$Family.name)))


# sample of 1.000 Jetz trees
trees = read.tree("protracted_sp/analysis/data/Jetz_etal/EricsonStage2_0001_1000/AllBirdsEricson1_reduced.tre")
jetz = trees[[1]]
found = lapply(refBirdLife, function(x){
  na.exclude(match(x$Scientific.name, table = jetz$tip.label))
})
# TOTAL number of species found
sum(sapply(found, length)) # 8373
Ntip(jetz) # 9993  -> diff=1620
# total NUMBER of species found per FAMILY
sapply(found, length)
# PROPORTION of species found per FAMILY
round(sapply(found, length)/sapply(refBirdLife, nrow), 3)




# EXTRACT the phylogeny for each family
bool = sapply(found, length) > 2
nodes = sapply(found[bool], getMRCA, phy = jetz)
phylos = mclapply(nodes, extract.clade, phy = jetz, mc.cores = 3)
# tips per family phylogeny
sapply(phylos, Ntip)
# number of species per family according to BirdLife
sapply(refBirdLife[bool], nrow)
# compare NTIP and number of SPECIES
sapply(phylos, Ntip) <= sapply(refBirdLife[bool], nrow)
# compare NTIP and number of species FOUND
(bool2 = sapply(phylos, Ntip) <= sapply(found[bool], length))
sapply(phylos[bool2], Ntip)









