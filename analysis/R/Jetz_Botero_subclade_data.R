
setwd("~/Desktop/gitHub/protracted_sp/")
library(ggplot2)
library(parallel)
library(geiger)
source("analysis/R/auxiliary_functions.R")


##### Adds taxonomic info to the Botero dataset
# loads the dataset from Botero 2014
load(file = "analysis/data/Botero14.RData")
# creates a discrete categorization of 'ENV.HARSHNESS'
birds$harshness = sapply(birds$ENV.HARSHNESS, bin_harshness, bins = 3)
birds$harshness5 = sapply(birds$ENV.HARSHNESS, bin_harshness, bins = 5)
#   Creates another classification of "latitudinal distrib" based on the centroid.
# Basically, eliminates the "Mixed" category and assigns it either to "Tropical" 
# or "Temperate" depending on where the centroid is located.
birds$LAT.CENTROID = as.character(ifelse(abs(birds$CENTROID) < 23, "Tropical", "Temperate"))

head(birds)
ggplot(birds, aes(ENV.HARSHNESS)) + geom_histogram() + facet_grid(LAT.RANGE ~.)
ggplot(birds, aes(x = abs(CENTROID), y = ENV.HARSHNESS, color = LAT.RANGE)) + geom_point()
ggplot(birds, aes(SUBSPECIES)) + geom_histogram(aes(y = ..density..)) + facet_grid(harshness ~.)
# add taxonomy columns to the dataset
birds$Order = birds$Family = birds$Family.name = birds$Genus = birds$Common.name = NA

# loads list of species per family (each entry a data.frame) from "SSP_DR/R/analysis_Family_level.R"
load(file = "SSP_DR/data/Reference_species_per_Family.RData")
# combine all data.frames into a big one
refBirdLife = do.call(rbind, refBirdLife)
# corrects the names and transforms all factor columns into characters
refBirdLife$Order = sapply(as.character(refBirdLife$Order), correct_name)
refBirdLife$Family.name = as.character(refBirdLife$Family.name)
refBirdLife$Genus = as.character(refBirdLife$Genus)
# assigns the species in Botero's dataset to their corresponding taxonomy (according to Birdlife)
for(i in 1:nrow(refBirdLife)){
  ind = match(refBirdLife$Scientific.name[i], birds$SPECIES)
  if(!is.na(ind)){
    birds[ind, "Common.name"] = refBirdLife$Common.name[i]
    birds[ind, "Genus"] = refBirdLife$Genus[i]
    birds[ind, "Family.name"] = refBirdLife$Family.name[i]
    birds[ind, "Family"] = refBirdLife$Family[i]
    birds[ind, "Order"] = refBirdLife$Order[i]
  }
}
head(birds)
str(birds)

write.csv(birds, file = "analysis/data/Botero14_plus_taxonomy.csv", row.names = FALSE)

birds=read.csv(file = "analysis/data/Botero14_plus_taxonomy.csv", as.is = TRUE)

##### Plots to visualize how the different taxonomic categories are distributed
## Latitude
order = get_prop("Order", type = "latitude")
order2gg = plot_prop(order, minimum = 30, type = "latitude")
ggsave(order2gg, filename = "analysis/output/Botero14_proportion_Lat.range_order.pdf")
family = get_prop("Family.name", type = "latitude")
family2gg = plot_prop(family, minimum = 30, type = "latitude")
ggsave(family2gg, filename = "analysis/output/Botero14_proportion_Lat.range_family.pdf")
genus = get_prop("Genus", type = "latitude")
genus2gg = plot_prop(genus, minimum = 30, type = "latitude")
ggsave(genus2gg, filename = "analysis/output/Botero14_proportion_Lat.range_genus.pdf")

## Harshness
orderH = get_prop("Order", type = "harshness")
order2ggH = plot_prop(orderH, minimum = 30, type = "harshness")
ggsave(order2ggH, filename = "analysis/output/Botero14_proportion_harshness_order.pdf")
familyH = get_prop("Family.name", type = "harshness")
family2ggH = plot_prop(familyH, minimum = 30, type = "harshness")
ggsave(family2ggH, filename = "analysis/output/Botero14_proportion_harshness_family.pdf")
genusH = get_prop("Genus", type = "harshness")
genus2ggH = plot_prop(genusH, minimum = 30, type = "harshness")
ggsave(genus2ggH, filename = "analysis/output/Botero14_proportion_harshness_genus.pdf")





##### Finds all subclades from the phylogenies that meet certain criteria (#tips, proportion of "Mixed", etc.)
jetz = read.tree("analysis/data/Jetz_etal/EricsonStage2_0001_1000/AllBirdsEricson1_reduced.tre")
Ncores = 3
# gets the tips that do NOT have data (appears in Botero14 dataset) and are NOT in islands
tips2drop = lapply(jetz, function(x) which( (!x$tip.label %in% birds$SPECIES) & (birds$ISLAND.DWELLING == 0) ))
# prunes phylogenies accordingly
clean_jetz = mcmapply(phy = jetz, tip = tips2drop, FUN = drop.tip,
                      SIMPLIFY = FALSE, mc.cores = Ncores)

####### LATITUDE
#   Gets the ordered (to match 'phy' tip.label) character vector indicating the
# latitudinal distribution of the species (ie, tipLabels) for each phylogeny.
# [if un-comment, then it's based on theoriginal latitude definition]
tipLatitude = mclapply(clean_jetz, 
                       #function(x) birds$LAT.RANGE[match(x$tip.label, birds$SPECIES)],
                       function(x) birds$LAT.CENTROID[match(x$tip.label, birds$SPECIES)],
                       mc.cores = Ncores)
# finally, gets the subclades
these_subclades = mcmapply(FUN = get_matching_clades,
                           phy = clean_jetz,
                           tipLabels = tipLatitude,
                           MoreArgs = list(clade.size = 30, minimum = TRUE,
                                           mixed = 0.3, multiLat = 0),
                           SIMPLIFY = FALSE, mc.cores = Ncores)
lapply(these_subclades, sapply, length)
# now extracts the actual subclades from the Jetz phylogeny
subclades = mcmapply(FUN = function(ppp, tsub){
  foo2 = sapply(tsub, FUN = getMRCA, phy = ppp)
  lapply(foo2, FUN = extract.clade, phy = ppp)
}, ppp = clean_jetz, tsub = these_subclades,
SIMPLIFY = FALSE, mc.cores = Ncores)
# tests if the extraction was sucessfull
all.equal(lapply(subclades, sapply, Ntip), lapply(these_subclades, sapply, length))
# transforms the phylogenies into vector of branching times
branches = lapply(subclades, lapply, branching.times)
# creates another set of branching times scaled to 1
branches.unity = lapply(branches, lapply, function(x) x / max(x))


####### HARSHNESS
#   Gets the ordered (to match 'phy' tip.label) character vector indicating the
# harshness category of the species (ie, tipLabels) for each phylogeny.
tipHarshness = mclapply(clean_jetz, 
                        function(x) birds$harshness[match(x$tip.label, birds$SPECIES)],
                        mc.cores = Ncores)
# finally, gets the subclades
these_subcladesHarsh = mcmapply(FUN = get_matching_cladesHarsh,
                                phy = clean_jetz,
                                tipLabels = tipHarshness,
                                MoreArgs = list(clade.size = 30, minimum = TRUE,
                                                mixed = 0.3),
                                SIMPLIFY = FALSE, mc.cores = Ncores)
lapply(these_subcladesHarsh, sapply, length)
# now extracts the actual subclades from the Jetz phylogeny
subcladesHarsh = mcmapply(FUN = function(ppp, tsub){
  foo2 = sapply(tsub, FUN = getMRCA, phy = ppp)
  lapply(foo2, FUN = extract.clade, phy = ppp)
}, ppp = clean_jetz, tsub = these_subcladesHarsh,
SIMPLIFY = FALSE, mc.cores = Ncores)
# tests if the extraction was sucessfull
all.equal(lapply(subcladesHarsh, sapply, Ntip), lapply(these_subcladesHarsh, sapply, length))
# transforms the phylogenies into vector of branching times
branchesHarsh = lapply(subcladesHarsh, lapply, branching.times)
# creates another set of branching times scaled to 1
branchesHarsh.unity = lapply(branchesHarsh, lapply, function(x) x / max(x))


####### "clean" HARSHNESS
# keep only the clades with only one category of harshness
# drops the "mixed" species from the harshness subclades
subcladesHarshClean = mcmapply(FUN = function(ppp, tsub){
  harsh = lapply(tsub, attr, which = "harshness")
  bool = mapply(FUN = function(atr, tip) tip[which(atr != max(atr))],
                atr = harsh, tip = tsub, SIMPLIFY = FALSE)
  mapply(FUN = drop.tip, phy = ppp, tip = bool, SIMPLIFY = FALSE)
}, ppp = subcladesHarsh, tsub = these_subcladesHarsh,
SIMPLIFY = FALSE, mc.cores = Ncores)
sapply(subcladesHarshClean, length) # number of phylogenies
sapply(subcladesHarshClean, sapply, Ntip) # number of tips per phylogeny
# transforms the phylogenies into vector of branching times
branchesHarshClean = lapply(subcladesHarshClean, lapply, branching.times)
# creates another set of branching times scaled to 1
branchesHarshClean.unity = lapply(branchesHarshClean, lapply, function(x) x / max(x))


# SAVES
save(these_subclades, subclades, these_subcladesHarsh, subcladesHarsh,
     file = "analysis/data/Jetz_Botero14_30min_0.3mixed_0multilat_subclades.RData")
save(branches, branches.unity, # temperate
     branchesHarsh, branchesHarsh.unity, # harshness
     branchesHarshClean, branchesHarshClean.unity, # "pure" harshness
     file = "analysis/data/Jetz_Botero14_30min_0.3mixed_0multilat_branches.RData")



# Let's check the proportion of Temp in the Trop clades (and vice-versa)
foo = lapply(these_subclades, sapply, table_prop_phylo, categories = c("Temperate", "Tropical"))
propMixed = data.frame(do.call(rbind, lapply(foo, t)))
cntrl = mapply(x = 1:10, times = sapply(these_subclades, length), FUN = rep)
propMixed$phy = unlist(cntrl)
propMixed$N = unlist(sapply(these_subclades, sapply, length))
head(propMixed)
propMixed.melt = reshape2::melt(propMixed, id.vars = 3:4)
head(propMixed.melt)
write.csv(propMixed, file = "analysis/data/Jetz_Botero14_propMixed.csv", row.names = FALSE)

ggplot(propMixed.melt, aes(x = phy, y = N, group = phy, color = variable)) + geom_jitter(height = 0)
propMixed.melt %>% filter(variable == "Temperate") %>% ggplot(aes(x = N, y = value, group = phy)) + geom_jitter(height = 0)
propMixed.melt %>% filter(phy == 2) %>% ggplot(aes(x = N, y = value, group = phy)) + geom_jitter(height = 0)+facet_grid(variable ~ .)

# checks for possible AGE x SpeciesNumber correlation
get_age = function(x) max(x$edge.length)
get_harshness = function(x) mean(birds$ENV.HARSHNESS[match(x$tip.label, birds$SPECIES)])
testSubclades = data.frame(tips = unlist(sapply(subclades, sapply, Ntip)),
                           age = unlist(sapply(subclades, sapply, get_age)),
                           phy = unlist(cntrl),
                           harshness = unlist(sapply(subclades, sapply, get_harshness)))
ggplot(testSubclades, aes(x = age, y = tips)) + geom_point() + geom_smooth(method='lm')
ggplot(testSubclades, aes(x = harshness, y = tips)) + geom_point() + geom_smooth(method='lm')



# gets number of subspecies per species for each subclade
spp = lapply(subclades, sapply, get_Botero_spp)
# uses the number of subspecies per species to create an exponencial prior for b1
priorb1exp = lapply(spp, lapply, get_prior_exp, output = "function", log = TRUE)
save(priorb1exp, file = "analysis/data/Jetz_Botero14_subclades_priorb1exp.RData")




# Let's check the proportion of Stable and Regular in the Variable clades (and vice-versa)
fooH = lapply(these_subcladesHarsh, sapply, table_prop_phylo, categories = c("Stable", "Regular", "Variable"))
propMixedHarsh = data.frame(do.call(rbind, lapply(fooH, t)))
cntrlH = mapply(x = 1:10, times = sapply(these_subcladesHarsh, length), FUN = rep)
propMixedHarsh$phy = unlist(cntrlH)
propMixedHarsh$N = unlist(sapply(these_subcladesHarsh, sapply, length))
propMixedHarsh$estimate = rownames(propMixedHarsh)
propMixedHarsh$max = apply(propMixedHarsh[ , 1:3], 1, max)
propMixedHarsh$which.max = factor(x = apply(propMixedHarsh[ , 1:3], 1, which.max),
                                  labels = c(colnames(propMixedHarsh))[1:3])
head(propMixedHarsh)
propMixed.meltHarsh = reshape2::melt(propMixedHarsh, id.vars = 4:8)
head(propMixed.meltHarsh)
write.csv(propMixedHarsh, file = "analysis/data/Jetz_Botero14_propMixedHarsh.csv", row.names = FALSE)

ggplot(propMixed.meltHarsh, aes(x = phy, y = N, group = phy, color = variable)) +
  geom_jitter(height = 0)
propMixed.meltHarsh %>%
  filter(phy == 1) %>%
  ggplot(aes(x = estimate, y = value, fill = variable)) +
  geom_bar(stat = "identity")
ggplot(propMixedHarsh, aes(x = N, y = max)) +
  geom_point() + facet_grid(which.max ~ .)


# gets number of subspecies per species for each subclade
sppH = lapply(subcladesHarsh, sapply, get_Botero_spp)
# uses the number of subspecies per species to create an exponencial prior for b1
priorb1expHarsh = lapply(sppH, lapply, get_prior_exp, output = "function", log = TRUE)
save(priorb1expHarsh, file = "analysis/data/Jetz_Botero14_subclades_priorb1expHarsh.RData")




# gets number of subspecies per species for each subclade
sppHclean = lapply(subcladesHarshClean, sapply, get_Botero_spp)
# uses the number of subspecies per species to create an exponencial prior for b1
priorb1expHarshClean = lapply(sppHclean, lapply, get_prior_exp, output = "function", log = TRUE)
save(priorb1expHarshClean, file = "analysis/data/Jetz_Botero14_subclades_priorb1expHarshClean.RData")
