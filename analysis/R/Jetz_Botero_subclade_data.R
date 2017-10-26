
setwd("~/Desktop/gitHub/protracted_sp/")
library(ggplot2)


##### Adds taxonomic info to the Botero dataset
# loads the dataset from Botero 2014
load(file = "analysis/data/Botero14.RData")
# creates a discrete categorization of 'ENV.HARSHNESS'
bin_harshness = function(x){
  if(x < -1.5){return("ExtStable")}
  if(x < -0.5){return("Stable")}
  if(x < 0.5){return("Regular")}
  if(x < 1.5){return("Variable")}
  return("ExtVariable")
}
birds$harshness = sapply(birds$ENV.HARSHNESS, bin_harshness)
head(birds)
ggplot(birds, aes(ENV.HARSHNESS)) + geom_histogram() + facet_grid(LAT.RANGE ~.)
ggplot(birds, aes(x = abs(CENTROID), y = ENV.HARSHNESS, color = LAT.RANGE)) + geom_point()
# add taxonomy columns to the dataset
birds$Order = birds$Family = birds$Family.name = birds$Genus = birds$Common.name = NA

# loads list of species per family (each entry a data.frame) from "SSP_DR/R/analysis_Family_level.R"
load(file = "SSP_DR/data/Reference_species_per_Family.RData")
# combine all data.frames into a big one
refBirdLife = do.call(rbind, refBirdLife)
correct_name = function(vec){
  # transforms the given names first letter into uppercase and all other into lowercase
  # just because.....
  
  out = strsplit(vec, split = "")[[1]]
  paste0(toupper(out[1]), paste0(tolower(out[-1]), collapse = ""))
}
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

save(birds, file = "analysis/data/Botero14_plus_taxonomy.RData")




##### Functions to analyse how the different taxonomic categories are distributed
get_prop_latitude = function(column, ref = birds){
  # estimates the number of species that are "Temperate", "Tropical", or "Mixed"
  #   - column, numeric or character indicating the column of the taxonomy level to be used
  #   - ref, data.frame with the information about latitudinal distribution. MUST have a column named "LAT.RANGE"
  
  require(dplyr)
  ref = ref[!is.na(ref[[column]]), ]
  groups = unique(ref[[column]])
  
  out = data.frame(matrix(data = NA, nrow = length(groups), ncol = 4,
                          dimnames = list(groups, c("Temperate", "Tropical", "Mixed", "N"))))
  
  if(class(column) == "numeric"){
    column = colnames(ref)[column]
  }
  
  for(i in 1:length(groups)){
    this_row = ref %>% 
      filter_(paste(column, "==", shQuote(groups[i]))) %>%
      summarise(Temp = sum(LAT.RANGE == "Temperate"),
                Trop = sum(LAT.RANGE == "Tropical"),
                Mixed = sum(LAT.RANGE == "Mixed"),
                N = length(LAT.RANGE))
    if(sum(this_row[1:3]) == this_row[4]){
      out[i, ] = as.numeric(this_row)
    } else{
      warning(paste0("Problems talling the latitudinal distribution for group ", groups[i]))
    }
  }
  
  return(out)
}
plot_prop_latitude = function(out, minimum = 1){
  # creates a stacked bar plot from the output of "get_prop_latitude"
  #   - minimum, is the minimum number of species a group must have to be plotted
  require(ggplot2)
  require(reshape2)
  out$taxa = rownames(out)
  
  gg = out %>% 
    filter(N >= minimum) %>%
    mutate(Temperate = 100 * (Temperate / N),
           Tropical = 100 * (Tropical / N),
           Mixed = 100  * (Mixed / N) ) %>%
    arrange(desc(N)) %>%
    mutate(taxa = factor(taxa, levels = taxa)) %>%
    melt(id.var = 4:5) %>%
    ggplot(aes(x = taxa, y = value, fill = variable)) + 
    geom_bar(stat = "identity") + 
    geom_text(aes(y = -1, label = N), size = 3)
  
  return(gg)
}

order = get_prop_latitude("Order")
order2gg = plot_prop_latitude(order, minimum = 30)
ggsave(order2gg, filename = "analysis/output/Botero14_proportion_Lat.range_order.pdf")
family = get_prop_latitude("Family.name")
family2gg = plot_prop_latitude(family, minimum = 30)
ggsave(family2gg, filename = "analysis/output/Botero14_proportion_Lat.range_family.pdf")
genus = get_prop_latitude("Genus")
genus2gg = plot_prop_latitude(genus, minimum = 30)
ggsave(genus2gg, filename = "analysis/output/Botero14_proportion_Lat.range_genus.pdf")






##### Finds all subclades from the phylogenies that meet certain criteria (#tips, proportion of "Mixed", etc.)
library(geiger)
jetz = read.tree("analysis/data/Jetz_etal/EricsonStage2_0001_1000/AllBirdsEricson1_reduced.tre")
# tree = phy = rcoal(10);plot(phy)
# tipLabels = sample(c("Temperate", "Tropical", "Mixed"), 10, TRUE)
get_matching_clades = function(phy, clade.size, tipLabels, minimum = TRUE, mixed = 1, multiLat = 0){
  # Finds all subclades from the phylogenies that meet certain criteria:
  #   - phy, phylogeny that'll serve as base
  #   - clade.size, minimum (or maximum) clade size
  #   - tipLabels, ordered (to match 'phy' tip.label) character vector indicating the latitudinal distribution of the species, one of: "Temperate", "Tropical", "Mixed"
  #   - minimum = TRUE, is 'clade.size' the minimum species number desired?
  #   - mixed = 1, proportion of "Mixed" species permited in the subclades
  #   - multiLat = 0, proportion of 'Tropical' species allowed in 'Temperate' subclades (and vice-versa)
  #
  # Returns one of the three: 
  # - NULL, if there is NO clade that meets 'clade.size' (and 'minimum') parameters
  # - NA, if there is NO clade that meets 'mixed' and 'multiLat' parameters
  # - a list of numeric vectors indicating the tips that meet the criteria. Each vector has 2 attributes: "species", with the corresponding tip.label; and "latitude", with the latitudinal distribution ("Temperate", "Tropical", "Mixed")
  
  if(mixed < 0 | mixed > 1){
    stop("'mixed' indicates proportion of 'Mixed' species in the subclades,\nso it must be between 0 and 1!")
  }
  if(multiLat < 0 | multiLat > 1){
    stop("'multiLat' indicates proportion of 'Tropical' species allowed in 'Temperate' subclades (and vice-versa),\nso it must be between 0 and 1!")
  }
  
  # get species names from phylogeny
  real.nms = phy$tip.label
  # assigns the latitudinal distribution to the tip.labels
  phy$tip.label = tipLabels
  
  # function to get all subclades from a phylogeny
  foo = ape:::prop.part(phy)
  # 'prop.part' returns numeric vectors with the tip numbers and tip.labels as an attribute
  nms = attr(foo, which = "labels")
  
  # filters and tests if any subclade meets 'clade.size' criteria
  if(minimum){
    foo2 = foo[sapply(foo, length) >= clade.size]
  } else{
    foo2 = foo[sapply(foo, length) <= clade.size]
  }
  if(length(foo2) == 0){
    return(NULL)
  }
  
  # gets the latitudinal distribution for all subclades of interest
  these.nms = lapply(foo2, function(x) nms[x])
  # gets the species names for all subclades of interest
  sp.nms = lapply(foo2, function(x) real.nms[x])
  # assigns "species" and "latitude" as attributes
  for(i in 1:length(foo2)){
    attr(foo2[[i]], which = "latitude") = these.nms[[i]]
    attr(foo2[[i]], which = "species") = sp.nms[[i]]
  }
  
  # creates a function to evaluate the subclades according to 'mixed' and 'multiLat' parameters
  if(multiLat == 0){
    fun = function(x){
      bool1 = ("Tropical" %in% x) != ("Temperate" %in% x)
      bool2 = (sum(x == "Mixed") / length(x)) <= mixed
      return(bool1 & bool2)
    }
  } else if(multiLat == 1){
    fun = function(x){
      bool = (sum(x == "Mixed") / length(x)) <= mixed
      return(bool)
    }
  } else{
    fun = function(x){
      p = c(sum(x == "Tropical"),
            sum(x == "Temperate"),
            sum(x == "Mixed")
      ) / length(x)
      bool1 = (p[1] <= multiLat) | (p[2] <= multiLat)
      bool2 = p[3] <= mixed
      return(bool1 & bool2)
    }
  }
  
  # filters and tests if any subclade follows the 'mixed' and 'multiLat' parameters
  bool = sapply(these.nms, FUN = fun)
  if(sum(bool) == 0){
    return(NA)
  }
  out = foo2[bool]
  
  return(out)
}

library(parallel)
Ncores = 3
# gets the tips for which we have data (appears in Botero14 dataset)
tips2keep = lapply(jetz, function(x) which(!x$tip.label %in% birds$SPECIES) )
# prunes phylogenies accordingly
clean_jetz = mcmapply(phy = jetz, tip = tips2keep, FUN = drop.tip,
                      SIMPLIFY = FALSE, mc.cores = Ncores)
#   Creates another classification of "latitudinal distrib" based on the centroid.
# Basically, eliminates the "Mixed" category and assigns it either to "Tropical" 
# or "Temperate" depending on where the centroid is located.
birds$LAT.CENTROID = as.character(ifelse(abs(birds$CENTROID) < 23, "Tropical", "Temperate"))
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

save(subclades, file = "analysis/data/Jetz_Botero14_subclades_30min_0.3mixed_0multilat.RData")



# # Let's check the percentage of "Mixed" species in the subclades for each phylogeny
# foo=lapply(these_subclades, function(yyy)sapply(yyy, function(x){x=attr(x, "latitude");(sum(x == "Mixed") / length(x))}))
# lat = lapply(these_subclades, sapply, function(x) ifelse("Temperate" %in% attr(x, "latitude"), "Temp", "Trop"))
# ggplot(propMixed, aes(x = phy, y = prop)) + geom_jitter(height = 0) + ylim(0,0.3) + facet_grid(latitude ~ .)

# Let's check the proportion of Temp in the Trop clades (and vice-versa)
foo = lapply(these_subclades, function(yyy) sapply(yyy, function(x){
  x = attr(x, "latitude")
  x = sum(x == "Temperate") / length(x)
  ifelse(x > 0.5, 1-x, x)
}))
cntrl = mapply(x = 1:10, times = sapply(these_subclades, length), FUN = rep)
lat = lapply(these_subclades, sapply, function(x) max(attr(x, which = "latitude")))
propMixed = data.frame(prop = unlist(foo),
                       phy = unlist(cntrl),
                       latitude = unlist(lat),
                       N = unlist(sapply(these_subclades, sapply, length)))

ggplot(propMixed, aes(x = phy, y = N, group = phy, color = latitude)) + geom_jitter(height = 0)
ggplot(propMixed, aes(x = N, y = prop, group = phy)) + geom_jitter(height = 0)+facet_grid(latitude ~ .)

sum(propMixed$latitude=="Temp")/nrow(propMixed)

# checks for possible AGE x SpeciesNumber correlation
get_age = function(x) max(x$edge.length)
get_harshness = function(x) mean(birds$ENV.HARSHNESS[match(x$tip.label, birds$SPECIES)])
testSubclades = data.frame(tips = unlist(sapply(subclades, sapply, Ntip)),
                           age = unlist(sapply(subclades, sapply, get_age)),
                           prop = unlist(foo),
                           phy = unlist(cntrl),
                           latitude = unlist(lat),
                           harshness = unlist(sapply(subclades, sapply, get_harshness)))
ggplot(testSubclades, aes(x = age, y = tips, color = latitude)) + geom_point() + geom_smooth(method='lm')
ggplot(testSubclades, aes(x = harshness, y = tips)) + geom_point() + geom_smooth(method='lm') + facet_grid(.~latitude, scales = "free")



get_Botero_spp = function(these, col.sp = "SPECIES", col.subsp = "SUBSPECIES", ref = birds){
  if(class(these) == "phylo"){
    these = these$tip.label
  }
  if(class(these) != "character"){
    stop("'these' must be either a character vector or a phylogeny with the species names.")
  }
  
  ind = match(these, ref[[col.sp]])
  if(any(is.na(ind))){
    warning("Some species names were not found in the 'ref' and will be excluded from analysis!")
    ind = na.exclude(ind)
  }
  
  out = ref[ind, col.subsp]
  
  return(out)
}
get_prior_exp = function(subspp, output = "function", log = TRUE){
  LLexp = function(r) {
    out = suppressWarnings(dexp(x = subspp, rate = r, log = TRUE))
    -sum(out)
  }
  par = stats4:::mle(LLexp, start = list(r = .1))
  if(output == "fit"){
    return(par)
  }
  if(output == "coef"){
    return(coef(par))
  }
  priorb1exp = function(x) dexp(x, rate = coef(par), log = log)
  return(priorb1exp)
}

spp = lapply(subclades, sapply, get_Botero_spp)
priorb1exp = lapply(spp, lapply, get_prior_exp, output = "function", log = TRUE)

save(priorb1exp, file = "analysis/data/Jetz_Botero14_subclades_priorb1exp.RData")

