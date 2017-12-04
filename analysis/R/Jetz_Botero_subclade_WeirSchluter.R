
# library(gdata) # for excel files -> "read.xls"
library(ape)
library(ggplot2)
library(gridExtra)
library(tidyverse)

#1) latitudinal ranges and approximates dates of divergence for sister species of birds
# wsBirds = read.xls("~/Dropbox/PhD/Projects/protracted_sp/Weir_Schluter_2007_SM.xls", 
#                sheet = 2, header = TRUE, skip = 1, as.is = TRUE)[ , 1:21]
wsBirds = read_csv("~/Dropbox/PhD/Projects/protracted_sp/Weir_Schluter_1 latitudinal ranges and dates of divergence for sister sp of birds.csv")
head(wsBirds)

#2) latitudinal ranges and approximates dates of divergence for sister species of mammals
# wsMammals = read.xls("~/Dropbox/PhD/Projects/protracted_sp/Weir_Schluter_2007_SM.xls", 
#                sheet = 3, header = TRUE, skip = 1, as.is = TRUE)[ , 1:20]
wsMammals = read_csv("~/Dropbox/PhD/Projects/protracted_sp/Weir_Schluter_2 latitudinal ranges and dates of divergence for sister sp of mammals.csv")
head(wsMammals)

#3) latitudinal ranges and maximum haplotype divergence within bird species
# wsBirdsHaplo = read.xls("~/Dropbox/PhD/Projects/protracted_sp/Weir_Schluter_2007_SM.xls", 
#                    sheet = 4, header = TRUE, as.is = TRUE)[ , 1:13]
wsBirdsHaplo = read_csv("~/Dropbox/PhD/Projects/protracted_sp/Weir_Schluter_3 latitudinal ranges and maximum haplotype divergence within bird sp.csv")
head(wsBirdsHaplo)

#4) latitudinal ranges and maximum haplotype divergence within mammal species
# wsMammalsHaplo = read.xls("~/Dropbox/PhD/Projects/protracted_sp/Weir_Schluter_2007_SM.xls", 
#                         sheet = 5, header = TRUE, as.is = TRUE)[ , 1:13]
wsMammalsHaplo = read_csv("~/Dropbox/PhD/Projects/protracted_sp/Weir_Schluter_4 latitudinal ranges and maximum haplotype divergence within mammal sp.csv")
head(wsMammalsHaplo)


#Phylogroup heading indicates wether haplotype variation within a species is divided into geographically segregated phylogroups.
lm_eqn <- function(m){
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq))
}
get_jetz_names = function(phy, sp1, sp2){
  # finds the corresponding phylogeny tree node to extract ages
  nodeBirds = matrix(c(match(sp1, phy$tip.label), match(sp2, phy$tip.label)), 
                     ncol = 2, byrow = FALSE)
  # function to get ages by checking the 'edge' and 'edge.length' objects
  fun = function(x){
    if(all(is.na(x))){ # no sp of the sp-pair was found
      return( structure(c(NA, NA), names = c("age", "sister.sp")) )
    }
    
    if(any(is.na(x))){ # only one sp was found
      x = na.exclude(x)
      age = phy$edge.length[match(x, phy$edge[ , 2])]
      out = structure(c(age, NA), names = c("age", "sister.sp"))
    } else{ # both sp were found
      age = mean(phy$edge.length[match(x, phy$edge[ , 2])])
      sister = diff(phy$edge[match(x, phy$edge[ , 2]), 1]) == 0
      out = structure(c(age, sister), names = c("age", "sister.sp"))
    }
    return(out)
  }
  out = t(apply(nodeBirds, 1, fun))
}
get_botero_ssp = function(x){
  y = which(botero$SPECIES == x)
  if(length(y) > 0){
    return(mean(botero$SUBSPECIES[y], na.rm = TRUE))
  } else{
    return(NA)
  }
}

# replaces any parethesis in sp names
spBirds1 = gsub(x = wsBirds$"Species 1", pattern = " *\\(.*?\\) *", replacement = "")
spBirds2 = gsub(x = wsBirds$"Species 2", pattern = " *\\(.*?\\) *", replacement = "")
# formats the sp names according to Jetz
spBirds1 = gsub(x = spBirds1, pattern = " ", replacement = "_")
spBirds2 = gsub(x = spBirds2, pattern = " ", replacement = "_")
# loads one of Jetz trees
jetz = read.tree("~/Desktop/gitHub/protracted_sp/analysis/data/Jetz_etal/EricsonStage2_0001_1000/AllBirdsEricson1_reduced.tre")
# gets ages of sister species in jetz trees
agesFULL = lapply(jetz, get_jetz_names, sp1 = spBirds1, sp2 = spBirds2)
# gets the average ages for the sister species
ages = do.call(rbind, agesFULL) %>% 
  as.data.frame() %>% 
  mutate(id = rep(1:(length(spBirds1)), times = length(jetz))) %>%
  group_by(id) %>% 
  summarise_all(c("mean", "sd"))
# load Botero '14 dataset on subspecies
botero = read_csv("~/Desktop/gitHub/protracted_sp/analysis/data/Botero14_plus_taxonomy.csv")
# matches the sp in WeirSchluter dataset and Botero and extracts ssp numbers
ssp = apply(matrix(c(spBirds1, spBirds2), ncol = 2, byrow = FALSE), 1, get_botero_ssp)
# make the dataset
timeBirds = tibble(sp1 = spBirds1,
                   sp2 = spBirds2,
                   tree = ages$age_mean,
                   tree.sd = ages$age_sd,
                   sister = ages$sister.sp_mean,
                   ws = wsBirds$Time,
                   ssp = ssp)
lmBirds = lm(tree ~ ws, data = timeBirds)
corBirds = cor(timeBirds$ws, timeBirds$tree, use = "na.or.complete")
ggFull = ggplot(timeBirds, aes(ws, tree)) + 
  geom_point() + 
  #geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  geom_text(x = 1.5, y = 17, label = paste0("cor = ", round(corBirds, 3))) +
  #geom_text(x = 1.5, y = 19, label = lm_eqn(lmBirds), parse = TRUE) + 
  labs(x = "Divergence time (Weir Schluter)", y = "Divergence time (avg. Jetz)")
ggFull




# replaces any parethesis in sp names
spBirdsHaplo = gsub(x = wsBirdsHaplo$Species, pattern = " *\\(.*?\\) *", replacement = "")
# formats the sp names according to Jetz
spBirdsHaplo = paste(wsBirdsHaplo$Genus, spBirdsHaplo, sep = "_")
# gets ages of sister species in jetz trees
agesFULLHaplo = lapply(jetz, get_jetz_names,
                       sp1 = spBirdsHaplo, sp2 = rep(NA, length(spBirdsHaplo)) )
# gets the average ages for the sister species
agesHaplo = lapply(agesFULLHaplo, function(x)x[ , 1]) %>%
  do.call(cbind, .) %>% 
  apply(., 1, mean)
# matches the sp in WeirSchluter dataset and Botero and extracts ssp numbers
sspHaplo = sapply(spBirdsHaplo, get_botero_ssp)
# make the dataset
timeBirdsHaplo = tibble(sp = spBirdsHaplo,
                        tree = agesHaplo,
                        haplotype = wsBirdsHaplo$Time,
                        ssp = sspHaplo)
lmBirdsHaplo = lm(tree ~ haplotype, data = timeBirdsHaplo)
corBirdsHaplo = cor(timeBirdsHaplo$haplotype, timeBirdsHaplo$tree, use = "na.or.complete")
ggFullHaplo = ggplot(timeBirdsHaplo, aes(haplotype, tree)) + 
  geom_point() + 
  #geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  geom_text(x = 1.5, y = 12, label = paste0("cor = ", round(corBirdsHaplo, 3))) +
  #geom_text(x = 1.5, y = 13, label = lm_eqn(lmBirdsHaplo), parse = TRUE) + 
  labs(x = "Divergence time (Haplotypes)", y = "Species age (avg. Jetz)")
ggFullHaplo
ggsave(plot = grid.arrange(ggFull, ggFullHaplo),
       filename = "analysis/output/Jetz_WeirSchluter.pdf")




lmBotero = lm(ssp ~ ws, data = timeBirds)
corBotero = cor(timeBirds$ws, timeBirds$ssp, use = "na.or.complete")
ggFullBotero = ggplot(timeBirds, aes(ws, ssp)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  geom_text(x = 6, y = 19, label = paste0("cor = ", round(corBotero, 3))) +
  geom_text(x = 6, y = 17, label = lm_eqn(lmBotero), parse = TRUE) + 
  labs(x = "Divergence time (Weir Schluter)", y = "Number of subspecies (Botero)")
ggFullBotero

lmBoteroHaplo = lm(ssp ~ haplotype, data = timeBirdsHaplo)
corBoteroHaplo = cor(timeBirdsHaplo$haplotype, timeBirdsHaplo$ssp, use = "na.or.complete")
ggHaploBotero = ggplot(timeBirdsHaplo, aes(haplotype, ssp)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  geom_text(x = 3.5, y = 35, label = paste0("cor = ", round(corBoteroHaplo, 3))) +
  geom_text(x = 3.5, y = 31, label = lm_eqn(lmBoteroHaplo), parse = TRUE) + 
  labs(x = "Divergence time (Haplotypes)", y = "Number of subspecies (Botero)") 
ggHaploBotero
ggsave(plot = grid.arrange(ggFullBotero, ggHaploBotero),
       filename = "analysis/output/Jetz_WeirSchluter_ssp.pdf")



ggSisterHist = ggplot(timeBirds, aes(sister, y = ..count../sum(..count..))) +
  geom_histogram(bins = 10) +
  labs(x = "Sister species recovery in Jetz", y = "Proportion")
ggSisterHist
ggSister = timeBirds %>%
  mutate(dif = ws - tree) %>%
  ggplot(., aes(sister, dif)) + 
  geom_jitter(width = 0.01, height = 0) +
  geom_smooth(method = "lm") +
  labs(x = "Proportion of sister species recovery in Jetz", y = "Diference in age (ws-Jetz)")
ggSister
ggSisterSD = ggplot(timeBirds, aes(sister, tree.sd)) + 
  scale_y_log10() +
  geom_jitter(width = 0.01, height = 0) +
  geom_smooth(method = "lm") +
  labs(x = "Proportion of sister species recovery in Jetz", y = "LOG age sd (Jetz)")
ggSisterSD

ggsave(plot = grid.arrange(ggSisterHist, ggSister, ggSisterSD),
       filename = "analysis/output/Jetz_WeirSchluter_sisterSp_comparison.pdf")

