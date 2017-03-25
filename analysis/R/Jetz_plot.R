



library(ape)

load(file = "analysis/data/Botero14.RData")

jetz = read.tree("analysis/data/Jetz_etal/EricsonStage2_0001_1000/AllBirdsEricson1_reduced.tre")


# get the list of species that are NOT in the list
absent = which(!jetz[[1]]$tip.label %in% birds$SPECIES)
# only for the first phylogeny
jj = drop.tip(jetz[[1]], tip = absent)
# match the regions to the tips
sp = lapply(b1names[1:3], function(xx) birds$SPECIES[birds$LAT.RANGE %in% xx])
edg = lapply(sp, function(x) sapply(x, function(x) which(jj$tip.label == x)))
tip.col = rep("black", length(jj$edge))
tip.col[which(jj$edge[ , 2] %in% edg[[2]])] = "red" # Tropical
tip.col[which(jj$edge[ , 2] %in% edg[[3]])] = "blue" # Temperate
# plot the phylogeny
png("analysis/output/Jetz_regions.png", width = 750, height = 750)
par(mai = rep(0, 4))
plot(jj, type = "fan", show.tip.label = FALSE, edge.color = tip.col)
legend("topright", legend = c("Tropical", "Temperate"), lty = 1, col = c("red", "blue"), bty = "n")
dev.off()


# by region
miss = lapply(sp, function(x) which(!jetz[[1]]$tip.label %in% x))
jj = lapply(miss, drop.tip, phy = jetz[[1]])