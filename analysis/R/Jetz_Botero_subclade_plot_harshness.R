

library(parallel)
library(geiger)
library(coda)
library(dplyr)
library(ggplot2)

# setwd("~/Desktop/gitHub/protracted_sp/")
# source("analysis/R/auxiliary_functions.R") 
file_path = "analysis/data/subclade_harshness/" # path to the files

# name of the files
files = list.files(path = file_path, pattern = "*PBD_Bayes.txt")
# creates a data frame with the info about the runs - eg, category and chain.
cntrl = do.call(rbind, sapply(files, strsplit, split = "_"))
cntrl = data.frame(type = cntrl[ , 6],
                   phy = cntrl[ , 5],
                   category = cntrl[ , 7],
                   est = cntrl[ , 9],
                   chain = cntrl[ , 11],
                   row.names = NULL)

# Ncores = 3
# # read all the files, already skips the burn-in to accelarate the process
# raw = mclapply(paste0(file_path, files), FUN = read.csv, sep = "\t", skip = ceiling(5000000 * 0.3), mc.cores = Ncores)
# # assigns the correct colnames
# raw = lapply(raw, function(x){colnames(x) = c("logLik", "prior", "posterior", "b", "mu1", "la1", "mu2", "accepted", "proposed.par", "converged"); x})
# # discards the in between values; ie, keeps only the real samples from posterior
# coda.reads = lapply(lapply(raw, function(x)x[ , 1:7]), mcmc, thin = 1000)
# # creates clean data frames, with only the real samples and the 4 parameters!
# ind = seq(1, 3500001, by = 1000)
# reads = do.call(rbind, lapply(raw, function(x)x[ind, 1:7]))
# nr = length(ind)
# reads$type = rep(cntrl$type, each = nr)
# reads$phy = rep(cntrl$phy, each = nr)
# reads$category = rep(cntrl$category, each = nr)
# reads$est = rep(cntrl$est, each = nr)
# reads$chain = rep(cntrl$chain, each = nr)
# 
# save(coda.reads, file = "analysis/data/subclade_harshness/subclade_harshness.RData")
# write.csv(file = "analysis/data/subclade_harshness/subclade_harshness_reads.csv", reads, row.names = FALSE)

load("analysis/data/subclade_harshness/subclade_harshness.RData")
reads = read.csv(file = "analysis/data/subclade_harshness/subclade_harshness_reads.csv", header = TRUE)

# subspecies formation rate
reads %>%
  filter(type == "harshness") %>%
  ggplot(aes(b, fill = category)) + geom_histogram(bins = 100)
# time to speciation
reads %>%
  filter(type == "harshness") %>%
  ggplot(aes(la1, fill = category)) + geom_histogram(bins = 100)
# species extinction rate
reads %>%
  filter(type == "harshness") %>%
  ggplot(aes(mu1, fill = category)) + geom_histogram(bins = 100)
# subspecies extinction rate
reads %>%
  filter(type == "harshness") %>%
  ggplot(aes(mu2, fill = category)) + geom_histogram(bins = 100)

# species diversification rate
reads %>%
  filter(type == "harshness") %>%
  ggplot(aes(b - mu1, fill = category)) + geom_histogram(bins = 100)
# subspecies diversification rate
reads %>%
  filter(type == "harshness") %>%
  ggplot(aes((1/la1) - mu2, fill = category)) + geom_histogram(bins = 100)





# subspecies formation rate
reads %>%
  filter(type == "harshnessClean") %>%
  ggplot(aes(b, fill = category)) + geom_histogram(bins = 100)
# time to speciation
reads %>%
  filter(type == "harshnessClean") %>%
  ggplot(aes(la1, fill = category)) + geom_histogram(bins = 100)
# species extinction rate
reads %>%
  filter(type == "harshnessClean") %>%
  ggplot(aes(mu1, fill = category)) + geom_histogram(bins = 100)
# subspecies extinction rate
reads %>%
  filter(type == "harshnessClean") %>%
  ggplot(aes(mu2, fill = category)) + geom_histogram(bins = 100)

# species diversification rate
reads %>%
  filter(type == "harshnessClean") %>%
  ggplot(aes(b - mu1, fill = category)) + geom_histogram(bins = 100)
# subspecies diversification rate
reads %>%
  filter(type == "harshnessClean") %>%
  ggplot(aes((1/la1) - mu2, fill = category)) + geom_histogram(bins = 100)
