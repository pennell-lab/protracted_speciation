
library(reshape2)
library(ggplot2)
library(ape)
library(gridExtra)

plot.prior.post = function(coda, priors, par, ...){
  lcoda = length(coda)
  coda.median = sapply(coda,
                       FUN = function(x){
                         apply(x, MARGIN = 2, FUN = median)
                       })
  for(i in 1:4){
    plot(density(coda[[1]][, i]), main = colnames(coda[[1]])[i], ...)
    for(k in 2:lcoda){
      par(new = TRUE)
      plot(density(coda[[k]][, i]), lty = k,
           axes = FALSE, xlab = "", ylab = "",
           main = "", ...)
    }
    abline(v = coda.median[i, ], lty = 1:lcoda)
    abline(v = par[i], col = "green")
    curve(exp(priors[[i]](x)), add = TRUE, col = "red")
    legend("topright", legend = c("prior", paste0("chain",1:lcoda)),
           lty = c(1, 1:lcoda), col = c("red", rep("black", lcoda)),
           cex = 0.7)
  }
}




#############################################################
#   Function to inspect for possible correlations between
# parameters of a Bayesian model. Better than the function
# coda::pairs.
#
#   Developed by Florian Hartig
#############################################################
library(IDPmisc)
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="blue4", ...)
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, method = "spearman"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex * r)
}

betterPairs <- function(YourData){
  return(pairs(YourData, lower.panel=function(...) {par(new=TRUE);ipanel.smooth(...)}, diag.panel=panel.hist, upper.panel=panel.cor))
}