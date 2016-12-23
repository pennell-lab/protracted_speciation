context("opt_logLik")

library(geiger)
data(amphibia)
data(caniformia)
data(carnivores)
data(caudata)
data(chelonia)
data(geospiza)
data(whales)

phy = list(amphibia[[3]], caniformia[[1]], carnivores[[1]],
           caudata[[1]], chelonia[[1]], geospiza[[1]], whales[[1]])
bran = lapply(phy, branching.times)
bran[[8]] = 1:10
input = matrix(data = c(0.2,0.1,1,0.1,
                        0.18,0.1,1,0.15,
                        0.18,0.07,1.2,0.1), 
               ncol = 4, byrow = TRUE)
pars2 = data.frame(c(0, 1, 1, 2, 1, 1),
                   c(1, 0, 1, 1, 1, 1),
                   c(2, 2, 1, 1, 2, 2),
                   rep(0, 6),
                   c(rep("lsoda", 4), "lsodes", "vode"),
                   rep(0, 6),
                   rep(0, 6),
                   stringsAsFactors = FALSE)

testthat("opt_logLik is identical to pbd_logLik",{
  for(i in 1:nrow(pars2)){
    for(k in 1:length(branch)){
      # MY function
      my_function = opt_loglik(brts = branch[[k]], pars2 = pars2[i, ])
      my = apply(input, MARGIN = 1, FUN = my_function)
      
      # PBD's function
      pbd = apply(input, MARGIN = 1, FUN = pbd_loglik,
                  brts = branch[[k]],
                  pars2 = as.vector(pars2[i, ], mode = "character"))
      
      expect_that(my, is_identical_to(pbd))
    }
  }
})











