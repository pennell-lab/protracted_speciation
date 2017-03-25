

# https://susanejohnston.wordpress.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/
# Posted on August 9, 2012 by susanejohnston
ggplotRegression <- function (fit) {
  require(ggplot2)
  if(length(class(fit)) == 1){ # it's a regular "lm"
    R2 = summary(fit)$adj.r.squared
  } else{
    R2 = cor(fit$y,predict(fit))^2
  }
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(R2, 5),
                       "Intercept =",signif(coef(fit)[1],5 ),
                       " Slope =",signif(coef(fit)[2], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5))) 
}