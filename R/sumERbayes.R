#' Summarize evaluation of estimated exposure-response function for Bayesian methods 
#'
#' @param fit fit or summary from one of the 5 Bayesian methods 
#' @param h true exposure-response function
#'
#' @return summary of estimated exposure-response function


sumERbayes <- function(fit,h){
  # fit is a summary typically
  hfit <- lm(fit$risk ~ h) # est. exposure response regressed on true 
  MSE <- sqrt(mean((h - fit$risk)^2))
  
  # h.hat 
  postsum <- (cbind(h, fit$risk.summary))
  # postsum has: true h, post mean, post sd, lwr .025, upper .975
  cvg <- mean(as.numeric((postsum[,1] > postsum[,4]) & (postsum[,1] < postsum[,5])))
  
  h.sum <- c(hfit$coefficients[1], hfit$coefficients[2], summary(hfit)$r.squared, MSE, cvg)
  return(h.sum)
}
