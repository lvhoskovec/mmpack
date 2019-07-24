#' Summarize evaluation of estimated exposure-response function for linear models 
#'
#' @param fit object of type lm
#' @param h true exposure-response function
#' @param interaction logical; was LM fit with interactions?
#' @param X exposure data
#' @param W covariate data
#'
#' @import stats
#' @return summary of estimated exposure-response function


sumERlinearmodels <- function(fit, h, interaction, X, W){
  
  ints <- combn(1:ncol(X), 2)
  Z <- apply(ints, 2, FUN = function(x) {
    X[,x[1]] * X[,x[2]]
  })
  p <- ncol(X) # main
  d <- ncol(Z) # ints
  q <- ncol(W) # covariates
  
  lm.sum <- summary(fit)
  
  if(interaction == TRUE){
    X.all <- cbind(1, X, Z) # exposure main effects and interactions
    beta.main <- lm.sum$coefficients[,1][1:(p+1)] # include intercept
    beta.int <- lm.sum$coefficients[,1][(p+q+2):(p+q+d+1)]
    beta.all <- c(beta.main, beta.int)
    h.hat <- X.all%*%beta.all
    y.hat <- cbind(1, X, W, Z)%*%lm.sum$coefficients[,1]
    
    # covariance matrix for exposure response 
    Xmat <- model.matrix(fit)
    # var.delta.hat <- lm.sum$sigma^2 * solve(t(Xmat) %*% Xmat) # diagonal elements are se for coefficients
    var.delta.hat <- vcov(fit)
    var.beta.hat.main <- var.delta.hat[(1:(p+1)),(1:(p+1))]
    var.beta.hat.int <- var.delta.hat[(p+q+2):(p+q+d+1), (p+q+2):(p+q+d+1)]
    cov.beta.main.int <- var.delta.hat[(1:(p+1)), (p+q+2):(p+q+d+1)]
    var.allbeta.hat <- cbind(rbind(var.beta.hat.main, t(cov.beta.main.int)), 
                             rbind(cov.beta.main.int, var.beta.hat.int))
    
    # sqrt(diag(var.allbeta.hat)) = se for regression coefficients for exposures
    var.h.hat <- X.all %*% var.allbeta.hat %*% t(X.all)
    
    
    
  }else{
    X.all <- cbind(1,X)
    beta.all <- lm.sum$coefficients[,1][1:(p+1)]
    h.hat <- X.all%*%beta.all
    y.hat <- cbind(1, X, W)%*%lm.sum$coefficients[,1]
    
    # covariance matrix for exposure response 
    Xmat <- model.matrix(fit)
    var.delta.hat <- lm.sum$sigma^2 * solve(t(Xmat) %*% Xmat) # all coefficients 
    var.beta.hat <- var.delta.hat[(1:(p+1)),(1:(p+1))] # only exposures 
    var.h.hat <- X.all %*% var.beta.hat %*% t(X.all)
    
  }
  
  hfit <- lm(h.hat ~ h)
  # MSE
  MSE <- sqrt(mean((h - h.hat)^2))
  
  # confidence interval for exposure-response function is h.hat +/- t_crit*se(h.hat) for each i
  sd.h <- sqrt(diag(var.h.hat)) 
  # critical value
  # crit <- qt(.975, lm.sum$df[2]) # not correct - multiple testing 
  # Working-Hotelling for 95% confidence curves 
  # df = (p, n-p) for estimating sigma^2
  crit <- sqrt(2*qf(.95, lm.sum$df[1], lm.sum$df[2]))
  
  lwr.h <- h.hat - crit*sd.h
  upr.h <- h.hat + crit*sd.h
  postsum <- cbind(h, h.hat, sd.h, lwr.h, upr.h)
  cvg <- mean(as.numeric((postsum[,1] > postsum[,4]) & (postsum[,1] < postsum[,5])))
  h.sum <- c(hfit$coefficients[1], hfit$coefficients[2], summary(hfit)$r.squared, MSE, cvg)
  
  return(h.sum)
  
}
