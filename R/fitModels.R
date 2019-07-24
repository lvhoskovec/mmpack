#' Fit all 5 Bayesian methods and summarize output 
#'
#' @param names optional name for the simulation (i.e. "linear")
#' @param simnum simulation number
#' @param niter number of total iterations for each method
#' @param nburn number of burn-in iterations for each method
#' @param X exposure data
#' @param W covariate data
#' @param data list returned from function simLinearResponse, simNonlinearResponse, or simProfilesResponse
#' @param seed seed
#'
#' @return data frame summarizing each method
#' @export
#'
fitModels <- function(names=NULL, simnum, niter, nburn, X, W,
                       data, seed = NULL){
  scaleY <- TRUE
  p <- ncol(X)
  q <- ncol(W)
  Y <- data$Y
  h <- data$h
  active <- data$active
  active.ints <- data$active.ints
  priors.npb <- NULL
  priors.prem <- NULL
  priors.bpr <- NULL
  set.seed(seed) 
  fit.npbr <- npb(niter = niter, nburn = nburn, X = X, Y = Y, W = W, 
                  scaleY = scaleY, priors = priors.npb, interact = FALSE)
  set.seed(seed) 
  fit.npb <- npb(niter = niter, nburn = nburn, X = X, Y = Y, W = W, 
                 scaleY = scaleY, priors = priors.npb, interact = TRUE)
  set.seed(seed)  
  fit.bkmr <- bkmr_wrapper(niter = niter, nburn = nburn, Y = Y, X = X, W = W, 
                           varsel = TRUE) 
  set.seed(seed)
  fit.premium <- premium_wrapper(niter = niter, nburn = nburn,
                                 Y = Y, X = X, W = W, scaleY = scaleY, varSelectType = "BinaryCluster",
                                 simnum = simnum, priors = priors.prem, seed = seed)
  set.seed(seed) 
  fit.upr <- profileReg(niter = niter, nburn = nburn, X = X, Y = Y, W = W, scaleY = scaleY, 
                        varsel = TRUE, priors = priors.bpr, sup = FALSE)
  
  # normal linear models
  fit.lm.main <- lm(Y ~ X + W)
  fit.lm.int <- lm(Y ~ .^2 + W, data = data.frame(X))
  
  # get exposure response function estimate
  NPBr.h <- sumERbayes(fit=summary(fit.npbr), h=h)
  NPB.h <- sumERbayes(fit=summary(fit.npb), h=h)
  bkmr.h <- sumERbayes(fit=fit.bkmr, h=h)
  prem.h <- sumERbayes(fit=fit.premium, h=h)
  upr.h <- sumERbayes(fit=summary(fit.upr), h=h)
  lm.main.h <- sumERlinearmodels(fit=fit.lm.main, h=h, interaction = FALSE, X = X, W = W)
  lm.int.h <- sumERlinearmodels(fit=fit.lm.int, h=h, interaction = TRUE, X = X, W = W)
  
  NPBr.pip <- c(mean(colMeans(fit.npbr$pip.beta)[active] > 0.5),
                mean(colMeans(fit.npbr$pip.beta)[-active] > 0.5),
                brier(X, active, colMeans(fit.npbr$pip.beta)))
  NPB.pip <- c(mean(colMeans(fit.npb$pip.beta)[active] > 0.5),
               mean(colMeans(fit.npb$pip.beta)[-active] > 0.5),
               brier(X, active, colMeans(fit.npb$pip.beta)))
  BKMR.pip <- c(mean(fit.bkmr$pips[active] > 0.5), 
                mean(fit.bkmr$pips[-active] > 0.5),
                brier(X, active, fit.bkmr$pips))
  upr.pip <- c(mean(colMeans(fit.upr$rho)[active] > 0.5), 
               mean(colMeans(fit.upr$rho)[-active] > 0.5), 
               brier(X, active, colMeans(fit.upr$rho)))
  prem.pip <- c(mean(colMeans(fit.premium$rho)[active] > 0.5),
                mean(colMeans(fit.premium$rho)[-active] > 0.5),
                brier(X, active, colMeans(fit.premium$rho)))
  
  is.integer0 <- function(x)
  {
    is.integer(x) && length(x) == 0L
  }
   # LM
  lm.main.actors <- which(summary(fit.lm.main)$coefficients[-1,4] <= .05)
  lm.main.actors <- lm.main.actors[which(lm.main.actors <= p)]
  if(is.integer0(lm.main.actors)){
    tsr.lm.main <- 0
    fsr.lm.main <- 0
  }else{
    tsr.lm.main <- mean(active %in% lm.main.actors)
    fsr.lm.main <- mean(seq(1:p)[-active] %in% lm.main.actors)
  }
  lm.main.pip <- c(tsr.lm.main, fsr.lm.main, NA) 
  
  # LM-int
  lm.int.actors <- which(summary(fit.lm.int)$coefficients[-1,4] <= .05)
  lm.int.actors <- lm.int.actors[which(lm.int.actors <= p)]
  if(is.integer0(lm.int.actors)){
    tsr.lm.int <- 0
    fsr.lm.int <- 0
  }else{
    tsr.lm.int <- mean(active %in% lm.int.actors)
    fsr.lm.int <- mean(seq(1:p)[-active] %in% lm.int.actors)
  }
  lm.int.pip <- c(tsr.lm.int, fsr.lm.int, NA) 
  
  pips <- rbind(NPBr.pip, NPB.pip, upr.pip,
                prem.pip, BKMR.pip, lm.main.pip, lm.int.pip)
  
  colnames(pips) <- c("trueSelectRate", "falseSelectRate", "brierScore")
  
  ints <- combn(1:ncol(X), 2)
  Z <- apply(ints, 2, FUN = function(x) {
    X[,x[1]] * X[,x[2]]
  })
  
  if(!is.null(active.ints)){
    NPB.pip.ints <- c(mean(colMeans(fit.npb$pip.zeta)[active.ints] > 0.5),
                      mean(colMeans(fit.npb$pip.zeta)[-active.ints] > 0.5),
                      brier(Z, active.ints, colMeans(fit.npb$pip.zeta)))
    
    pip.ints <- matrix(NA, 7, 3)
    pip.ints[2,] <- NPB.pip.ints
    colnames(pip.ints) <- c("tsr", "fsr", "brier")
    lm.int.actints <- which(summary(fit.lm.int)$coefficients[-1,4] <= .05)
    lm.int.actints <- lm.int.actints[which(lm.int.actints > p + q)] - (p+q)
    
    if(is.integer0(lm.int.actints)){
      tsr.lm.int.ints <- 0
      fsr.lm.int.ints <- 0
    }else{
      tsr.lm.int.ints <- mean(active.ints %in% lm.int.actints)
      fsr.lm.int.ints <- mean(seq(1:ncol(Z))[-active.ints] %in% lm.int.actints)
    }
    
    lm.int.ints.pip <- c(tsr.lm.int.ints,  fsr.lm.int.ints, NA) # not real pips
    pip.ints[7,] <- lm.int.ints.pip
    
    
  }else {
    pip.ints <- matrix(NA, 7, 3)
    colnames(pip.ints) <- c("tsr", "fsr", "brier")
  }
    
  
  df <- data.frame(method = c(paste0("AME.",names), paste0("EMM.",names), paste0("UDR.",names),
                              paste0("PReMiuM.",names), paste0("BKMR.", names),
                              paste0("main.lm.",names), paste0("int.lm.", names)),
                   intercept = c(NPBr.h[1], NPB.h[1], upr.h[1], prem.h[1], bkmr.h[1], lm.main.h[1], lm.int.h[1]),
                   slope = c(NPBr.h[2], NPB.h[2], upr.h[2], prem.h[2], bkmr.h[2], lm.main.h[2], lm.int.h[2]),
                   R2 = c(NPBr.h[3], NPB.h[3], upr.h[3], prem.h[3], bkmr.h[3], lm.main.h[3], lm.int.h[3]),
                   MSE = c(NPBr.h[4], NPB.h[4], upr.h[4], prem.h[4], bkmr.h[4], lm.main.h[4], lm.int.h[4]),
                   Cvg = c(NPBr.h[5], NPB.h[5], upr.h[5], prem.h[5], bkmr.h[5], lm.main.h[5], lm.int.h[5]))
  
  df.vs <- cbind(df, pips, pip.ints)
  rownames(df.vs) <- NULL
  
  return(df.vs)
}


