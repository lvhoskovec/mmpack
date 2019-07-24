#' Wrapper function to fit bkmr model in package bkmr.
#'
#' @param niter number of iterations of sampler
#' @param nburn number of burn-in iterations 
#' @param Y n x 1 vector of response data
#' @param X n x p matrix of exposure data
#' @param W n x q matrix of covariate data
#' @param varsel logical; if TRUE variable selection is implemented, default is FALSE
#' @param groups optional vector of length p indicating group membership to conduct hierarchical variable selection 
#'  
#'
#' @return a list with components 
#' \itemize{
#'     \item fit: the original model fit 
#'     \item risk: the posterior mean risk for each subject
#'     \item pips: component-wise (if groups = NULL) or conditional (if groups is specified) posterior inclusion probabilities for exposures
#'     \item group.pips: group posterior inclusion probabilities (only if groups is specified))
#'     \item risk.summary: mean, SD, and 0.95 CI for risk for each subject 
#'     \item preds: predicted values for each subject based on exposures and covariates
#'     \item risk.distn: posterior distribution of risk for each subject
#' }
#' @import bkmr
#' @export


bkmr_wrapper <- function(niter, nburn, Y, X, W, varsel = FALSE, groups = NULL){
  
  if(nburn >= niter) stop("Number of iterations (niter) must be greater than number of burn-in iteractions (nburn)")
  
  
  fit.bayes <- kmbayes(y = Y, Z = X, X = W, iter = niter, 
                       varsel = varsel, groups = groups, est.h=TRUE)
  

  # the following generates posterior samples from h(x)
  h.star <- SamplePred(fit = fit.bayes, Znew = NULL, Xnew = NULL, 
                        Z = X, X = NULL, sel = seq(nburn+1,niter))
  
  posterior.h <- cbind(apply(h.star,2,mean), 
                       apply(h.star,2,sd),
                       apply(h.star,2,FUN = function(x) quantile(x, .025)),
                       apply(h.star,2,FUN = function(x) quantile(x, .975)))
  
  risk <- posterior.h[,1] # posterior mean, estimated exposure response function h
  
  if(is.null(groups)){
    group.pips <- NULL
    pips <- CalcPIPs(fit.bayes, sel = seq(nburn+1, niter))
  }else{
    group.pips <- CalcGroupPIPs(fit.bayes, sel = seq(nburn+1, niter))
    pips <- CalcWithinGroupPIPs(fit.bayes, sel = seq(nburn+1, niter))
  }

  samps <- SamplePred(fit.bayes, Znew = X, Xnew = W) # original exposures
  preds <- apply(samps, 2, mean) # posterior predicted values
  
  return(list(fit = fit.bayes, risk = risk, pips = pips, 
              group.pips = group.pips,
              risk.summary = posterior.h, preds = preds,
              risk.distn = h.star))
  
}
