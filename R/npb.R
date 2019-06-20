#' Function to fit nonparametric Bayes shrinkage
#'
#' Fits nonparametric Bayes shrinkage model for continuous response data using Markov chain Monte Carlo (MCMC) methods
#'
#' @param niter number of iterations
#' @param nburn number of burn-in iterations
#' @param X n by p matrix of predictor data
#' @param Y n by 1 vector of continuous response data
#' @param W n by q matrix of covariate data 
#' @param scaleY logical; if TRUE response will be centered and scaled before model fit, default is FALSE
#' @param priors list of prior hyperparameters, see package documentation for details 
#' @param interact logical; if TRUE (default) include all pairwise interactions of predictors, if FALSE include main effects only
#'
#' @return an object of class "npb", which has the associated methods:
#' \itemize{
#' \item \code{\link{print}} (i.e. \code{\link{print.npb}} )
#' \item \code{\link{summary}} (i.e. \code{\link{summary.npb}} )
#' \item \code{\link{predict}} (i.e. \code{\link{predict.npb}} )
#' }
#' 
#' 
#' @return a list with components
#' \itemize{
#'    \item X: predictor data matrix
#'    \item Y: response data vector
#'    \item Z: matrix of pairwise multiplicative interactions of predictor data (only if interact = TRUE)
#'    \item W: covariate data
#'    \item alpha: DP parameter for main effects estimates 
#'    \item alpha.2: DP parameter for interactions estimates 
#'    \item beta: unscaled main effect regression coefficient estimates 
#'    \item zeta: unscaled interaction regression coefficient estimates 
#'    \item mu: mu estimates 
#'    \item mu.2: mu2 estimates 
#'    \item phi2inv: phi2inv estimates
#'    \item phi2inv.2: phi2inv.2 estimates 
#'    \item gamma: unscaled covariate regression coefficient estimates 
#'    \item sig2inv: sig2inv estimates 
#'    \item pip.beta: inclusion indicators for main effects
#'    \item pip.zeta: inclusion indicators for interactions
#' }
#' @export
#'

npb <- function(niter, nburn, X, Y, W, scaleY = FALSE, priors, interact = FALSE){
  
  if(interact == FALSE){
    fit <- npb_main(niter = niter, nburn = nburn, X = X, Y = Y, W = W, 
                    scaleY = scaleY, priors = priors)
  }else{
    fit <- npb_int(niter = niter, nburn = nburn, X = X, Y = Y, W = W, 
                   scaleY = scaleY, priors = priors)
  }
  
  return(fit)

}
