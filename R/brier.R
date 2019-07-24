#' Calculate brier scores for Bayesian methods 
#'
#' @param X exposure data
#' @param active true active exposures 
#' @param pips posterior inclusion probabilities from model fit
#'
#' @return brier score

brier <- function(X, active, pips){
  b <- rep(0,ncol(X))
  b[active] <- 1
  return(mean((b - pips)^2))
}