#' Predict method for class "npb"
#'
#' @param object object of class "npb"
#' @param ... ignored
#'
#' @return list with components
#' \itemize{
#'      \item fitted.vals: posterior mean fitted values for each subject
#'      \item fitted.distn: posterior distribution of fitted values for each subject  
#' }
#'
#' @export


predict.npb <- function(object, ...){
  
  # distribution of fitted values
  # add covariates, but not intercept as it is already in risk.distn
  fitted.distn <- summary(object)$risk.distn + object$gamma[,-1] %*% t(object$W[,-1]) 
  
  # posterior mean fitted values for each subject
  fitted.vals <- apply(fitted.distn,2,mean)

  # gamma0 + xbeta + wgamma
  return(list(fitted.vals = fitted.vals, fitted.distn = fitted.distn))
  
}
