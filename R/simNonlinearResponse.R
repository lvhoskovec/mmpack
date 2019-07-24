#' Function to simulate response data as a nonlinear function of predictors
#'
#' @param X predictor data
#' @param W covariate data
#'
#' @return list with components
#' \itemize{ 
#'        \item Y: response data
#'        \item h: evaluation of exposure-response function for each individual 
#'        \item active: main effects selected to be in the exposure-response function
#'        \item active.ints: interactions selected to be in the exposure-response function
#' }
#' @export
#' 
#' 
simNonlinearResponse <- function(X, W){
  
  if(ncol(X)<3)stop("need at least 3 predictors")
  
  gamma <- rnorm(ncol(W),0,1)
  e.vec <- sample(1:ncol(X), 3, replace = FALSE)
  a <- e.vec[1]
  b <- e.vec[2]
  c <- e.vec[3]
  h <- 2/(1 + exp(-3*X[,a])) - 2/(1 + exp(-5*X[,b])) + 
    2/(1 + exp(-5*X[,c])) - .4*X[,a]*X[,b]
  Y <- h + W %*% gamma + rnorm(nrow(X),0,1)
  active <- sort(e.vec)
  ints <- combn(1:ncol(X), 2)
  active.ints <- which(colSums(apply(ints, 2, function(x) {x %in% c(a,b)})) == 2)
  
  return(list(Y = Y, h = h, active = active, active.ints = active.ints))
  
}