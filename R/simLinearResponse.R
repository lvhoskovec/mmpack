#' Function to simulate response data as a linear function of predictors
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


simLinearResponse <- function(X, W){
  # need at least 2 predictors 
  if(ncol(X) < 4)stop("need at least 4 predictors")
  
  gamma <- rnorm(ncol(W), 0, 1)
  e.vec <- sample(1:ncol(X), 4, replace = FALSE)
  a <- e.vec[1]
  b <- e.vec[2]
  c <- e.vec[3]
  d <- e.vec[4]
  h <-  1*X[,a] - 1*X[,b] + 1*X[,c] - 1*X[,d] + 
    .7*X[,a]*X[,b] - .5*X[,c]*X[,d]
  Y <- h + W%*%gamma + rnorm(nrow(X),0,1) 
  active <- sort(e.vec)
  ints <- combn(1:ncol(X), 2)
  active.ints1 <- which(colSums(apply(ints, 2, function(x) {x %in% c(a,b)})) == 2)
  active.ints2 <- which(colSums(apply(ints, 2, function(x) {x %in% c(c,d)})) == 2)
  active.ints <- sort(c(active.ints1, active.ints2))
  
  return(list(Y=Y, h = h, active = active, active.ints = active.ints))
}