#' Simulate Data for MixModelPack
#'
#' function to simulate exposure, covariate, and response data
#'
#' @param n sample size 
#' @param datX exposure data
#'
#' @return list with components
#' \itemize{
#'    \item Y: response data
#'    \item X: exposure data
#'    \item W: covariate data
#'    \item h: exposure-response function
#'    \item active: active main effects
#'    \item active.ints: active interactions 
#' }
#'
#' @export


simexpodat <- function(n, datX){

  W <- matrix(rnorm(n*10), n, 10)
  samps <- sample(1:1000, n, replace = FALSE)
  X <- datX[samps,]
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
  
  return(list(Y=Y, X = X, W = W, 
              h = h, active = active, active.ints = active.ints))
  
}

