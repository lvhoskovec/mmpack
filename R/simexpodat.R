#' Simulate Data for MixModelPack
#'
#' Function to simulate exposure, covariate, and response data for up to 1000 observations 
#'
#' @param n sample size 
#' @param Xdat Xdat exposure data matrix loaded from package mmpack
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


simexpodat <- function(n, Xdat){
  
  if(n > 1000) stop("sample size cannot be larger than 1000")

  W <- matrix(rnorm(n*10), n, 10)
  samps <- sample(1:1000, n, replace = FALSE)
  X <- Xdat[samps,]
  X <- apply(X,2,scale)
  gamma <- rnorm(ncol(W), 0, 1)
  e.vec <- c(3,4,5,7)
  a <- e.vec[1]
  b <- e.vec[2]
  c <- e.vec[3]
  d <- e.vec[4]
  h <-  3*X[,a] - 2*X[,b] + 2.5*X[,c] - 4*X[,d] + 
    .3*X[,a]*X[,b] - .6*X[,c]*X[,d]
  Y <- h + W%*%gamma + rnorm(nrow(X),0,1) 
  
  active <- sort(e.vec)
  ints <- combn(1:ncol(X), 2)
  active.ints1 <- which(colSums(apply(ints, 2, function(x) {x %in% c(a,b)})) == 2)
  active.ints2 <- which(colSums(apply(ints, 2, function(x) {x %in% c(c,d)})) == 2)
  active.ints <- sort(c(active.ints1, active.ints2))
  
  return(list(Y = Y, X = X, W = W, 
              h = h, active = active, active.ints = active.ints))
  
}

