#' Function to simulate response data as a piece-wise constant function of predictor data (fixed profiles)
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
#'        \item silhouette: silhouette statistic for clustering observations
#'        \item gapWidth: gap width statistic for clustering 
#' }
#' @export

simProfilesResponse <- function(X, W){
  
  e.vec <- sample(1:ncol(X), 2, replace = FALSE)
  a <- e.vec[1]
  b <- e.vec[2]
  gamma <- rnorm(ncol(W),0,1)
  n <- nrow(X)
  low <- which(X[,a] <= median(X[,a]) & X[,b] <= median(X[,b]))
  med <- which(X[,a] <= median(X[,a]) & X[,b] > median(X[,b]))
  med.hi <- which(X[,a] > median(X[,a]) & X[,b] <= median(X[,b])) 
  hi <- which(X[,a] > median(X[,a]) & X[,b] > median(X[,b]))
  h <- rep(NA, n)
  h[low] <- -2
  h[med] <- -1
  h[med.hi] <- 0
  h[hi] <- 2
  Y <- rep(NA, n)
  Y[low] <- h[low] + W[low,] %*% gamma + rnorm(length(low),0,1)
  Y[med] <- h[med] + W[med,] %*% gamma + rnorm(length(med),0,1)
  Y[med.hi] <- h[med.hi] + W[med.hi,] %*% gamma + rnorm(length(med.hi),0,1)
  Y[hi] <-  h[hi] + W[hi,] %*% gamma + rnorm(length(hi),0,1)
  active <- sort(e.vec)
  ints <- combn(1:ncol(X), 2)
  active.ints <- which(colSums(apply(ints, 2, function(x) {x %in% c(a,b)})) == 2)
  
  # silhouette
  clust <- recode.Z(h)
  d <- dist(X)
  sil <- silhouette(clust, dist = d)
  save.sil <- sil[,3]
  
  # gap width 
  FUNcluster <- function(X, k=4){
    clust = list()
    clust$cluster = recode.Z(h)
    return(clust)
  }
  FUNcluster(X)
  gaptab <- clusGap(x=X, FUNcluster = FUNcluster, K.max = 10, B = 500) 
  gap <- gaptab$Tab[,3]
  
  
  return(list(Y=Y, h=h, active = active, active.ints = active.ints,
              silhouette = save.sil, gapWidth = gap))
  
}