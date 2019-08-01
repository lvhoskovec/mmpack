#' Best Clustering by Least Squared Distance to Probability Matrix
#'
#' @param bpr object of class bpr; output from function profileReg
#'
#' @return vector denoting the best clustering of the data using least squares distance to a similarity matrix as described in Dahl (2006)
#' @references Dahl, D. B. (2006). Model-Based Clustering for Expression Data via a Dirichlet Process Mixture Model. Bayesian Inference for Gene Expression and Proteomics, 201–218. https://doi.org/10.1017/CBO9780511584589.011
#' 
#'


bestcluster <- function(bpr){
  
  Z_keep <- bpr$Z
  C <- bpr$C

  n <- ncol(Z_keep)
  niter <- nrow(Z_keep)

  # calculate probability matrix P (symmetric)
  Prob <- matrix(NA, n, n)
  for (i in 1:n){
    for (j in 1:i){
      Prob[i, j] <- Prob[j, i] <- mean(Z_keep[,i]==Z_keep[,j])
    }
  }

  # create similarity matrix S for each clustering in Z_keep
  # calculate squared distance from S to P
  LSdist <- rep(NA, niter)
  for (s in 1:niter){
    S <- matrix(0, n, n)
    for (c in 1:C){
      S[which(Z_keep[s,] == c), which(Z_keep[s,] == c)] <- 1
    }
    LSdist[s] <- sum( (S - Prob)^2 )

  }

  # select Z_keep that corresponds to the S that minimizes LS distance to P
  best <- which.min(LSdist)
  return(Z_keep[best,])


}
