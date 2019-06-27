#' Model Average Cluster-Specific Intercepts for Profile Regression
#'
#' computes model averaged cluster-specific intercepts (theta) from supervised or unsupervised Bayesian profile regression model 
#'
#' @param bpr object of class bpr 
#'
#' @importFrom stats sd quantile
#' @return summary of model averaged posterior estimates for cluster-specific parameters
#'


# adjusted for no intercept in model 6/28/18

modelaves <- function(bpr){
  
  Zbest <- bestcluster(bpr)
  Z_keep <- bpr$Z
  delta_keep <- bpr$delta
  C <- bpr$C
  
  theta_keep <- delta_keep[,1:C]
  
  niter <- nrow(theta_keep)

  thetaStar <- matrix(0, niter, C)

  for (c in 1:C){
    if(c %in% Zbest){
      whoBest <- which(Zbest == c)
      for (s in 1:niter){
        thetaStar[s,c] <- mean(theta_keep[s, Z_keep[s,whoBest]])
      }
    }
  }

  # whoBest is the indices of individuals who ended up in cluster c in the best clustering
  # we want, looking at only individuals whoBest, the theta values
  # for the clusters they were assigned to at each iteration

  # Z_keep[s,whoBest] is the cluster that the individuals who ended
  # up in the cluster c in the best clustering were assigned to in the
  # s iteration

  # thetaStar[s,c] is the mean of the theta values from the s^th iteration
  # of the people assigned to cluster c in the best clustering

  # then we grab the theta value assigned to these cluster(s) in the s iteration

  theta.star <- thetaStar[,c(which(colMeans(thetaStar)!=0))]
  # same order as thetaStar just re-numbered
  # numbering corresponds to reodered Z in recode.Z
  theta.center <- theta.star - apply(matrix(theta.star, nrow = niter), 1, mean)
  
  return(list(theta.star = theta.star, theta.center = theta.center))
}





