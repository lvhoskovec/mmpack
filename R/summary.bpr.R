#' Summary method for class "bpr"
#'
#' @param object object of class "bpr"
#' @param ... ignored
#' @return list of summary data with components
#' \itemize{
#'      \item exposure.response: cluster means, SD, 0.95 CI, and size, plus mean/SD of repsonse for subjects in each cluster 
#'      \item risk: mean risk for each subject
#'      \item cluster.summary: mean, SD of each pollutant in each cluster
#'      \item cluster.exposures: array of estimated exposure means for each pollutant in each cluster
#'      \item risk.summary: mean, SD, and 0.95 CI of cluster-intercept (risk) for each subject
#'      \item clusters: vector of cluster assignments for each subject from the best clustering  
#'      \item risk.distn: distribution of model-averaged cluster-intercepts for each cluster in the best clustering
#'      \item groupList: list of which subjects belong to which cluster in the best clustering 
#'      \item rho: posterior mean probability of inclusion for each exposure 
#' }
#' 
#' @method summary bpr
#' 
#' @export
#'
#'
summary.bpr <- function(object, ...){

  #### model average health effects (i.e. exposure-response function) ####
  niter <- length(object$alpha)
  Y <- object$Y
  C <- object$C # maximum clusters
  p <- ncol(object$X) # number of exposures
  Z_best <- bestcluster(object)
  Z <- recode.Z(Z_best)

  # (risk.distn)
  theta.star <- matrix(modelaves(object)$theta.star, ncol = length(unique(Z)))
  
  if(object$scaleY == TRUE){
    theta.star <- theta.star*sd(object$Y) + mean(object$Y) 
  }
  
  count <- rep(0, length(unique(Z)))
  for(c in 1:length(unique(Z))){
    count[c] <- length(which(Z == c))
  }
  
  
  mean.response <- rep(0, length(unique(Z)))
  sd.response <- rep(0, length(unique(Z)))
  for(c in 1:length(unique(Z))){
    mean.response[c] <- mean(Y[which(Z == c)])
    sd.response[c] <- sd(Y[which(Z == c)])
  }
  
  sum.fun <- function(x){
    cbind(mean(x), sd(x), quantile(x, .025), quantile(x, .975))
  }
  
  # summary of exposure-response functions for each unique cluster
  er <- data.frame(t(apply(theta.star, 2, sum.fun)), count, mean.response, sd.response)
  rownames(er) <- paste("cluster", seq(length(unique(Z))))
  colnames(er) <- c("mean risk", "SD", "95% CI Lower", "95% CI Upper",
                    "n", "mean(Y)", "SD(Y)")
  
  # estimated mean risk (h.hat)
  risk <- er$`mean risk`[Z]
  
  ## summary table of risk
  risk.summary <- er[Z,(1:4)]
  rownames(risk.summary) <- NULL
  
  # mean and SD for exposures in each cluster
  clusters <- array(NA, dim = c(p, 2, length(unique(Z))), 
                    dimnames = list(paste("exposure", seq(1:p)), c("mean", "SD"),
                                    paste("cluster", seq(1:length(unique(Z))))))
  
  # distribution of estimated exposure means for each cluster
  c <- 1
  clust.expo.means <- array(NA, dim = c(niter, p, length(unique(Z_best))))
  for(i in unique(Z_best)){
    est.sd <- matrix(NA,niter,p)
    for(s in 1:niter){
      est.sd[s,] <- sqrt(diag(chol2inv(chol(object$SigInv[s,,,i]))))
    }
    clust.expo.means[,,c] <- object$mu[,i,]
    clusters[,,c] <- matrix(c(apply(object$mu[, i,  ], 2, mean), 
                              apply(est.sd,2,mean)), ncol = 2, nrow = p)
    c <- c+1
  }
  
  groupList <- list()
  for (c in 1:length(unique(Z))){
    groupList[[c]] <- which(Z == c)
  }
  
  rho <- apply(object$rho, 2, mean)
  
  
  list1 <- list(exposure.response = er, 
                risk = risk, cluster.summary = clusters,
                cluster.exposures = clust.expo.means,
                risk.summary = risk.summary, 
                clusters = Z, risk.distn = theta.star,
                groupList = groupList,
                rho = rho)
  
  class(list1) <- "summary.bpr"
  return(list1)
  
}
