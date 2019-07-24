#' Summary method for class "bpr"
#'
#' @param object object of class "bpr"
#' @param ... ignored
#' @return list of summary data with components
#' \itemize{
#'      \item exposure.response: cluster means, SD, 0.95 CI, and size, plus mean/SD of repsonse for subjects in each cluster 
#'      \item risk: mean risk for each subject
#'      \item cluster.summary: array of empirical mean and SD of exposures for individuals assigned to each cluster
#'      \item risk.summary: mean, SD, and 0.95 CI of cluster-intercept (risk) for each subject
#'      \item clusters: vector of most optimal clustering  
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
  niter <- length(object$alpha) # number of iterations after burn-in
  Y <- object$Y # response
  X <- object$X
  n <- length(Y)
  C <- object$C # maximum clusters
  p <- ncol(X) # number of exposures
  Z_best <- bestcluster(object) # best clustering of subjects
  Z <- recode.Z(Z_best) # recoding best clustering to sequential order
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
  
  
  # what might be more interesting is the empirical distribution of exposures 
  # mean,sd of exposures for the individuals assigned to each cluster 
  
  clusters <- array(NA, dim = c(p, 2, length(unique(Z))), 
                    dimnames = list(paste("exposure", seq(1:p)), c("mean", "SD"),
                                    paste("cluster", seq(1:length(unique(Z))))))
  for(c in unique(Z)){
    
    Xdata <- matrix(X[which(Z==c),],nrow = length(which(Z==c)),ncol = p)
    
    clusters[,1,c] <- apply(Xdata, 2, mean)
    clusters[,2,c] <- apply(Xdata, 2, sd)
  }
  
  groupList <- list()
  for (c in 1:length(unique(Z))){
    groupList[[c]] <- which(Z == c)
  }
  
  rho <- apply(object$rho, 2, mean)
  
  
  list1 <- list(exposure.response = er, 
                risk = risk, cluster.summary = clusters,
                risk.summary = risk.summary, 
                clusters = Z, risk.distn = theta.star,
                groupList = groupList,
                rho = rho)
  
  class(list1) <- "summary.bpr"
  return(list1)
  
}
