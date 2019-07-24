#' Wrapper function to fit profile regression model from package PReMiuM
#'
#' @param niter number of iterations through sampler
#' @param nburn number of burn-in iterations
#' @param Y response data
#' @param X exposure data
#' @param W covariate data
#' @param scaleY logical; if TRUE response is centered and scaled before model fit
#' @param varSelectType type of variable selection to be used: "None" or "BinaryCluster", default is "None" 
#' @param simnum number of output file (.txt files will be saved in the file path paste0("Premium_output/output",simnum)
#' @param priors prior hyperparameters
#' @param seed random number seed
#'
#'
#' @return list with components
#' \itemize{
#'     \item fit: model fit, an object of type "runInfoObj". See profRegr documentation in R package PReMiuM
#'     \item risk: mean risk for each subject
#'     \item rho: posterior distribution of probability of inclusion for each exposure 
#'     \item risk.summary: mean, SD, and 0.95 CI of cluster-intercept (risk) for each subject
#'     \item exposure.response: cluster means, SD, 0.95 CI, and size, plus mean/SD of repsonse for subjects in each cluster  
#'     \item fitted vals: fitted values for each subject
#'     \item clusters: vector of most optimal clustering 
#'     \item cluster.summary: array of empirical mean and SD of exposures for individuals assigned to each cluster
#'     \item risk distn: distribution of model-averaged cluster-intercepts for each cluster in the best clustering
#'     \item groupList: list of which subjects are in which group in the most optimal clustering 
#'     \item riskProfileObj: object of type "riskProfileObj". See PReMiuM documentation for function "calcAvgRiskAndProfile"
#'   
#'}
#' 
#' @import PReMiuM
#' @importFrom utils read.table
#' @export
#'

premium_wrapper <- function(niter, nburn, Y, X, W, scaleY = FALSE, varSelectType = "None", simnum = NULL,
                            priors, seed = NULL){
  
  if(nburn >= niter) stop("Number of iterations (niter) must be greater than number of burn-in iteractions (nburn)")
  
  p <- ncol(X)
  
  if(scaleY == TRUE){
    Y.save <- Y
    Y <- scale(Y)
  }else{
    Y.save <- Y
  }

  # create data frame for PReMiuM package
  input.dat <- data.frame(cbind(Y, X, W))
  
  colnames(input.dat) <- c("outcome", sprintf("var%d", seq(1,ncol(X))), 
                           sprintf("fixed%d", seq(1, ncol(W))))
  
  inputs.test <- list(inputData = input.dat, covNames = sprintf("var%d", seq(1,ncol(X))),
                      xModel = "Normal", yModel = "Normal", nCovariates = ncol(X), 
                      fixedEffectNames = sprintf("fixed%d", seq(1, ncol(W))))
  
  preds <- inputs.test$inputData[,-1]

  # hyper CANNNOT be NULL in profRegr, set some default values if priors is NULL
  if(is.null(priors)) priors <- list(aRho = 0.5, bRho = 0.5)
  
  # change: let nClusInit be default 7/7/18
  fit.prem <- profRegr(yModel=inputs.test$yModel,
                         xModel=inputs.test$xModel, sampler = "Truncated", 
                         nSweeps=niter-nburn, 
                         nBurn=nburn, data=inputs.test$inputData, 
                         output=paste0("Premium_output/output",simnum),
                         covNames = inputs.test$covNames,
                         fixedEffectsNames = inputs.test$fixedEffectNames,
                         hyper = priors, seed = seed, predict = preds,
                         varSelectType = varSelectType,
                         useHyperpriorR1 = FALSE)
  
  # get h.hat from premium
  # postprocessing, check these post-processing steps and parameter options
  dissimObj <- calcDissimilarityMatrix(fit.prem, onlyLS = TRUE)
  clusObj <- calcOptimalClustering(dissimObj, maxNClusters = 40, useLS = TRUE)
  clusters <- clusObj$clustering
  riskProfileObj <- calcAvgRiskAndProfile(clusObj)

  beta <- read.table(paste0("Premium_output/output",simnum,"_beta.txt"), sep = " ", header = FALSE)
  Wbeta <- t(W %*% t(beta)) # distribution of Wbeta for each individual
  
  if(scaleY == TRUE){
    Wbeta <- Wbeta*sd(Y.save)
  }
  
  Wbeta.hat <- apply(Wbeta, 2, mean)
  
  # (risk.distn) estimated risk at each sweep for each cluster 
  theta.star <- matrix(riskProfileObj$risk, ncol = length(unique(clusters)))
  
  if(scaleY == TRUE){
    theta.star <- theta.star*sd(Y.save) + mean(Y.save)
  }
  
  theta <- apply(theta.star, 2, mean)
  
  fitted.vals <- theta[clusters] + Wbeta.hat # predictions 

  count <- rep(0, length(unique(clusters)))
  for(c in 1:length(unique(clusters))){
    count[c] <- length(which(clusters == c))
  }
  
  mean.response <- rep(0, length(unique(clusters)))
  sd.response <- rep(0, length(unique(clusters)))
  for(c in 1:length(unique(clusters))){
    mean.response[c] <- mean(Y.save[which(clusters == c)])
    sd.response[c] <- sd(Y.save[which(clusters == c)])
  }
  
  # (exposure.response)
  er <- data.frame(apply(theta.star, 2, mean), 
                   apply(theta.star, 2, sd), 
                   apply(theta.star, 2, function(x) quantile(x, .025)), 
                   apply(theta.star, 2, function(x) quantile(x, .975)),
                   count, mean.response, sd.response)

  rownames(er) <- paste("cluster", seq(length(unique(clusters))))
  colnames(er) <- c("mean risk", "SD", "95% CI Lower", "95% CI Upper", 
                                   "n", "mean(Y)", "SD(Y)")

  # (risk)
  risk <- er$`mean risk`[clusters] # average risk for each individual
  
  # (risk.summary)
  risk.summary <- er[clusters,(1:4)]
  rownames(risk.summary) <- NULL
  
  if (varSelectType == "BinaryCluster"){
    rho <- summariseVarSelectRho(fit.prem)$rho
  }else rho <- NULL
  
  # empirical exposure summary for each cluster 
  # mean, sd of exposures for individuals assigned to each cluster 
  cluster.summary <- array(NA, dim = c(p, 2, length(unique(clusters))), 
                    dimnames = list(paste("exposure", seq(1:p)), c("mean", "SD"),
                                    paste("cluster", seq(1:length(unique(clusters))))))
  for(c in unique(clusters)){
    
    Xdata <- matrix(X[which(clusters==c),],nrow = length(which(clusters==c)),ncol = p)
    
    cluster.summary[,1,c] <- apply(Xdata, 2, mean)
    cluster.summary[,2,c] <- apply(Xdata, 2, sd)
  }
  
  
  groupList <- list()
  for (c in 1:length(unique(clusters))){
    groupList[[c]] <- which(clusters == c)
  }
  
  
  return(list(fit = fit.prem, risk = risk, 
              rho = rho, risk.summary = risk.summary,
              exposure.response = er,
              fitted.vals = fitted.vals, clusters = clusters,
              cluster.summary = cluster.summary,
              risk.distn = theta.star,
              groupList = groupList,
              riskProfileObj = riskProfileObj))
  
  
}
