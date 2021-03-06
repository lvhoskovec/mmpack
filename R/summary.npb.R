#' Summary method for class "npb"
#'
#' @param object object of class "npb"
#' @param ... ignored
#' @return list of summary data with components
#' \itemize{
#'     \item main.effects: mean, sd, and 0.95 CI for main effects
#'     \item interactions: mean, sd, and 0.95 CI for interactions
#'     \item selected.main: mean, sd, and 0.95 CI for main effects with PIPs > 0.5
#'     \item selected.int: mean, sd, and 0.95 CI for interactions with PIPs > 0.5
#'     \item risk.summary: mean, sd and 0.95 CI for exposure-response function "risk"
#'     \item risk: mean risk for each individual
#'     \item risk.distn: matrix, columns are the posterior distribution of risk for each individual 
#' } 
#'
#' @method summary npb
#' 
#' @export
#' 
 
summary.npb <- function(object,...){
  
  prob <- .5 # PIP threshold
  sum.fun <- function(x){
    cbind(mean(x), sd(x), quantile(x, .025), quantile(x, .975))
  }
  
  n <- length(object$Y)

  # main.effects
  main.effects <- cbind(t(apply(object$beta, 2, sum.fun)), apply(object$pip.beta, 2, mean)) 
  colnames(main.effects) <- c("Posterior Mean", "SD", "95% CI Lower", "95% CI Upper", "PIP")
  
  # interactions
  if(!is.null(object$zeta)){
    interactions <- cbind(t(apply(object$zeta, 2, sum.fun)), apply(object$pip.zeta, 2, mean)) 
    colnames(interactions) <- c("Posterior Mean", "SD", "95% CI Lower", "95% CI Upper", "PIP")
  }else interactions <- NULL
  
  # selected.main
  if(any(index <- which(colMeans(object$pip.beta) > prob))){
    index <- which(colMeans(object$pip.beta) > prob)
    num.maineffects <- length(index)
    output <- matrix(object$beta[,index], ncol = length(index))
    dat <- cbind(t(apply(output, 2, sum.fun)), apply(as.matrix(object$pip.beta[,index]), 2, mean))             
    rownames(dat) <- index
    colnames(dat) <- c("Posterior Mean", "SD", "95% CI Lower", "95% CI Upper","PIP")
  }else {
    dat <- NULL
    num.maineffects <- 0 
  }
  
  # selected.int
  if(!is.null(object$zeta)) {
    if(any(index.int <- which(colMeans(object$pip.zeta) > prob) )){
      index.int <- which(colMeans(object$pip.zeta) > prob) 
      num.int <- length(index.int)
      output.int <- matrix(object$zeta[,index.int], ncol = length(index.int)) 
      dat.int <- cbind(t(apply(output.int, 2, sum.fun)), apply(as.matrix(object$pip.zeta[,index.int]), 2, mean)) 
      rownames(dat.int) <- index.int
      colnames(dat.int) <- c("Posterior Mean", "SD", "95% CI Lower", "95% CI Upper", "PIP")
    }else {
      dat.int <- NULL
      num.int <- 0
    }
  }else {
    dat.int <- NULL
    num.int <- 0
  }
  
  # risk.distn for each subject
  if(!is.null(object$zeta)){
    post.h <- object$gamma[,1] + t(object$X %*% t(object$beta)) + t(object$Z %*% t(object$zeta))
  }else{
    post.h <- object$gamma[,1] + t(object$X %*% t(object$beta))
  }
  
  # risk.summary, mean risk for each subject
  posterior.h <- t(apply(post.h, 2, sum.fun))
  colnames(posterior.h) <- c("Posterior Mean", "SD", "95% CI Lower", "95% CI Upper")
  h.hat <- posterior.h[,1]

  # covariate posterior mean and credible intervals
  cov.sum <- t(apply(object$gamma, 2, sum.fun))
  colnames(cov.sum) <- c("Posterior Mean", "SD", "95% CI Lower", "95% CI Upper")
  
  
  list1 <- list(main.effects = main.effects, interactions = interactions,
                selected.main = dat, selected.int = dat.int,
                risk.summary = posterior.h, risk = h.hat,
                risk.distn = post.h, covariates = cov.sum)
  class(list1) <- "summary.npb"
  return(list1)

}



