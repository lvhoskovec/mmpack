#' Predict method for class "bpr"
#'
#' @param object object of class "bpr"
#' @param ... ignored
#'
#' @return list with components
#' \itemize{
#'      \item fitted.vals: posterior mean fitted values for each subject
#'      \item fitted.distn: posterior distribution of fitted values for each subject  
#' }
#' @export
#'

predict.bpr <- function(object, ...){
  
  C <- object$C
  W <- object$W 
  q <- ncol(W)
  
  Z.best.bpr <- bestcluster(object)
  Z.bpr <- recode.Z(Z.best.bpr)
  theta.star <- modelaves(object)$theta.star 
  
  if(object$scaleY == TRUE){
    theta.star <- theta.star*sd(object$Y) + mean(object$Y) 
    Wgamma <- object$delta[,C+1:q] %*% t(object$W) * sd(object$Y)
  }else{
    Wgamma <- object$delta[,C+1:q] %*% t(object$W)
  }
  
  preds.distn <- matrix(theta.star, ncol = length(unique(Z.bpr)))[,Z.bpr] + Wgamma
  mean.preds <- apply(preds.distn,2,mean)
  
  return(list(fitted.vals = mean.preds, fitted.distn = preds.distn))
  
}