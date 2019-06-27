#' Recode Z for BPR
#'
#' @param Z.best output from bestcluster function
#'
#' @return recoded Z
#'

recode.Z <- function(Z.best){
  
  # retains the order of Z
  # appropriate to use on theta.star in bestcluster
  
  levels <- sort(unique(Z.best))
  for (i in 1:length(levels)){
    Z.best[which(Z.best == levels[i])] <- i
  }
  return(Z.best)
}
