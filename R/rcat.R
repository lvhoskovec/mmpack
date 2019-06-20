#' Random categorical generator
#'
#' @param probs n by K matrix of probabilities n subjects to be assigned to K categories 
#'
#' @return return vector of length nrow(probs) of categorical indicators
#'
#'
#'


rcat <- function(probs){
  return( apply ( t(apply(probs,1,cumsum)) >
                    runif(nrow(probs),0,rowSums(probs)),1,
                  function(x) min(which(x)) ) )
}
