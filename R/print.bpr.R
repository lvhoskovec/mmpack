#' Print method for class "bpr"
#'
#' @param x object of class bpr
#' @param ... ignored
#' @return call to summary.bpr
#'
#' @method print bpr
#' 
#' @export

print.bpr <- function(x,...){
  
  return(summary.bpr(x))
  
  
}