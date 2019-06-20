#' Print method for class "npb"
#'
#' @param x object of class "npb"
#' @param ... ignored
#' @return list of output from summary.npb
#' 
#' @method print npb
#' 
#' @export
#'
#'
print.npb <- function(x,...){
  
  return(summary.npb(x))
  
}