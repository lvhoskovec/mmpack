#' Log of sum of exponentiated log-likelihood
#'
#' @param x vector of log-likelihoods
#'
#' @return return log of the sum of the exponentiated log-likelihoods



logsum <- function(x) {
  log(sum(exp(x - max(x)))) + max(x)
}
