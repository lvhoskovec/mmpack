#' Ordinary Least Squares
#'
#' @param dat simulated data for npb model
#'
#' @return least squares estimates of beta 
#' 
#' @importFrom stats lm
#' @export
#'
#'

ols <- function(dat){

  # data
  X <- dat$X # pollutant coefficients
  W <- dat$W # covariate coefficients
  # health response
  Y1 <- dat$Y1 # sigma = sd(xbeta)
  Y2 <- dat$Y2 # sigma = 1
  Y3 <- dat$Y3 # sigma = 100

  # least squares estimates
  beta.ols1 <- lm(Y1 ~ X)
  beta.ols2 <- lm(Y2 ~ X)
  beta.ols3 <- lm(Y3 ~ X)


  ols.sum1 <- summary(beta.ols1)
  betahat.ols1 <- as.numeric(ols.sum1$coefficients[,1])
  ols.sum2 <- summary(beta.ols2)
  betahat.ols2 <- as.numeric(ols.sum2$coefficients[,1])
  ols.sum3 <- summary(beta.ols3)
  betahat.ols3 <- as.numeric(ols.sum3$coefficients[,1])


  return(list(betahat.ols1 = betahat.ols1,
              betahat.ols2 = betahat.ols2,
              betahat.ols3 = betahat.ols3))

}
