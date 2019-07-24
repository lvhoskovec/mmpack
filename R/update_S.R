#' Update allocation variables in NPB
#'
#' Function to update the allocations of regression coefficients to clusters in NPB model
#'
#' @param S allocation variable to update
#' @param p number of parameters
#' @param alpha.pi shape1 parameter for beta prior on pi0
#' @param beta.pi shape2 parameter for beta prior on pi0
#' @param beta regression coefficients corresponding to S
#' @param theta clusters corresponding to S
#' @param zeta other regression coefficients
#' @param Y response
#' @param W covariates
#' @param gamma covariate regression coefficients
#' @param X data for regression coefficients being updated
#' @param Z data for regression coefficients NOT being updated
#' @param sig2inv error precision
#' @param phi2inv precision for theta
#' @param mu mean for theta
#' @param alpha DP parameter for theta
#' @param iter MCMC iteration for debugging purposes
#'
#' @importFrom stats dnorm rnorm
#' @return vector of updated allocation variables that assign each regression coefficient (beta) to a specific cluster (theta)
#'

update_S <- function(S, p,
                     alpha.pi, beta.pi,
                     beta, theta, zeta,
                     Y, W, gamma,
                     X, Z,
                     sig2inv, phi2inv, mu,
                     alpha, iter){

  n <- length(Y)
  K <- length(theta)

  for(j in 1:p){

    # recalculate pcwoj from the updated S
    # Keep null cluster even if nothing is assigned to it
    pcwoj <- table(S)
    if(!any(beta==0)){
      pcwoj <- c(0, pcwoj)
      names(pcwoj)[1] <- "0"
    }

    # find the j^th cluster
    c <- S[j]
    # take out the j^th element
    pcwoj[c] <- pcwoj[c] - 1
    

    # relabel pcwoj, theta, K, and S if beta[j] was in a non-null cluster by itself
    if (pcwoj[c] == 0 & c != 1){
      pcwoj <- pcwoj[-c]
      theta <- theta[-c]
      K <- length(theta)
      for (i in 1:p){
        if (i == j) S[i] <- NA
        else S[i] <- which(theta == beta[i])
      }
    }

    loglik <- rep(NA,K+1)

    # calculate log-likelihoods #

    # set Ystar for this j
    if(!is.null(Z)){
    Ystar <- Y - (X[,-j] %*% beta[-j]) - (W %*% gamma) - (Z %*% zeta)
    }else Ystar <- Y - (X[,-j] %*% beta[-j]) - (W %*% gamma)


    # 1) Pr(S_j = 0): null group
    loglik[1] <- log(pcwoj[1]+alpha.pi) - log(p+alpha.pi + beta.pi - 1) +
      sum(dnorm(Ystar, 0, sqrt(1/sig2inv), log = T))


    if(K > 1){
      # 2) Pr(S_j = c): existing group
      for(c in 2:K){
        loglik[c] <- (log(p - pcwoj[1] + beta.pi - 1) - log(p + alpha.pi + beta.pi -1)  +
                        log(pcwoj[c]) - log(p + alpha - pcwoj[1] - 1)  +
                        sum(dnorm(Ystar, X[,j] * theta[c],
                                  sqrt(1/sig2inv), log = T)))
      }

    }

    # 3) Pr(S_j = c*): new group #
    loglik[K+1] <- (log(p - pcwoj[1] + beta.pi - 1) -
                      log(p + alpha.pi + beta.pi - 1) + log(alpha)
                    - log(p + alpha - pcwoj[1] - 1)
                    + (
                      (-n/2)*log(2*pi) + (n/2)*log(sig2inv) + .5*log(phi2inv) -
                        .5*log(sig2inv * sum(X[,j]^2) + phi2inv) -
                        .5*( sig2inv*sum(Ystar^2) + phi2inv*mu^2 -
                               ( ((sig2inv*sum(Ystar*X[,j]) + phi2inv*mu)^2)/
                                   (sig2inv*sum(X[,j]^2) + phi2inv) ) )
                    ) )


    # calculate likelihood
    lik <- exp(loglik - logsum(loglik))
    
    # update the allocation for beta_j
    S[j] <- sample(K+1, 1, prob = lik)

    # handle case 3)
    if(S[j] == K+1){
      # calculate mean and variance
      v <- 1/(sig2inv * sum( X[,j]^2 ) + phi2inv)
      m <- v * (sig2inv * sum( Ystar * X[,j] ) + phi2inv * mu)
      # draw a new theta value
      thetanew <- rnorm(1, m, sqrt(v))
      # add this new theta value to the theta vector
      theta <- c(theta, thetanew)
    }

    # update number of thetas
    K <- length(theta)
    # update beta given new updated S and theta
    beta <- theta[S]

  }

  return(list(S = S, beta = beta, theta = theta))

}
