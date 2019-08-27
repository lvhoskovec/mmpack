#' Fit nonparametric Bayes shrinkage model with main effects only
#'
#' @param niter number of iterations
#' @param nburn number of burn-in iterations
#' @param X n by p matrix of predictor data
#' @param Y n by 1 vector of response data
#' @param W n by q matrix of covariate data 
#' @param scaleY logical; if TRUE response will be centered and scaled before model fit
#' @param priors list of prior hyperparameters
#' @param intercept logical; indicates if an overall intercept should be estimated with covariates
#'
#' @importFrom stats var rnorm rbeta runif rgamma 
#' 
#' 
#' @return list of model estimates, see npb documentation
#'

npb_main <- function(niter, nburn, X, Y, W, scaleY = FALSE, priors, intercept = TRUE){
  
  if(nburn >= niter) stop("Number of iterations (niter) must be greater than number of burn-in iteractions (nburn)")
  

  ##########
  # Priors #
  ##########

  if(missing(priors)) priors <- NULL
  if(is.null(priors$a.sig)) priors$a.sig <- 1 # shape param for gamma prior on sig2inv
  if(is.null(priors$b.sig)) priors$b.sig <- 1 # rate param for gamma prior on sig2inv
  if(is.null(priors$a.phi1)) priors$a.phi1 <- 1 # shape param for gamma prior on phi2inv, precision of base distribution for main effects
  if(is.null(priors$b.phi1)) priors$b.phi1 <- 1 # rate param for gamma prior on phi2inv, precision of base distribution for interactions
  if(is.null(priors$alpha.a)) priors$alpha.a <- 2 # shape param for gamma prior on alpha, DP parameter for main effects
  if(is.null(priors$alpha.b)) priors$alpha.b <- 1 # rate param for gamma prior on alpha, DP parameter for main effects
  # intercept
  if(is.null(priors$mu.0)) priors$mu.0 <- 0 # mean parameter for normal prior on gamma_0, intercept
  if(is.null(priors$kappa2inv.0)) priors$kappa2inv.0 <- 1 # precision parameter for normal prior on gamma_0, intercept
  # covariates
  if(is.null(priors$mu.gamma) & !is.null(W)) priors$mu.gamma <- rep(0, ncol(W)) # mean vector for normal prior on covariates
  if(is.null(priors$kap2inv)  & !is.null(W)) priors$kap2inv <- rep(1, ncol(W)) # precision vector for normal prior on covariates 
  # # #
  if(is.null(priors$sig2inv.mu1)) priors$sig2inv.mu1 <- 1 # precision param for normal prior on mean of base distribution for main effects
  if(is.null(priors$alpha.pi)) priors$alpha.pi <- 1 # shape1 parameter for beta prior on pi0 for main effects
  if(is.null(priors$beta.pi)) priors$beta.pi <- 1 # shape2 parameter for beta prior on pi0 for main effects
  
  # if intercept = TRUE add a column of 1's to the covariates
  if(intercept == TRUE){
    if(is.null(W)) {
      W <- matrix(1, length(Y), 1)
    }else {W <- cbind(1, W)}
    priors$mu.gamma <- c(priors$mu.0, priors$mu.gamma)
    priors$kap2inv <- c(priors$kappa2inv.0, priors$kap2inv)
  }
  # if intercept = FALSE then we use W as is 
  
  if(scaleY == TRUE){
    Y.save <- Y
    Y <- scale(Y)
  }else{
    Y.save <- Y
  }
  
  #######################
  ### Starting values ###
  #######################

  ################
  # main effects #
  ################
  
  X <- as.matrix(X)
  p <- ncol(X) # number of pollutants
  n <- length(Y) # sample size
  beta <- rnorm(p) # regression coefficients
  theta <- unique(c(0, beta)) # theta[1] = 0 for null group
  K <- length(theta) # unique clusters

  S <- rep(NA, p) # allocation variable
  for (j in 1:p){
    S[j] <- which(theta == beta[j])
  }

  phi2inv <- 1 # precision for base distribution G1 
  mu <- mean(X %*% beta) # mean for base distribution G1
  alpha <- 1 # DP parameter for main effects

  ##############
  # covariates #
  ##############

  if(!is.null(W)){
    W <- as.matrix(W)
    q <- ncol(W) # number of covariates, including overall intercept
    gamma <- rnorm(q) # regression coefficients 
  }else{
    q <- 0
    gamma <- NULL
  }


  #########
  # error #
  #########

  sig2inv <- 1 # error precision

  #####################
  ### Storage space ###
  #####################

  K_keep <- rep(NA, niter) # number of unique clusters
  alpha_keep <- matrix(NA, niter) # DP parameter
  beta_keep <- matrix(NA, niter, p) # main effects
  mu_keep <- rep(NA, niter) # mean of base distribution
  phi2inv_keep <- rep(NA, niter) # precision of base distribution
  if(!is.null(W)) gamma_keep <- matrix(NA, niter, q) # covariates
  sig2inv_keep <- rep(NA, niter) # error precision
  pip.beta <- matrix(NA, niter, p) # beta PIPs

  ############
  ### MCMC ###
  ############

  for (s in 1:niter){

    ################
    # update alpha #
    ################

    eta <- rbeta(1, alpha + 1, p)
    pi.eta <- (priors$alpha.a + K - 1)/(p*(priors$alpha.b - log(eta)))/
      ((priors$alpha.a + K - 1)/(p*(priors$alpha.b - log(eta))) + 1)
    if (runif(1) < pi.eta){
      alpha <- rgamma(1, priors$alpha.a + K, priors$alpha.b - log(eta))
    }else{
      alpha <- rgamma(1, priors$alpha.a + K - 1, priors$alpha.b - log(eta))
    }

    ##################################################
    # update S: allocation variable for main effects #
    ##################################################

    up.S <- update_S(S = S, p = p, alpha.pi = priors$alpha.pi, beta.pi = priors$beta.pi,
                     beta = beta, theta = theta, zeta = NULL,
                     Y = Y, W = W, gamma = gamma,
                     X = X, Z = NULL, sig2inv = sig2inv,
                     phi2inv = phi2inv, mu = mu, alpha = alpha, iter = s)


    S <- up.S$S
    beta <- up.S$beta
    theta <- up.S$theta
    K <- length(theta)

    ################################
    # block update theta and gamma #
    ################################

    if(length(theta)>1){
      
      # reparameterize main effects to remove nulls
      # create the matrix T1 indicating to which theta each beta belongs 
      # first column for theta = 0
      # rows for beta, columns for theta
      # a single 1 in each row
      
      T1 <- matrix(0, p, K)
      for(c in 1:K){
        T1[which(S==c),c] <- 1
      }
      
      # clear out the empty clusters/columns in T1 and the null cluster
      T1 <- as.matrix(T1[which(beta!=0),-1]) 
      # clear out the columns in X with null betas
      X0 <- as.matrix(X[, which(beta != 0)]) 
      Ttheta <- T1 %*% theta[-1]
      Xbeta <- X0 %*% Ttheta # equivalent to full matrix X %*% beta
      XT <- X0 %*% T1 # n by length(theta)-1 design matrix for updating theta

      # start precision vector
      precision.vector <- rep(phi2inv, K-1)

    }else {
      # if length(theta) = 1, all betas are 0
      Xbeta <- 0
      XT <- NULL
      precision.vector <- NULL
    }

    # add covariate precisions to precision.vector
    precision.vector <- c(precision.vector, priors$kap2inv)

    # update delta = c(theta, gamma)
    A <- cbind(XT, W)
    if(!is.null(A)){
      m.star <- c(rep(mu, K-1), priors$mu.gamma)    
      v <- chol2inv(chol(sig2inv * crossprod(A) + 
                           diag(x=precision.vector, nrow = length(precision.vector))))
      m <- v %*% (sig2inv * (t(A)%*%Y) + precision.vector * m.star)
      delta <- drop(m + t(chol(v)) %*% rnorm(K+q-1))
      gamma <- delta[K:(K+q-1)]
      Adelta <- A %*% delta
    }else{
      Adelta <- 0 # no covariates, no intercept, and all betas are 0 
    }

    if(K > 1){
      theta <- delta[1:(K-1)]
      theta <- c(0, theta)
    }


    # update beta deterministically
    beta <- theta[S]

    ##################
    # update sig2inv #
    ##################

    sig2inv <- rgamma(1, priors$a.sig + n/2, priors$b.sig + .5*sum( (Y - Adelta)^2)  )

    ##################
    # update phi2inv #
    ##################

    if(K > 1){
      phi2inv <- rgamma(1, priors$a.phi1 + (K-1)/2,
                        priors$b.phi1 + .5*sum( (theta[-1] - mu)^2 ) )
    }else{
      phi2inv <- rgamma(1, priors$a.phi1, priors$b.phi1)
    }

    #############
    # update mu #
    #############

    if(K > 1){
      mu <- rnorm(1, (phi2inv*mean(theta[-1]))/((K-1)*phi2inv + priors$sig2inv.mu1),
                  sqrt(1/((K-1)*phi2inv + priors$sig2inv.mu1)))
    }else{
      mu <- rnorm(1, 0, sqrt(1/priors$sig2inv.mu1))
    }


    #############
    # inclusion #
    #############

    # indicator if beta is NOT in the null cluster
    pip.beta[s,] <- as.numeric(beta != 0)

    #################
    # store results #
    #################

    K_keep[s] <- K # number of main effect clusters
    alpha_keep[s] <- alpha # DP1 param
    beta_keep[s,] <- beta # main effects
    mu_keep[s] <- mu # base mean for D1
    phi2inv_keep[s] <- phi2inv # base precision for D1
    if(!is.null(W)) gamma_keep[s,] <- gamma # covariates
    sig2inv_keep[s] <- sig2inv # error
    
  }
  
  # unscale estimates, Y.save is what went into function 
  if(scaleY == TRUE){
    gamma_keep.unscaled <- gamma_keep*sd(Y.save) 
    if(is.null(W)) {
      gamma_keep.unscaled <- matrix(mean(Y.save), length(Y), 1)
    }else {
      gamma_keep.unscaled[,1] <- gamma_keep.unscaled[,1] + mean(Y.save) # overall intercept
    }
    beta_keep.unscaled <- beta_keep*sd(Y.save)
  }else{
    if(is.null(W)){
      gamma_keep.unscaled <- matrix(0, length(Y), 1) # to avoid later problems in summary and predict
    }else{
      gamma_keep.unscaled <- gamma_keep
    }
    beta_keep.unscaled <- beta_keep
  }

  list1 <- list(X = X, Y = Y.save, W = W,
               alpha = alpha_keep[-(1:nburn)],
               beta = beta_keep.unscaled[-(1:nburn),],
               mu = mu_keep[-(1:nburn)],
               phi2inv = phi2inv_keep[-(1:nburn)],
               gamma = gamma_keep.unscaled[-(1:nburn),],
               sig2inv = sig2inv_keep[-(1:nburn)],
               pip.beta = pip.beta[-(1:nburn),])
  
  class(list1) <- "npb"
  
  return(list1)
}
