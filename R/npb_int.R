#' Fit nonparametric Bayes shrinkage model with interactions
#'
#' @param niter number of iterations
#' @param nburn number of burn-in iterations
#' @param X n by p matrix of predictor data
#' @param Y n by 1 vector of response data
#' @param W n by q matrix of covariate data 
#' @param scaleY logical; if TRUE response will be centered and scaled before model fit
#' @param priors list of prior hyperparameters
#' @param intercept logical; indicates if an overall intercept should be estimated with covariates
#' @param XWinteract logical; indiciates in X and W can interact
#' 
#' @importFrom stats var rnorm rbeta runif rgamma
#' @importFrom utils combn
#' 
#' @return list of model estimates, see npb documentation

npb_int <- function(niter, nburn, X, Y, W, scaleY = FALSE, priors, intercept = TRUE,
                    XWinteract = FALSE){
  
  if(nburn >= niter) stop("Number of iterations (niter) must be greater than number of burn-in iteractions (nburn)")
  

  ##########
  # Priors #
  ##########
  
  if(missing(priors)) priors <- NULL
  if(is.null(priors$a.sig)) priors$a.sig <- 1 # shape param for gamma prior on sig2inv
  if(is.null(priors$b.sig)) priors$b.sig <- 1 # rate param for gamma prior on sig2inv
  if(is.null(priors$a.phi1)) priors$a.phi1 <- 1 # shape param for gamma prior on phi2inv, precision of base distribution for main effects
  if(is.null(priors$b.phi1)) priors$b.phi1 <- 1 # rate param for gamma prior on phi2inv, precision of base distribution for main effects
  if(is.null(priors$a.phi2)) priors$a.phi2 <- 1 # shape param for gamma prior on phi2inv.2, precision of base distribution for interactions
  if(is.null(priors$b.phi2)) priors$b.phi2 <- 1 # rate param for gamma prior on phi2inv.2, precision of base distribution for interactions
  if(is.null(priors$alpha.a)) priors$alpha.a <- 2 # shape param for gamma prior on alpha, DP parameter for main effects
  if(is.null(priors$alpha.b)) priors$alpha.b <- 1 # rate param for gamma prior on alpha, DP parameter for main effects
  if(is.null(priors$alpha.2.a)) priors$alpha.2.a <- 2 # shape param for gamma prior on alpha.2, DP parameter for interactions
  if(is.null(priors$alpha.2.b)) priors$alpha.2.b <- 1 # rate param for gamma prior on alpha.2, DP parameter for interactions
  # intercept
  if(is.null(priors$mu.0)) priors$mu.0 <- 0 # mean parameter for normal prior on gamma_0, intercept
  if(is.null(priors$kappa2inv.0)) priors$kappa2inv.0 <- 1 # precision parameter for normal prior on gamma_0, intercept
  # covariates
  if(is.null(priors$mu.gamma) & !is.null(W)) priors$mu.gamma <- rep(0, ncol(W)) # mean vector for normal prior on covariates
  if(is.null(priors$kap2inv)  & !is.null(W)) priors$kap2inv <- rep(1, ncol(W)) # precision vector for normal prior on covariates 
  # # #
  if(is.null(priors$sig2inv.mu1)) priors$sig2inv.mu1 <- 1 # precision param for normal prior on mean of base distribution for main effects
  if(is.null(priors$sig2inv.mu2)) priors$sig2inv.mu2 <- 1 # precision param for normal prior on mean of base distribution for interactions
  if(is.null(priors$alpha.pi)) priors$alpha.pi <- 1 # shape1 parameter for beta prior on pi0 for main effects
  if(is.null(priors$beta.pi)) priors$beta.pi <- 1 # shape2 parameter for beta prior on pi0 for main effects
  if(is.null(priors$alpha.pi2)) priors$alpha.pi2 <- 9 # shape1 parameter for beta prior on pi0 for interactions
  if(is.null(priors$beta.pi2)) priors$beta.pi2 <- 1 # shape2 parameter for beta prior on pi0 for interactions
  
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
    q <- ncol(W) # number of covariates, including overall intercept if there is one 
    gamma <- rnorm(q) # regression coefficients 
  }else{
    q <- 0
    gamma <- NULL
  }
  
  ################
  # interactions #
  ################

  ints <- combn(1:ncol(X), 2)
  Z <- apply(ints, 2, FUN = function(x) {
    X[,x[1]] * X[,x[2]]
  })
  
  ###################################################
  # if interaction between covariates and exposures #
  ###################################################
  
  # include interaction between exposures X and covariates W
  if(XWinteract){
    if(intercept){
      q.star <- q-1
    }else{
      q.star <- q
    }
    XWint <- numeric()
    for(j in 1:p){
      XWint <- cbind(XWint,   apply(W[,-1], 2, FUN = function(wcol) wcol*X[,j]))
    }
    Z <- cbind(Z, XWint)
  }

  r <- ncol(Z) # number of interactions
  zeta <- rnorm(r) # regression coefficients
  psi <- unique(c(0, zeta)) # psi[1] = 0 for null group
  M <- length(psi) # unique clusters

  Q <- rep(NA, r) # allocation variable
  for (j in 1:r){
    Q[j] <- which(psi == zeta[j])
  }

  phi2inv.2 <- 1 # precision for base distribution G2 
  mu.2 <- mean(Z %*% zeta) # mean for base distribution G2
  alpha.2 <- 1 # DP parameter for interactions

  #########
  # error #
  #########

  sig2inv <- 1 # error precision
  
  #####################
  ### Storage space ###
  #####################

  K_keep <- rep(NA, niter) # number of unique main effect clusters
  M_keep <- rep(NA, niter) # number of unique interaction clusters
  alpha_keep <- matrix(NA, niter) # DP for main effects
  alpha.2_keep <- matrix(NA, niter) # DP for interactions
  beta_keep <- matrix(NA, niter, p) # main effects
  zeta_keep <- matrix(NA, niter, r) # interactions
  mu_keep <- rep(NA, niter) # mean of D1
  mu.2_keep <- rep(NA, niter) # mean of D2
  phi2inv_keep <- rep(NA, niter) # precision of D1
  phi2inv.2_keep <- rep(NA, niter) # precision of D2
  if(!is.null(W)) gamma_keep <- matrix(NA, niter, q) # covariates
  sig2inv_keep <- rep(NA, niter) # error precison 
  pip.beta <- matrix(NA, niter, p) # beta PIPs
  pip.zeta <- matrix(NA, niter, r) # zeta PIPs

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

    ##################
    # update alpha.2 #
    ##################

    eta.2 <- rbeta(1, alpha.2 + 1, r)
    pi.eta.2 <- (priors$alpha.2.a + M - 1)/(r*(priors$alpha.2.b - log(eta.2)))/
      ((priors$alpha.2.a + M - 1)/(r*(priors$alpha.2.b - log(eta.2))) + 1)
    if (runif(1) < pi.eta.2){
      alpha.2 <- rgamma(1, priors$alpha.2.a + M, priors$alpha.2.b - log(eta.2))
    }else{
      alpha.2 <- rgamma(1, priors$alpha.2.a + M - 1, priors$alpha.2.b - log(eta.2))
    }

    ##################################################
    # update S: allocation variable for main effects #
    ##################################################

    up.S <- update_S(S = S, p = p, alpha.pi = priors$alpha.pi, beta.pi = priors$beta.pi,
                     beta = beta, theta = theta, zeta = zeta,
                     Y = Y, W = W, gamma = gamma,
                     X = X, Z = Z, sig2inv = sig2inv,
                     phi2inv = phi2inv, mu = mu, alpha = alpha, iter = s)


    S <- up.S$S
    beta <- up.S$beta
    theta <- up.S$theta
    K <- length(theta)

    ##################################################
    # update Q: allocation variable for interactions #
    ##################################################

    up.Q <- update_S(S = Q, p = r, alpha.pi = priors$alpha.pi2, beta.pi = priors$beta.pi2,
                     beta = zeta, theta = psi, zeta = beta,
                     Y = Y, W = W, gamma = gamma,
                     X = Z, Z = X, sig2inv = sig2inv,
                     phi2inv = phi2inv.2, mu = mu.2, alpha = alpha.2, iter = s)


    Q <- up.Q$S
    zeta <- up.Q$beta
    psi <- up.Q$theta
    M <- length(psi)

    ######################################
    # block update theta, psi, and gamma #
    ######################################

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

    if(length(psi)>1){
      
      # reparameterize interactions to remove nulls
      # create the matrix T2 indicating to which psi each zeta belongs 
      # first column for psi = 0
      # rows for zeta, columns for psi
      # a single 1 in each row

      T2 <- matrix(0, r, M)
      for(c in 1:M){
        T2[which(Q==c),c] <- 1
      }

      # clear out the empty clusters/columns in T2 and the null cluster
      T2 <- as.matrix(T2[which(zeta!=0),-1])
      # clear out the columns in Z with null zetas
      Z0 <- as.matrix(Z[, which(zeta != 0)])
      Tpsi <- T2 %*% psi[-1]
      Zzeta <- Z0 %*% Tpsi # equivalent to full matrix Z %*% zeta
      ZT <- Z0 %*% T2 # n by length(psi)-1 design matrix for updating psi

      # add to precision vector
      precision.vector <- c(precision.vector, rep(phi2inv.2, M-1))

    }else {
      # if length(psi) = 1, all zetas are 0
      Zzeta <- 0
      ZT <- NULL
      # add nothing to precision vector
    }

    # add covariate precisions to precision.vector
    precision.vector <- c(precision.vector, priors$kap2inv)

    # update delta = c(theta, psi, gamma)
    A <- cbind(XT, ZT, W)
    if(!is.null(A)){
      m.star <- c(rep(mu, K-1), rep(mu.2, M-1), priors$mu.gamma)
      v <- chol2inv(chol(sig2inv * crossprod(A) + 
                           diag(x=precision.vector, nrow = length(precision.vector))))
      m <- v %*% (sig2inv * (t(A)%*%Y) + precision.vector * m.star)
      delta <- drop(m + t(chol(v)) %*% rnorm(K+M+q-2))
      gamma <- delta[(K+M-1):(K+M+q-2)]
      Adelta <- A %*% delta
    }else{
      Adelta <- 0 # no covariates, no intercept, all regression coefficients are 0 
    }


    if(K > 1){
      theta <- delta[1:(K-1)]
      theta <- c(0, theta)
    }
    if(M > 1){
      psi <- delta[K:(K+M-2)]
      psi <- c(0, psi)
    }


    # update beta deterministically
    beta <- theta[S]

    # update zeta deterministically
    zeta <- psi[Q]

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

    ####################
    # update phi2inv.2 #
    ####################

    if(M > 1){
      phi2inv.2 <- rgamma(1, priors$a.phi2 + (M-1)/2,
                          priors$b.phi2 + .5*sum( (psi[-1] - mu.2)^2 ) )
    }else{
      phi2inv.2 <- rgamma(1, priors$a.phi2, priors$b.phi2)
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

    ###############
    # update mu.2 #
    ###############

    if(M > 1){
      mu.2 <- rnorm(1, (phi2inv.2*mean(psi[-1]))/((M-1)*phi2inv.2 + priors$sig2inv.mu2),
                    sqrt(1/((M-1)*phi2inv.2 + priors$sig2inv.mu2)))
    }else{
      mu.2 <- rnorm(1, 0, sqrt(1/priors$sig2inv.mu2))
    }

    #############
    # inclusion #
    #############

    # indicator if betas/zetas NOT in the null cluster
    pip.beta[s,] <- as.numeric(beta != 0)
    pip.zeta[s,] <- as.numeric(zeta != 0)

    #################
    # store results #
    #################

    K_keep[s] <- K # number of main effect clusters
    M_keep[s] <- M # number of interaction clusters
    alpha_keep[s] <- alpha # DP1 param
    alpha.2_keep[s] <- alpha.2 # DP2 param
    beta_keep[s,] <- beta # main effects
    zeta_keep[s,] <- zeta # interactions
    mu_keep[s] <- mu # base mean for D1
    mu.2_keep[s] <- mu.2 # base mean for D2
    phi2inv_keep[s] <- phi2inv # base precision for D1
    phi2inv.2_keep[s] <- phi2inv.2 # base precision for D2
    gamma_keep[s,] <- gamma # covariates
    sig2inv_keep[s] <- sig2inv # error

  }
  
  # unscale the estimates for predict and summary functions
  if(scaleY == TRUE){
    gamma_keep.unscaled <- gamma_keep*sd(Y.save) 
    if(is.null(W)) {
      gamma_keep.unscaled <- matrix(mean(Y.save), length(Y), 1)
    }else {
      gamma_keep.unscaled[,1] <- gamma_keep.unscaled[,1] + mean(Y.save) # overall intercept
    }
    beta_keep.unscaled <- beta_keep*sd(Y.save)
    zeta_keep.unscaled <- zeta_keep*sd(Y.save)
  }else{
    if(is.null(W)){
      gamma_keep.unscaled <- matrix(0, length(Y), 1) # to avoid later problems in summary and predict
    }else{
      gamma_keep.unscaled <- gamma_keep
    }
    beta_keep.unscaled <- beta_keep
    zeta_keep.unscaled <- zeta_keep
  }

  list1 <- list(X = X, Y = Y.save, Z = Z, W = W, 
               alpha = alpha_keep[-(1:nburn)],
               alpha.2 = alpha.2_keep[-(1:nburn)],
               beta = beta_keep.unscaled[-(1:nburn),], 
               zeta = zeta_keep.unscaled[-(1:nburn),], 
               mu = mu_keep[-(1:nburn)],
               mu.2 = mu.2_keep[-(1:nburn)],
               phi2inv = phi2inv_keep[-(1:nburn)],
               phi2inv.2 = phi2inv.2_keep[-(1:nburn)],
               gamma = gamma_keep.unscaled[-(1:nburn),], 
               sig2inv = sig2inv_keep[-(1:nburn)],
               pip.beta = pip.beta[-(1:nburn),], 
               pip.zeta = pip.zeta[-(1:nburn),])
  
  class(list1) <- "npb"
  
  return(list1)
}

