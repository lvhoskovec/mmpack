#' Fit Bayesian profile regression 
#'
#' Fits the Bayesian profile regression model using MCMC methods 
#'
#' @param niter number of total iterations
#' @param nburn number of burn-in iterations
#' @param X n by p matrix of predictor data
#' @param Y n by 1 vector of continuous response data
#' @param W n by q matrix of covariate data
#' @param C maximum number of clusters allowed
#' @param scaleY logical; if TRUE response will be centered and scaled before model fit
#' @param DPgamma logical; if TRUE (default) alpha has a gamma prior, else alpha has Unif(.03, 10) prior
#' @param varsel logical; if TRUE binary cluster variable selection is implemented (see Chung and Dunson 2009), default is FALSE
#' @param priors list of prior hyperparameters, see package documentation for details 
#' @param sup logical; if TRUE (default) fits supervised model, else fits unsupervised model
#'
#' @importFrom stats pgamma qgamma rWishart rbinom rgamma
#' @importFrom mvtnorm dmvnorm 
#' @importFrom matrixcalc is.positive.definite
#'
#'
#' @return list with components
#' \itemize{
#'    \item X: predictor data matrix
#'    \item W: covariate data matrix
#'    \item Y: response data vector
#'    \item C: maximum number of allowable clusters
#'    \item alpha: DP parameter estimates
#'    \item mu: array of cluster means estimates
#'    \item psi: cluster weights estimates
#'    \item Z: cluster indicators at each iteraction
#'    \item delta: regression coefficient estimates for theta (risk) and gamma (fixed effects)
#'    \item sig2inv: error precision estimates
#'    \item kap2inv: cluster intercept (risk) precision estimates
#'    \item phi2inv: fixed effect precision estimates
#'    \item rho: rho estimates
#' }
#' @export


# could probably remove DPgamma from the package...

profileReg <- function(niter, nburn, X, Y, W, C = 20,
                       scaleY = FALSE, DPgamma = TRUE,
                       varsel = FALSE,
                       priors, 
                       sup = TRUE){
  
  if(nburn >= niter) stop("Number of iterations (niter) must be greater than number of burn-in iteractions (nburn)")
  
  
    ##############
    ### priors ###
    ##############
    X <- as.matrix(X)
    W <- as.matrix(W)
  
  
    if(missing(priors)) priors <- NULL
    if(is.null(priors$alpha.sig)) priors$alpha.sig <- 2.5 # shape parameter for gamma prior on sig2inv
    if(is.null(priors$beta.sig)) priors$beta.sig <- 2.5 # rate parameter for gamma prior on sig2inv
    if(is.null(priors$alpha.kap)) priors$alpha.kap <- 7/2 # shape parameter for gamma prior on kap2inv, random precision for thetas
    if(is.null(priors$beta.kap)) priors$beta.kap <- 43.75/2 # rate parameter for gamma prior on kap2inv, random precision for thetas
    if(is.null(priors$alpha.phi)) priors$alpha.phi <- 7/2 # shape parameter for gamma prior on phi2inv, random precision for fixed effects
    if(is.null(priors$beta.phi)) priors$beta.phi <- 43.75/2 # rate parameter for gamma prior on phi2inv, random precision for fixed effects
    if(is.null(priors$alpha.alpha)) priors$alpha.alpha <- 2 # shape parameter for gamma prior on alpha
    if(is.null(priors$beta.alpha)) priors$beta.alpha <- 1 # rate parameter for gamma prior on alpha
    if(is.null(priors$nu)) priors$nu <- colMeans(X) # mean parameter for normal prior on mu, vector of exposure profile means
    if(is.null(priors$R)) priors$R <- (1/ncol(X))*(chol2inv(chol(var(X)))) # scale matrix hyperparameter for Wishart prior in cluster precision matrix
    if(is.null(priors$r)) priors$r <- ncol(X) # degrees of freedom parameter for Wishart prior on cluster precisions
    if(is.null(priors$alpha.rho)) priors$alpha.rho <- .5 # shape1 parameter on beta dist for rho
    if(is.null(priors$beta.rho)) priors$beta.rho <- .5 # shape2 parameter on beta dist for rho
    
    # Lambda is hyperparameter covariance matrix for multivariate mu
    # reflects prior knowledge on how correlated the pollutants are
    # if Lambda is not specified, set it to a diagonal matrix
    
    p <- ncol(X) # number of predictor variables
    if(is.null(priors$Lambda)) {
      priors$Lambda <- apply(X, 2, function(X) (max(X) - min(X))^2)
      if (p > 1) priors$Lambda <- diag(priors$Lambda)
    } else {
      if (!isSymmetric(priors$Lambda) | !is.positive.definite(priors$Lambda))
        stop("Lambda must be a symmetric positive definite matrix")
    }
    
    LamInv <- chol2inv(chol(priors$Lambda)) # precision matrix 
    
    if(!isSymmetric(priors$R) | !is.positive.definite(priors$R)){
      stop("scale matrix R must be symmetric and positive definite")
    }else{
      Rinv <- chol2inv(chol(priors$R)) # inverse scale parameter for Wishart distribution
    }
    
    ##############
    ### params ###
    ##############
    
    if(scaleY == TRUE){
      Y.save <- Y
      Y <- scale(Y)
    }else{
      Y.save <- Y
    }
    
    n <- length(Y) # sample size
    pw <- ncol(W) # number of covariates
    
    #######################
    ### starting values ###
    #######################
    
    alpha <- 1 # DP parameter
    psi <- rep(1/C, C) # cluster weights
    V <- rbeta(C, 1, alpha) # conditional weights
    V[C] <- 1 # apply DP truncation to C classes
    
    Z <- rcat(matrix(runif(n*C),n,C)) # initial category assignment
    loglik_c <- matrix(NA,n,C) # log-likelihood of being in each cluster
    
    nc <- unlist(lapply(seq(1:C), FUN = function(c) length(which(Z==c))))
    mu <- matrix(rnorm(C*p), C, p) # cluster means
    
    SigInv <- array(diag(p), c(p, p, C)) # cluster precisions
    gamma <- rnorm(pw) # regression coefficients for covariates
    theta <- rnorm(C) # cluster-specific parameters
    delta <- c(theta, gamma) # block parameters
    sig2inv <- 1 # error term precision
    kap2inv <- rep(1, C) # theta precisions
    phi2inv <- rep(1, pw) # gamma precisions
    tau2inv <- c(kap2inv, phi2inv) # precisions for theta and gamma, theta[c] exists for each cluster
    
    # variable selection staring values 
    
    pi.vs <- matrix(1,C,p, byrow = T)
    p.w <- .5 # prior selection probability
    rho <- rep(.5, p)
    PI <- array(NA, dim = c(p, p, C))
    for (c in 1:C){
      PI[,,c] <- diag(pi.vs[c,])
    }
    
    mu.0 <- matrix(0, C, p) # starting values for mu_c
    omega <- rep(1,p)
    
    #####################
    ### Storage Space ###
    #####################
    
    alpha_keep <- rep(NA, niter)
    mu_keep <- array(NA, c(niter, C, p))
    SigInv_keep <- array(NA, c(niter, p, p, C))
    psi_keep <- matrix(NA, niter, C)
    Z_keep <- matrix(NA, niter, n)
    delta_keep <- matrix(NA, niter, C + pw)
    sig2inv_keep <- rep(NA, niter)
    kap2inv_keep <- matrix(NA, niter, C)
    phi2inv_keep <- matrix(NA, niter, pw)
    pi.vs_keep <- matrix(NA, niter, p)
    rho_keep <- matrix(NA, niter, p)
    # omega_keep <- matrix(NA, niter, p)
    
    ############
    ### MCMC ###
    ############
    
    for (s in 1:niter){
      
      ##################################
      ### update cluster assignments ###
      ##################################
      
      if(sup == FALSE){
        for (c in 1:C){
          tryCatch({
            loglik_c[,c] <- log(psi[c]) +
              dmvnorm(X, mu[c,], chol2inv(chol(SigInv[,,c])), log = T) }, 
            error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        }
      }else{
        for(c in 1:C){
          tryCatch({
            loglik_c[,c] <- log(psi[c]) +
              dmvnorm(X, mu[c,], chol2inv(chol(SigInv[,,c])), log = T) +
              dnorm(Y, delta[c] + W %*% delta[C + 1:pw],
                    sqrt(1/sig2inv), log = T) }, 
            error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        }
      }
      
      # calculate likelihood
      lik_c <- exp(loglik_c - apply(loglik_c, 1, logsum))
      
      # update cluster assignment
      Z <- rcat(lik_c) 
      nc <- unlist(lapply(seq(1:C), FUN = function(c) length(which(Z==c))))
      
      
      ##############################
      ### update cluster weights ###
      ##############################
      
      # update conditional weights for stick-breaking process
      for (c in 1:(C-1)) V[c] <- rbeta(1, nc[c] + 1, alpha + sum(nc[(c+1):C]))
      
      # handle overflow
      V[which(V == 1)] <- 1-10e-8
      V[C] <- 1
      
      # update mixing weights 
      psi[1] <- V[1]
      cumV <- cumprod(1-V)
      for (c in 2:C) psi[c] <- V[c]*cumV[c-1]
      
      
      ############################################
      ### update mu and Sigma for each cluster ###
      ############################################
      
      if (varsel == FALSE){
        
        for(c in 1:C){
          
          # update mu
          v <- chol2inv(chol(nc[c] * SigInv[,,c] + LamInv)) 
          m <- v %*% (SigInv[,,c] %*% colSums(matrix(X[Z==c,], ncol = ncol(X))) + LamInv %*% priors$nu) 
          mu[c,] <- m + t(chol(v)) %*% rnorm(p)
          
          # update SigInv
          M <- Rinv # iteratively sum matrices to get new scale matrix
          for(i in which(Z==c)) {
            M <- M + matrix(apply(matrix(X[i,], ncol = ncol(X)), 1, 
                                  FUN = function(x) tcrossprod(x - mu[c,])), ncol = ncol(X), byrow = TRUE)
          }
          
          SigInv[,,c] <- matrix(rWishart(1, nc[c] + priors$r, chol2inv(chol( M ))),
                                ncol(X), ncol(X))
          
          
        }
        
      }else{
        
        # empty clusters: update from prior 
        for(c in which(nc == 0)){
          
          # update pi.vs
          pi.vs[c,] <- 0 # no variables important for membership into an empty cluster
          PI[,,c] <- diag(pi.vs[c,])
          
          # update mu from prior
          mu.0[c,] <- priors$nu + t(chol(priors$Lambda)) %*% rnorm(p) 
          mu[c,] <- mu.0[c,]
          
          # update SigInv from prior
          SigInv[,,c] <- matrix(rWishart(1, priors$r, priors$R), p, p)
          
        }
        
        # active clusters: update from posterior 
        for (c in which(nc>0)){
          
          # calculate rstar
          # integrate out mu to update pi.vs
          sig2 <- diag(chol2inv(chol(SigInv[,,c])))
          
          lam0 <- diag(priors$Lambda)
          g <- nc[c]/sig2 + 1/lam0
          
          log.gam1 <- (-nc[c]/2)*log(2*pi) - (nc[c]/2)*log(sig2) + log(rho) -.5*log(lam0)-
            
            (apply(matrix(X[Z==c,], ncol = p),2,FUN = function(x) x%*%x)/sig2 + priors$nu^2/lam0 )/2 +
            
            (( apply(matrix(X[Z==c,], ncol = p),2,sum)/sig2 + priors$nu/lam0 )^2/g)/2 +
            
            log(g^(-1/2))
          
          # this is correct, but try to make it faster  
          log.gam0.save <- rep(NA, p) 
          
          for(j in 1:p){
            log.gam0.save[j] <- sum(dnorm(X[Z==c,j], colMeans(X)[j], sqrt(sig2)[j], log = TRUE))
          }
          
          log.gam0 <- log.gam0.save
          
          rstar <- 1 - 1/(1 + exp(log.gam1 - log.gam0))
          
          # update pi.vs 
          pi.vs[c,] <- rbinom(p,1,rstar) 
          pi.vs[c, which(omega == 0)] <- 0 # sparsity inducing prior
          
          PI[,,c] <- diag(pi.vs[c,]) 
          
          # covariate means for all individuals in cluster c
          xbar.c <- colMeans(matrix(X[Z==c,], ncol = p))
          
          # update mu.0 new update
          Sigma.tilde <- chol2inv(chol(LamInv + nc[c] * PI[,,c] %*% SigInv[,,c] %*% PI[,,c])   )
          mu.tilde <- Sigma.tilde %*% ( LamInv %*% priors$nu + nc[c] * PI[,,c] %*% SigInv[,,c]
                                        %*% ( xbar.c - ( (diag(p) - PI[,,c]) %*% colMeans(X) ) ) )
          mu.0[c,] <- mu.tilde + t(chol(Sigma.tilde)) %*% rnorm(p)
          
          mu[c,] <- pi.vs[c,] * mu.0[c,] + (1 - pi.vs[c,]) * colMeans(X)
          
          # update SigInv
          M <- Rinv # adding matrices together, take the sum iteratively to get new scale matrix
          for(i in which(Z==c)) {
            M <- M + matrix(apply(matrix(X[i,], ncol = ncol(X)), 1, 
                                  FUN = function(x) tcrossprod(x - mu[c,])), ncol = ncol(X), byrow = TRUE)
          }
          
          SigInv[,,c] <- matrix(rWishart(1, nc[c] + priors$r, chol2inv(chol( M ))), p, p)
          
          
        } # end for active clusters
        
        
        # update global variable selection parameters: omega then rho
        # this needs to be after the update for pi.vs
        C.star <- length(unique(Z))
        gam <- ((gamma(priors$alpha.rho + priors$beta.rho)*gamma(priors$beta.rho + C.star))/
                  (gamma(priors$beta.rho)*gamma(priors$alpha.rho + priors$beta.rho + C.star)))
        pw.star <- gam*p.w / ( gam*p.w + (1 - p.w)) 
        
        # update omega
        omega[which(apply(pi.vs,2,sum)>0)] <- 1
        omega[which(apply(pi.vs,2,sum)==0)] <- rbinom(length(which(apply(pi.vs,2,sum)==0)), 1, pw.star)
        
        # update rho
        rho[which(omega == 0)] <- 0 # option to get rid of omega, no sparsity inducing prior
        ncp <- colSums(pi.vs) # number of switches that equal 1 for each pollutant
        rho[which(omega == 1)] <- rbeta(length(which(omega==1)), ncp + priors$alpha.rho, C.star - ncp + priors$beta.rho)
        
        
      }# end if varsel = TRUE                           
      
      
      ###################################
      ### update alpha based on prior ###
      ###################################
      
      if (DPgamma == TRUE){
        # alpha has a gamma prior
        alpha <- rgamma(1, priors$alpha.alpha + C - 1, 
                        priors$beta.alpha - sum(log(1-V[1:(C-1)])))
        
      }else{
        # alpha has a uniform prior
        alpha <- qgamma(runif(1, pgamma(0.3, C, sum(-log(1-V[1:(C-1)]))),
                              pgamma(10, C, sum(-log(1-V[1:(C-1)]))))
                        , C, sum(-log(1-V[1:(C-1)])))
      }
      
      
      #########################################
      ### update theta and gamma as a block ###
      #########################################
      
      # design matrix for cluster indicators
      Zstar <- matrix(0, n, C)
      for (c in 1:C){
        Zstar[which(Z == c), c] <- 1
      }
      
      A <- cbind(Zstar, W)
      
      v <- chol2inv(chol(sig2inv * t(A)%*%A + diag(tau2inv)))
      m <- v %*% (sig2inv * t(A)%*%Y)
      delta <- m + t(chol(v)) %*% rnorm(ncol(A))
      
      
      ######################
      ### update tau2inv ###
      ######################
      
      # precision for gamma
      phi2inv <- rgamma(pw, priors$alpha.phi + .5, 
                              priors$beta.phi + .5*(delta[C+1:pw]^2))
      # precision for theta
      kap2inv <- rgamma(C, priors$alpha.kap + .5, priors$beta.kap + .5*delta[1:C]^2)
      
      # # precision for theta
      # if (tprior == FALSE){
      #   # all from same distribution (normal prior)
      #   kap2inv <- rgamma(C, priors$alpha.kap + C/2, priors$beta.kap + .5*sum((delta[1:C])^2))
      # }else {
      #   # each from different distribution (t prior)
      #   kap2inv <- rgamma(C, priors$alpha.kap + .5, priors$beta.kap + .5*delta[1:C]^2) # mean 0
      # }

      tau2inv <- c(kap2inv, phi2inv)
      
      ######################
      ### update sig2inv ###
      ######################
      
      sig2inv <- rgamma(1, priors$alpha.sig + n/2, priors$beta.sig + .5*sum((Y - A%*%delta)^2))
      
      #####################
      ### store results ###
      #####################
      
      alpha_keep[s] <- alpha
      mu_keep[s,,] <- mu
      SigInv_keep[s,,,] <- SigInv[,,]
      psi_keep[s,] <- psi
      Z_keep[s,] <- Z
      delta_keep[s,] <- delta
      sig2inv_keep[s] <- sig2inv
      kap2inv_keep[s,] <- tau2inv[1:C]
      # what would this be anyway? sum or mean? these are cluster-specific but PIPs are model-specific
      #pi.vs_keep[s,] <- as.numeric(apply(pi.vs, 2, FUN = function(x) sum(x) > 0 ))
      phi2inv_keep[s,] <- phi2inv
      rho_keep[s,] <- rho
      # omega_keep[s,] <- omega # don't need to save this
      
      
    }
    
    list1 <- list(X = X, W = W, Y = Y.save, C = C,
                  scaleY = scaleY,
                  alpha = alpha_keep[-(1:nburn)],
                  mu = mu_keep[-(1:nburn),,],
                  SigInv = SigInv_keep[-(1:nburn),,,],
                  psi = psi_keep[-(1:nburn),],
                  Z = Z_keep[-(1:nburn),],
                  delta = delta_keep[-(1:nburn),],
                  sig2inv = sig2inv_keep[-(1:nburn)],
                  kap2inv = kap2inv_keep[-(1:nburn),],
                  phi2inv = phi2inv_keep[-(1:nburn),],
                  rho = rho_keep[-(1:nburn),])
    
    class(list1) <- "bpr"
    
    return(list1)
    
  }
