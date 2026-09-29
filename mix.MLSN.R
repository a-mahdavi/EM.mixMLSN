#==============================================================
# ECM algorithm for finite mixtures of multivariate
# logit-skew-normal distributions (FM-MLSN)
# Sahu-Dey-Branco (SDB) parameterization
#
#==============================================================
 library(mvtnorm)
 library(tmvtnorm)

 mix.MLSN <- function(y, g = 1, pi = NULL, mu = NULL, S = NULL, Lambda = NULL,
                     iter.max = 500, tol = 1e-6, get.init = TRUE, 
                     diagonal.Lambda = TRUE, group = TRUE) {
  
  begin <- proc.time()[3]
  
  #==============================================================
  # Logit function
  #==============================================================
  logit <- function(y) {
    qlogis(y)
  } 
  #==============================================================
  # Multivariate MLSN density 
  #==============================================================
  dMLSN <- function(y, mu, Sigma, Lambda, log = FALSE) {
    if(is.vector(y)) y <- matrix(y, nrow = 1)
    n <- nrow(y)
    p <- ncol(y)
    u <- qlogis(y)
    
    if(!is.matrix(Lambda)) {
      Lambda <- if(is.vector(Lambda)) diag(Lambda) else as.matrix(Lambda)
    }
    
    Omega <- Sigma + tcrossprod(Lambda)
    Omega <- (Omega + t(Omega)) / 2
    Omega_inv <- chol2inv(chol(Omega))
    Delta <- diag(p) - crossprod(Lambda, Omega_inv) %*% Lambda
    Delta <- (Delta + t(Delta)) / 2
    phi1 <- dmvnorm(u, mean = mu, sigma = Omega)
    arg <- crossprod(Lambda, Omega_inv) %*% t(sweep(u, 2, mu))
    Phi <- numeric(n)
    
    for(i in 1:n) {
      Phi[i] <- pmvnorm(lower = rep(-Inf, p), upper = arg[, i], sigma = Delta)[1]
    }
    
    Phi <- pmax(Phi, 1e-300)
    dens <- (2^p) * phi1 * Phi * exp(-rowSums(log(y) + log1p(-y)))
    dens <- pmax(dens, 1e-300)
    
    if(log) return(log(dens)) else return(dens)
  }
  
  #==============================================================
  # Mixture MLSN density
  #==============================================================
  dmixMLSN <- function(y, pi, mu, Sigma, Lambda, log = FALSE) {
    g <- length(pi)
    n <- nrow(y)
    f <- numeric(n)
    
    for(j in 1:g) {
      f <- f + pi[j] * dMLSN(y, mu[[j]], Sigma[[j]], Lambda[[j]])
    }
    
    f[f < 1e-300] <- 1e-300
    if(log) return(log(f)) else return(f)
  }
  
  #==============================================================
  # Data preprocessing
  #==============================================================
  y <- as.matrix(y)
  colnames(y) <- NULL
  n <- nrow(y)
  p <- ncol(y)
  y <- apply(y, 2, function(x) ifelse(x >= 1, 0.999, x))
  y <- apply(y, 2, function(x) ifelse(x <= 0, 0.001, x))
  
  #==============================================================
  # Initialization
  #==============================================================
  if(get.init == TRUE) {
    init <- kmeans(logit(y), g, algorithm = "Hartigan-Wong")
    pi <- init$size / nrow(y)
    mu <- S <- Lambda <- list()
    
    for(j in 1:g) {
      x <- y[init$cluster == j, , drop = FALSE]
      if(nrow(x) == 1) x <- matrix(x, nrow = 1)
      u <- logit(x)
      mu[[j]] <- colMeans(u)
      S[[j]] <- cov(u) + diag(1e-4, p)
      S[[j]] <- (t(S[[j]]) + S[[j]]) / 2
      
      if(diagonal.Lambda) {
        Lambda_j <- rep(0, p)
        for(r in 1:p) {
          u_r <- u[, r]
          m3 <- mean((u_r - mean(u_r))^3)
          m2 <- mean((u_r - mean(u_r))^2)
          if(m2 > 0) {
            Lambda_j[r] <- min(max(0.1 * m3 / (m2^(3/2)), -0.5), 0.5)
          }
        }
        Lambda[[j]] <- diag(Lambda_j)
      } else {
        Lambda[[j]] <- diag(rep(0.05, p))
      }
    }
  }
  
  #==============================================================
  # ECM iterations
  #==============================================================
  LL.old <- 1
  iter <- 0
  conv <- FALSE
  
  while(!conv && iter < iter.max) {
    iter <- iter + 1
    
    #==============================
    # E-step
    #==============================
    dens.comp <- matrix(0, n, g)
    tau.hat <- vector("list", g)
    Psi.hat <- vector("list", g)
    
    for(j in 1:g) {   
      tau.hat[[j]] <- matrix(0, n, p)
      Psi.hat[[j]] <- array(0, c(p, p, n))
      dens.comp[, j] <- pi[j] * dMLSN(y, mu[[j]], S[[j]], Lambda[[j]])
    }
    
    mix.dens <- rowSums(dens.comp)
    z.hat <- dens.comp / mix.dens
    
    #==============================
    # Truncated-normal moments using MomTrunc
    #==============================
    U <- logit(y)
    
    for(j in 1:g) {
      Omega <- S[[j]] + Lambda[[j]] %*% t(Lambda[[j]])
      Omega <- (t(Omega) + Omega) / 2
      Omega.inv <- solve(Omega)
      Sigma.tau <- diag(p) - t(Lambda[[j]]) %*% Omega.inv %*% Lambda[[j]]
      Sigma.tau <- (t(Sigma.tau) + Sigma.tau) / 2
      
      for(i in 1:n) {
        mu.tau <- t(Lambda[[j]]) %*% Omega.inv %*% (U[i, ] - mu[[j]])
        mu.tau <- as.vector(mu.tau)
        mt <- MomTrunc::meanvarTMD(mu = mu.tau, Sigma = Sigma.tau, 
                                   lower = rep(0, p), upper = rep(Inf, p), 
                                   dist = "normal")    
        tau.hat[[j]][i, ] <- mt$mean
        Psi.hat[[j]][, , i] <- mt$EYY
      }
    }
    
    #==============================
    # CM-step 1: mixing proportions
    #==============================
    nj <- colSums(z.hat)
    pi <- nj / n
    
    #==============================
    # CM-step 2: update mu
    #==============================
    for(j in 1:g) {   
      num <- rep(0, p)
      for(i in 1:n) {
        num <- num + z.hat[i, j] * (U[i, ] - Lambda[[j]] %*% tau.hat[[j]][i, ])
      }
      if(nj[j] > 1e-10) {
        mu[[j]] <- num / nj[j]
      }
    }
    
    #==============================
    # CM-step 3: update Lambda
    #==============================
    if(diagonal.Lambda) {  
      for(j in 1:g) {     
        A <- matrix(0, p, p)
        B <- matrix(0, p, p)
        
        for(i in 1:n) {
          ui <- matrix(U[i, ] - mu[[j]], p, 1)
          taui <- matrix(tau.hat[[j]][i, ], p, 1)      
          A <- A + z.hat[i, j] * taui %*% t(ui)
          B <- B + z.hat[i, j] * Psi.hat[[j]][, , i]
        }
        
        Sigma.inv <- solve(S[[j]])   
        lhs <- Sigma.inv * B
        rhs <- Sigma.inv * A
        lambda.f <- (solve(lhs) %*% rhs) %*% rep(1, p)
        Lambda[[j]] <- diag(as.vector(lambda.f))
      } 
    } else {
      for(j in 1:g) {
        #==============================
        # FULL LAMBDA FORMULA - Section 3, Equation (20)
        # Lambda_j = A * B^{-1}
        #==============================
        A <- matrix(0, p, p)
        B <- matrix(0, p, p)
        
        for(i in 1:n) {
          ui <- matrix(U[i, ] - mu[[j]], p, 1)
          taui <- matrix(tau.hat[[j]][i, ], p, 1)      
          A <- A + z.hat[i, j] * ui %*% t(taui)
          B <- B + z.hat[i, j] * Psi.hat[[j]][, , i]
        }
        
        Lambda[[j]] <- A %*% solve(B)
      }
    }
    
    #==============================
    # CM-step 4: update Sigma 
    #==============================
    for(j in 1:g) {    
      Sigma.new <- matrix(0, p, p)
      
      for(i in 1:n) {
        ui <- matrix(U[i, ] - mu[[j]], p, 1)
        taui <- matrix(tau.hat[[j]][i, ], p, 1)
        Psii <- Psi.hat[[j]][, , i]
        Sigma.new <- Sigma.new + z.hat[i, j] * (
          ui %*% t(ui)
          - ui %*% t(taui) %*% t(Lambda[[j]])
          - Lambda[[j]] %*% taui %*% t(ui)
          + Lambda[[j]] %*% Psii %*% t(Lambda[[j]])
        )
      }
      
      S[[j]] <- Sigma.new / nj[j]
      S[[j]] <- (S[[j]] + t(S[[j]])) / 2
    }
    
    # Compute log-likelihood
    LL_new <- sum(log(dmixMLSN(y, pi, mu, S, Lambda)))
    diff <- abs((LL_new-LL.old) / LL.old)   
    
    cat("Iter:", iter, " LogLik:", round(LL_new, 4), 
        " Diff:", round(diff, 6), "\n")
    
    conv <- (diff < tol)
    LL.old <- LL_new
  }
  
  if(diagonal.Lambda) {
    n.par <- (g - 1) + g * p + g * p * (p + 1) / 2 + g * p
  } else {
    n.par <- (g - 1) + g * p + g * p * (p + 1) / 2 + g * p * p
  }
  
  AIC <- -2 * LL_new + 2 * n.par
  BIC <- -2 * LL_new + log(n) * n.par
  EDC <- -2 * LL_new + 0.2 * sqrt(n) * n.par
  
  #==============================
  # Posterior probabilities and classifications
  #==============================
  dens.comp <- matrix(0, n, g)
  
  for(j in 1:g) {
    dens.comp[, j] <- pi[j] * dMLSN(y, mu[[j]], S[[j]], Lambda[[j]])
  }
  
  posterior <- dens.comp / rowSums(dens.comp)
  classification <- apply(posterior, 1, which.max)
  elapsed <- proc.time()[3] - begin
  
  #==============================
  # Output
  #==============================
  fit <- list(
    family = "Finite Mixture MLSN (SDB)",
    g = g,
    pi = pi,
    mu = mu,
    Sigma = S,
    Lambda = Lambda,
    group = classification,
    loglik = LL_new,
    AIC = AIC,
    BIC = BIC,
    EDC = EDC,
    iter = iter,
    converged = conv,
    elapsed = elapsed,
    diagonal.Lambda = diagonal.Lambda
  )
  
  class(fit) <- "mixMLSN"
  return(fit)
 }














