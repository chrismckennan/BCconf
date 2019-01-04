library(MASS)

HubersLoss <- function(Y, X, ind.cov, Sigma, Gamma) {
  SXX.inv <- solve(t(X) %*% X)
  Y1 <- cbind(Y %*% X %*% SXX.inv)

  n <- nrow(X) - ncol(X); p <- nrow(Y); r <- ncol(Gamma)
  tmp.s.theory <- svd(t(Gamma) %*% Gamma * n/p)$d
  ols.s.theory <- t(Gamma/Sigma) %*% Gamma * n/p
  bias.corr.theory <- diag(tmp.s.theory / (tmp.s.theory - mean(Sigma)), nrow=r, ncol=r)
  bias.corr.ols <- solve(ols.s.theory - diag(1, nrow=r, ncol=r)) %*% ols.s.theory
  
  Omega.hub <- matrix(0, nrow=ncol(Gamma), ncol=length(ind.cov))
  Omega.hub2 <- matrix(0, nrow=ncol(Gamma), ncol=length(ind.cov))
  for (i in 1:length(ind.cov)) {
  	out.hub.i <- Pen.Huber(Y1[,ind.cov[i]], Sigma, Gamma)
  	#out.hub2.i <- Pen.Huber(Y1[,ind.cov[i]], Sigma / nrow(X), Gamma)
  	Omega.hub[,i] <- out.hub.i$alpha
  	#Omega.hub2[,i] <- out.hub.i$alpha
  }

  Omega.ols <- solve(t(Gamma / Sigma) %*% Gamma - diag(p/n, nrow=r, ncol=r), t(Gamma / Sigma) %*% Y1[,ind.cov])
  Omega.theory <- solve(t(Gamma)%*%Gamma - diag(p/n * mean(Sigma), nrow=r, ncol=r), t(Gamma)%*%Y1[,ind.cov])
  Omega.naive <- solve(t(Gamma)%*%Gamma, t(Gamma)%*%Y1[,ind.cov])

  return(list(Omega.ols=Omega.ols, Omega.theory=Omega.theory, Omega.hub=Omega.hub, Omega.naive=Omega.naive, out.hub=out.hub.i$out))
}

Pen.Huber <- function(Y, Sigma, Gamma, k.hub = 1.345, max.iter=500, tol=1e-8) {
  Y <- cbind(Y)
  Gamma <- cbind(Gamma)
  alpha <- solve(t(Gamma) %*% Gamma, t(Gamma) %*% Y)
  for (d in 1:ncol(Y)) {
    for (i in 1:max.iter) {
      u <- as.vector(Y[,d] - Gamma %*% alpha[,d]) / sqrt(Sigma)
      w <- MASS::psi.huber(u = u, k = k.hub)
      grad <- as.vector(t(Gamma / sqrt(Sigma)) %*% (w * u))
      if (max(abs(grad)) < tol) {break}
      alpha[,d] <- solve(t(Gamma * (w/Sigma)) %*% Gamma, t(Gamma * (w/Sigma)) %*% Y[,d])
    }
  }
  return(list(alpha=alpha))
}

Pen.Huber0 <- function(Y, Sigma, Gamma, Omega = 0, inv.V.Omega = 0, k.hub = 1.345, ind.use=NULL, q = 0, max.iter=500, tol=1e-8) {
  K <- ncol(Gamma)
  if (q == 0) {
    inv.V.Omega <- matrix(0, nrow=K, ncol=K)
    Omega <- matrix(0, nrow=K, ncol=1)
  }
  if (! is.null(ind.use)) {
    Y <- Y[ind.use]
    Sigma <- Sigma[ind.use]
    if (K > 1) {
      Gamma <- Gamma[ind.use,]
    } else {
      Gamma <- cbind(Gamma[ind.use])
    }
  }
  Gamma <- Gamma / sqrt(Sigma)
  Y <- Y / sqrt(Sigma) - Gamma %*% Omega
  
  alpha.tilde.0 <- solve( t(Gamma) %*% Gamma + q * inv.V.Omega, t(Gamma) %*% Y)       #Penalized least squares as a starting point
  
  for (iter in 1:max.iter) {
    psi.0 <- as.vector(psi.huber(Y - Gamma %*% alpha.tilde.0, k.hub))
    if (K == 1) {
      grad.0 <- q * inv.V.Omega %*% alpha.tilde.0 - sum(Gamma * psi.0)
      diff <- abs(grad.0)
    } else {
      grad.0 <- q * inv.V.Omega %*% alpha.tilde.0 - apply(Gamma * psi.0, 2, sum)
      diff <- max(abs(grad.0))
    }
    if (diff < tol) {
      return(list(alpha=alpha.tilde.0 + Omega, n.iter=iter, out=1, resid=Y - Gamma %*% alpha.tilde.0))
    }
    
    weights.0 <- as.vector(Weights.huber(Y - Gamma %*% alpha.tilde.0, k.hub))
    tmp.mat <- Gamma * weights.0
    alpha.tilde.0 <- solve( t(tmp.mat) %*% Gamma + q * inv.V.Omega, t(tmp.mat) %*% Y )
    if (K == 1) {
      alpha.tilde.0 <- matrix(alpha.tilde.0, nrow=1, ncol=1)
    }
    
  }
  return(list(alpha=alpha.0, n.iter=iter, out=0))
}



psi.huber <- function(x, k) {
  ind.small <- as.numeric(abs(x) <= k)
  return( ind.small * x/k + (1 - ind.small) * sign.vec(x) )
}

Weights.huber <- function(x, k) {
  ind.small <- as.numeric(abs(x) <= k)
  return( ind.small * 1/k + (1-ind.small) / abs(x) )
}

sign.vec <- function(x) {
  return( as.numeric( x > 0 ) - as.numeric( x < 0 ) )
}

rho.huber <- function(x, k) {
  ind.small <- abs(x) <= k
  return( sum(1/2/k * x[ind.small]^2) + sum(abs(x[!ind.small]) - k/2) )
}

fun.value <- function(x, k, diff, inv.V, q) {
  return( (rho.huber(x, k) + q/2 * t(diff) %*% inv.V %*% diff)/p )
}
