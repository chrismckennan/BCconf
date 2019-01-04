#####Simulate data with a t-distribution#####
require(cate)

SimulateData <- function(n, p, Info, Omega, L.sd, B.sd, B.pi, mean.Sigma=1, var.Sigma=0.5^2, t.df=Inf, B=NULL, L=NULL, X=NULL, LwithSigma=FALSE) {
  Sim <- list()
  
  K <- length(Info)
  m <- n - 1
  X <- t( qr.Q(qr(cbind(rep(1,n))), complete=T)[,2:n] ) %*% cbind(c( rep(1,n/2), rep(0,n/2) ))
  #multiplier <- sqrt((m-1) / sum(X^2))
  multiplier <- 1
  X <- multiplier * X
  Sim$X <- X
  
  Lambda <- Info / (m-1)
  beta <- mean.Sigma/var.Sigma; alpha <- beta * mean.Sigma
  Sigma <- rgamma(p, alpha, beta)
  Sim$Sigma <- Sigma
  if (t.df == Inf) {
    E <- matrix(rnorm(m*p), nrow=p, ncol=m)
  } else {
    E <- matrix(rt(m*p, df=t.df), nrow=p, ncol=m) * sqrt((t.df-2)/t.df)
  }
  E <- E * sqrt(Sigma)
  
  ##Simulate B##
  if (is.null(B)) {
    B <- rep(0, p); tmp <- rbinom(n = p, size = 1, prob = 1-B.pi)
    B[tmp == 1] <- rnorm(sum(tmp)) * B.sd / multiplier
    Sim$B <- B
  } else {
    Sim$B <- B / multiplier
  }
  
  ##Simulate L##
  Sim$L <- L
  if (is.null(L)) {
    for (k in 1:K) {
      lambda <- Lambda[k]
      delta.L <- lambda / L.sd^2   #Fraction of non-null L components
      if (k == 1) {
        if (delta.L > 1) {
          L <- cbind( rnorm(p) * sqrt(lambda) )
        } else {
          L <- cbind( rnorm(p) * L.sd * rbinom(p, 1, delta.L) )
        }
      } else {
        if (delta.L > 1) {
          L <- cbind(L, rnorm(p) * sqrt(lambda))
        } else {
          L <- cbind( L, rnorm(p) * L.sd * rbinom(p, 1, delta.L) )
        }
      }
    }
    if (LwithSigma) {
      L <- L * sqrt(Sigma)
    }
    svd.LtL <- svd(1/p * t(L)%*%L)
    L <- L %*% svd.LtL$u %*% diag( sqrt(Lambda/svd.LtL$d), nrow=K, ncol=K )
    Sim$L <- L
  }
    
  ##Simulate C##
  C.perp <- sqrt(m-1) * qr.Q(qr(cbind(X)), complete=T)[,2:m] %*% svd(matrix(rnorm((m-1)*K), nrow=m-1, ncol=K))$u
  C.ImX <- X %*% rbind(Omega + sqrt(1/sum(X^2)) * rnorm(K))
  C <- C.perp + C.ImX
  Sim$C <- C
  
  ##Simulate data##
  Y <- cbind(B) %*% t(X) + L %*% t(C) + E
  Sim$Y <- Y
  
  return(Sim)
}
