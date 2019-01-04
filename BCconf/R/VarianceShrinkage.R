#In this script I shrink the variances#
my.eb <- function(Sigma, n, tol=1e-8, max.iter=1e4) {
  s <- Sigma * n
  p <- length(Sigma)
  mu <- mean(1/Sigma)
  beta.0 <- mu / var(1/Sigma)
  alpha.0 <- mu * beta.0
  
  rm(Sigma)
  H.0 <- matrix(0, nrow=2, ncol=2)

  for (i in 1:max.iter) {
    grad.0 <- c( log(beta.0) - mean(log(beta.0 + 1/2*s)) + digamma(alpha.0+n/2) - digamma(alpha.0), alpha.0/beta.0 - (alpha.0 + n/2) * mean( 1/(beta.0 + s/2) ) )
    if (sqrt(sum(grad.0 * grad.0)) * p <= tol) {
      return(list( alpha=alpha.0, beta=beta.0, out=1, iter=i-1 ))
    }
    H.0[1,1] <- trigamma(alpha.0 + n/2) - trigamma(alpha.0)
    H.0[2,2] <- -alpha.0 / beta.0^2 + (alpha.0 + n/2) * mean(1/(beta.0 + s/2)^2)
    H.0[1,2] <- 1/beta.0 - mean(1/(beta.0 + s/2))
    H.0[2,1] <- H.0[1,2]
    theta <- c(alpha.0, beta.0) - solve(H.0, grad.0)
    
    alpha.0 <- theta[1]; beta.0 <- theta[2]
  }
}