library(cate)
library(irlba)

#' Estimate the number of latent factors
#' 
#' Uses the \code{est.factor.num} function from the \code{CATE} package to estimate the number of latent factors
#' 
#' @param Y a p x n data matrix, where p = #of units (i.e. genes) and n = #of samples
#' @param X an n x d model matrix, where d = total number of observed covariates to include in the model. This includes the covariates of interest (i.e. disease status), as well as nuisance covariates like the intercept, plate number, DNA concentration, etc.
#' @param method Either bi-cross validation ("bcv") or eigenvalue distance ("ed"), with "bcv" being the default. See \url{https://cran.r-project.org/web/packages/cate/cate.pdf} for more details.
#' @param max.r Maximum number of latent factors to consider
#' @param nRepeat Number of times to perform bcv
#' @param ... Additional parameters to include in \code{cate::est.factor.num}.
#' 
#' @return A \code{cate::est.factor.num} object.
#' @export
est.conf.num <- function(Y, X, method=c("bcv", "ed"), max.r, nRepeat=20, ...) {
  n <- nrow(X)
  d <- ncol(X)
  
  Q <- qr.Q(qr(X), complete=T)[,(d+1):n]
  r <- cate::est.factor.num(Y = t(Y %*% Q), method = method, rmax = max.r, nRepeat = nRepeat, ...)
  return(r)
}

svd.LandC <- function(Y, X=NULL, r, min.true.eig=0.1, max.iter.svd=1, svd.method="slow") {
  min.true.eig <- max(min.true.eig, 1e-2)
  p <- nrow(Y)
  n <- nrow(X)
  
  if (!is.null(X)) {
    d <- ncol(X)
    Q <- qr.Q(qr(X), complete = T )[,(d+1):n]
    Y.tilde <- Y %*% Q
  } else {
    d <- 0
    Q <- diag(n)
    Y.tilde <- Y
  }
  if (r > 0) {
    if (svd.method == "fast") {
      s <- irlba(A=Y.tilde, nv = r, tol = 1/sqrt(n) * 1e-4)
      C <- sqrt(n - d) * cbind(s$v)
      R <- Y.tilde - Y.tilde %*% C %*% solve(t(C)%*%C, t(C))
      Sigma <- 1/(n-d-r) * rowSums(R^2)
      if (max.iter.svd > 1) {
        for (i in 2:max.iter.svd) {
          s <- irlba(A=Y.tilde/sqrt(Sigma), nv = r, tol = 1/sqrt(n) * 1e-4)
          C <- sqrt(n - d) * cbind(s$v)
          R <- Y.tilde - Y.tilde %*% C %*% solve(t(C)%*%C, t(C))
          Sigma <- 1/(n-d-r) * rowSums(R^2)
        }
      }
    } else {
      s <- svd(Y.tilde, nv = r)
      C <- sqrt(n - d) * cbind(s$v[,1:r])
      R <- Y.tilde - Y.tilde %*% C %*% solve(t(C)%*%C, t(C))
      Sigma <- 1/(n-d-r) * rowSums(R^2)
      if (max.iter.svd > 1) {
        for (i in 2:max.iter.svd) {
          s <- svd(Y.tilde/sqrt(Sigma), nv = r)
          C <- sqrt(n - d) * cbind(s$v)
          R <- Y.tilde - Y.tilde %*% C %*% solve(t(C)%*%C, t(C))
          Sigma <- 1/(n-d-r) * rowSums(R^2)
        }
      }
    }
    L <- Y.tilde %*% C %*% solve(t(C)%*%C)
    s.L <- svd(L); L <- L %*% cbind(s.L$v); C <- C %*% cbind(s.L$v)
    rho <- mean(Sigma)
    lambda.min <- s$d[r]^2 / p / rho - 1
    if (lambda.min < min.true.eig) {
    	warning(paste("Last eigenvalue is too small. Try r.confound = ", as.character(r-1)))
    	return(0)
    }
    
    return(list(L=L, C=Q%*%C, Sigma=Sigma))
  } else {
    Sigma <- 1/(n-d) * rowSums(Y.tilde * Y.tilde)
    return(list(Sigma=Sigma))
  }
}

ml.LandC <- function(Y, X, r) {
  p <- nrow(Y)
  n <- nrow(X)
  d <- ncol(X) 
  
  Q <- qr.Q(qr(X), complete=T)[,(d+1):n]  
  Y.tilde <- Y %*% Q
  out.em <- fa.em(t(Y.tilde), r=r)
  return(list(L=out.em$Gamma, Sigma=out.em$Sigma, iter=out.em$niter))
}
