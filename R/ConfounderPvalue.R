#' Determine the significance of Omega.ols
#'
#' This tests is the null hypothesis Omega.ols = 0. If you accept the null hypothesis DO NOT replace Omega.ols with 0 in your estimates for the main effect B. Omega.ols is still needed to account for the correlation amongst test statistics.
#'
#' @param Corr.Object A list outputed by \code{Correction}
#' @param Omega.ols Tests the null hypothesis for this specified value of Omega.ols. The default is to use the estimate stored in \code{Corr.Object}
#' 
#' @return A list \item{p.values}{A \code{length(ind.cov)} vector of P values} \item{proxy.n}{A proxy for the sample size}
#' 
#' @export
Sig.Omega <- function(Corr.Object, Omega.ols=NULL) {
  X <- Corr.Object$X
  ind.x <- Corr.Object$ind.cov
  if (is.null(Omega.ols)) {
    Omega.ols <- cbind(Corr.Object$Omega)
  }
  Omega.ols <- cbind(Omega.ols)
  X <- cbind(X)
  d <- length(ind.x)
  if (d != ncol(Omega.ols)) {
    warning("Make sure Omega.ols is K x d, d = # of covariates of interest")
    return(0)
  }
  
  inv.SigmaX <- solve(t(X) %*% X)
  proxy.n <- rep(NA, d)
  for (i in 1:d) {
    proxy.n[i] <- 1/inv.SigmaX[ind.x[i], ind.x[i]]
  }
  sum.Omega.ols <- diag(t(Omega.ols) %*% Omega.ols)
  
  p.values <- pchisq( sum.Omega.ols * proxy.n, df=nrow(Omega.ols), lower.tail = F )
  return(list(p.values=p.values, proxy.n = proxy.n))
}

