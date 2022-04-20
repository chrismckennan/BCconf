#' Estimate the effects of interest in omic data while accounting for unobserved covariates
#' 
#' This function estimates the regression coefficients of interest while accounting for unobserved covariates (i.e. latent factors/confounding factors), like cell composition and technical effects. It assumes the number of latent factors in known.
#' 
#' @param Y a p x n data matrix, where p = #of units (i.e. genes) and n = #of samples
#' @param X an n x d model matrix, where d = total number of observed covariates to include in the model. This includes the covariates of interest (i.e. disease status), as well as nuisance covariates like the intercept, plate number, DNA concentration, etc.
#' @param ind.cov a vector indicating the columns of X that are of interest. For example, if the phenotype of interest is in the second column, ind.cov = c(2)
#' @param r.confound The number of latent factors to include in the model. This can be estimated using \code{est.conf.num} and must be >= 0.
#' @param max.iter.svd The number of SVD iterations. When p is large (i.e. DNA methlation or gene expression), the default of 1 will suffice. For smaller p (i.e. metabolomic data), a value of 3 is better.
#' @param method The confounder correction method to use. Both "ols" and "theory" use the bias-corrected method described in \url{https://arxiv.org/abs/1801.00865v2}, where "ols" weights units by their inverse variance. "iterative" iteratively estimates latent confounds by removing genomic units associated with the covariate of interest from the estimation procedure. If none of these methods is used, we recommend using "huber".
#' @param shrink.Sigma If \code{TRUE}, residual variances are shrunk using empirical Bayes with an inverse gamma prior. The default of \code{FALSE} is recommended.
#' @param min.true.eig Minimum eigenvalue of (n-d)/p * L'L to consider, where L are the p x r.confound confounding effects. Default is \code{0.1}.
#' @param svd.method SVD method. If "fast", \code{irlba} is used. The default is "slow".
#' 
#' @return A list. \item{B}{A \code{p} x \code{length(ind.cov)} matrix containing the estimated confounder-adjusted effects of interest for each unit (i.e. gene).} \item{Sigma}{A p-vector containing the estimated residual variances} \item{p.values}{A \code{p} x \code{length(ind.cov)} matrix containing the confounder adjusted P values of interest} \item{p.values.ols}{A \code{p} x \code{length(ind.cov)} matrix containing the unadjusted P values. These SHOULD NOT BE USED EXCEPT FOR COMPARISON} \item{p.values.f}{A length p-vector containing the confounder adjusted F-test p-values. These are only reported if \code{length(ind.cov)} > 1.} \item{z.scores}{A \code{p} x \code{length(ind.cov)} matrix of confounder adjusted Z-scores. These should be N(0,1) under the null hypothesis} \item{z.scores.ols}{A \code{p} x \code{length(ind.cov)} matrix of unadjusted Z-scores. These SHOULD NOT BE USED EXCEPT FOR COMPARISON} \item{C}{The estimated \code{n} x \code{r.confound} matrix of unobserved covariates.} \item{L}{A \code{p} x \code{r.confound} matrix of confounding effects} \item{Omega}{A \code{r.confound} x \code{length(ind.cov)} matrix that estimates Omega.ols = (X[,ind.cov]' X[,ind.cov])^\{-1\}X[,ind.cov]' C} \item{all.Omegas}{A list containing Omega.ols estimates for all methods} \item{Cov.total}{An \code{n} x \code{\{ncol(X) + r.confound\}} design matrix [X C] that can be used in downstream analysis} \item{method}{The confounder correction method used} \item{ind.cov}{The columns of X that are of interest}
#' @export
Correction <- function(Y, X, ind.cov, r.confound, max.iter.svd=1, method=c("ols", "theory", "huber", "iterative", "naive", "none"), shrink.Sigma=F, min.true.eig=0.1, svd.method = "slow", iterative.params=list(n.iter=5,qvalue.thresh=0.1)) {
	p <- nrow(Y)
	n <- ncol(Y)
	d <- ncol(X)
  
  method <- match.arg(method, c("ols", "theory", "huber", "iterative", "naive", "none"))  #Default is ols
	
	svd.out <- svd.LandC(Y=Y, X=X, r=r.confound, max.iter.svd=max.iter.svd, min.true.eig=min.true.eig, svd.method = svd.method)
	L.svd <- svd.out$L   #Counfounding Effects
	Sigma.svd <- svd.out$Sigma   #Sigma.hat
	if (r.confound > 0) {
	  Sigma.ols <- Sigma.svd * (n-d-r.confound) / (n-d) + rowSums(L.svd * L.svd)
	}
	
	if (shrink.Sigma) {
		out.shrink <- my.eb(Sigma=Sigma.svd, n=n - d - r.confound)
		Sigma.shrink <- ( (n - d - r.confound)*Sigma.svd + out.shrink$beta ) / ( n - d - r.confound + out.shrink$alpha - 1 )   #Shrunk using an inverse gamma prior for the variance. I will use this to do estimation
	} else {
		Sigma.shrink <- Sigma.svd
	}
	rm(Sigma.svd)

	Y1 <- (Y %*% X %*% solve(t(X) %*% X))[,ind.cov]
  if (method != "none" && r.confound > 0) {
    Omegas <- HubersLoss(Y=Y, X=X, ind.cov=ind.cov, Sigma=Sigma.shrink, Gamma=L.svd)   #Omega.ols using Huber's loss and bias-corrected estimate. Note the Huber's loss option is NOT bias corrected and may give spurious results when you have small confounders. I would compare analysis with the two.
    if ("qvalue" %in% installed.packages()) {
      Omegas$Omega.iterative <- IterativeOmega(Y1 = Y1, L = L.svd, Sigma = Sigma.shrink, v.e = solve(t(X) %*% X)[ind.cov,ind.cov], v.L = 1/(n-d), n.iter = iterative.params$n.iter, qvalue.thresh = iterative.params$qvalue.thresh)
      qvalue.exists <- T
    } else {
      if (method == "iterative") {
        warning("qvalue package is not installed. Switching method to `ols'")
        method <- "ols"
      }
    }
  }
	if (method == "iterative") {
	  Omega.ols <- Omegas$Omega.iterative
	}
	if (method == "theory") {
	  Omega.ols <- Omegas$Omega.theory
	}
	if (method == "huber") {
	  Omega.ols <- Omegas$Omega.hub
	}
	if (method == "naive") {
	  Omega.ols <- Omegas$Omega.naive
	}
	if (method == "ols") {
	  Omega.ols <- Omegas$Omega.ols
	}
	if (method == "none") {
	  Omegas <- NULL
	  Omega.ols <- 0
	}
	
	if (method == "none") {
		B <- Y1
	} else {
		B <- Y1 - L.svd %*% Omega.ols 
	}
	
	if (length(ind.cov) == 1) {
	  t <- B / sqrt(Sigma.shrink) / sqrt(solve(t(X) %*% X)[ind.cov,ind.cov] + 1/(n - d)*sum(Omega.ols * Omega.ols))
	  if (r.confound > 0) {
	    t.ols <- Y1 / sqrt(Sigma.ols) / sqrt(solve(t(X) %*% X)[ind.cov,ind.cov])
	  } else {
	    t.ols <- t
	  }
	} else {
		var.mat <- solve(t(X) %*% X)[ind.cov,ind.cov] + 1/(n - d)*t(Omega.ols) %*% Omega.ols
		t <- (B / sqrt(Sigma.shrink)) %*% diag(1/sqrt(diag(var.mat)))
		t.ols <- (B / sqrt(Sigma.ols)) %*% diag(1/sqrt(diag(solve(t(X) %*% X)[ind.cov,ind.cov])))
	}
	
	p.values <- 2 * pt( -abs(t), df=n-d-r.confound )
	p.values.ols <- 2 * pt( -abs(t.ols), df=n-d )
	if (method != "none") {
	  C <- svd.out$C+cbind(X[,ind.cov])%*%t(cbind(Omega.ols))
	} else {
	  C <- NULL
	}
	  if (length(ind.cov) > 1) {
	    X.0 <- cbind(X[,-ind.cov], C)
	    X.A <- cbind(X, C)
	    B.0 <- Y %*% X.0 %*% solve(t(X.0)%*%X.0)
	    B.A <- Y %*% X.A %*% solve(t(X.A)%*%X.A)
      
	    num.f <- rowSums( B.A * (Y %*% X.A) ) / length(ind.cov) - rowSums( B.0 * (Y %*% X.0) ) / length(ind.cov)
	    f.stats <- num.f / Sigma.shrink
	    pvalues.f <- pf(f.stats, df1=length(ind.cov), df2=n-d-r.confound, lower.tail=F)
	    rm(B.0)
	    rm(B.A)
	  } else {
	    pvalues.f <- NULL
	  }
	out <- list(B=B, Sigma=Sigma.shrink, p.values=p.values, p.values.ols=p.values.ols, p.values.f=pvalues.f, z.scores=qnorm( pt(t,df=n-d-r.confound) ), z.scores.ols=qnorm( pt(t.ols,df=n-d) ), C=C, L=L.svd, Omega=Omega.ols, all.Omegas=Omegas, Cov.total=cbind(X,C), method=method, X=X, ind.cov=ind.cov)
	out$Omega.pvalue <- Sig.Omega(out)
	return(out)
}


#' Change the confounding correction method
#' 
#' This function quickly re-estimates the confounding adjusted regression coefficients and C using a new user-specified method
#' 
#' @param Corr.Object Output from the \code{Correction} function.
#' @param new.method New method to use. Can be one of c("ols", "theory", "huber", "naive", "none") and the default is "ols".
#' 
#' @return A \code{Correction}-type list
#' @export
ChangeMethod <- function(Corr.Object, new.method=c("ols", "theory", "huber", "iterative", "naive", "none")) {
  X <- Corr.Object$X
  ind.cov <- Corr.Object$ind.cov
  Omega.old <- Corr.Object$Omega
  d <- ncol(X)
  n <- nrow(X)
  K <- nrow(cbind(Omega.old))
  
  new.method <- match.arg(new.method, c("ols", "theory", "huber", "iterative", "naive", "none"))  #Default is ols
  if (new.method == Corr.Object$method) {
    warning("New correction method is the same as the old correction method")
    return(Corr.Object)
  }
  
  if (new.method == "theory") {
    Omega.new <- Corr.Object$all.Omegas$Omega.theory
  } 
  if (new.method == "huber") {
    Omega.new <- Corr.Object$all.Omegas$Omega.hub
  }
  if (new.method == "naive") {
    Omega.new <- Corr.Object$all.Omegas$Omega.naive
  }
  if (new.method == "ols") {
    Omega.new <- Corr.Object$all.Omegas$Omega.ols
  }
  if (new.method == "iterative") {
    Omega.new <- Corr.Object$all.Omegas$Omega.iterative
  }
  
  B.old <- Corr.Object$B
  L <- Corr.Object$L
  Sigma.shrink <- Corr.Object$Sigma
  
  B <- B.old + L %*% (cbind(Omega.old) - cbind(Omega.new))
  
  if (length(ind.cov) == 1) {
    t <- B / sqrt(Sigma.shrink) / sqrt(solve(t(X) %*% X)[ind.cov,ind.cov] + 1/(n - d)*sum(Omega.new * Omega.new))
  } else {
    var.mat <- solve(t(X) %*% X)[ind.cov,ind.cov] + 1/(n - d)*t(Omega.new) %*% Omega.new
    t <- (B / sqrt(Sigma.shrink)) %*% diag(1/sqrt(diag(var.mat)))
  }
  
  p.values <- 2 * pt( -abs(t), df=n-d-K )
  C.new <- Corr.Object$C + cbind(X[,ind.cov])%*%t(cbind(Omega.new) - cbind(Omega.old))
  
  ##Remake out##
  out <- Corr.Object
  out$method <- new.method
  out$B <- B
  out$p.values <- p.values
  out$C <- C.new
  out$Cov.total <- cbind(X, C.new)
  out$Omega <- Omega.new
  out$z.scores <- qnorm( pt(t,df=n-d-K) )
  if (length(ind.cov) > 1) {
    X.0 <- out$Cov.total[,-ind.cov]
    X.use <- X[,ind.cov]
    f.stats <- rowSums( (B%*%( t(X.use)%*%X.use - t(X.use)%*%X.0%*%solve(t(X.0)%*%X.0,t(X.0)%*%X.use) )) * B ) / length(ind.cov) / Sigma.shrink
    out$p.values.f <- pf(f.stats, df1=length(ind.cov), df2=n-d-K, lower.tail=F)
  }
  
  return(out)
}

IterativeOmega <- function(Y1, L, Sigma, v.L, v.e, n.iter=5, qvalue.thresh=0.1) {
  L <- cbind(L)
  d <- NCOL(Y1)
  V <- diag(v.L, nrow = NCOL(L), ncol = NCOL(L))
  Omega <- solve(t(L/Sigma)%*%L - V,t(L/Sigma)%*%Y1)
  if (d==1) {
    Omega <- c(Omega); Y1 <- c(Y1); v.e <- c(v.e)
    for (i in 1:n.iter) {
      beta.i <- Y1 - c(L%*%Omega)
      var.i <- v.L*sum(Omega^2) + v.e
      q.i <- qvalue::qvalue(2*pnorm(-abs(beta.i)/sqrt(Sigma)/sqrt(var.i)))$qvalues
      ind.remove <- which(q.i<=qvalue.thresh)
      if (length(ind.remove)==0) {return(Omega)}
      Omega <- solve(t(cbind(L[-ind.remove,])/Sigma[-ind.remove])%*%cbind(L[-ind.remove,]) - V,t(cbind(L[-ind.remove,])/Sigma[-ind.remove])%*%Y1[-ind.remove])
    }
    return(Omega)
  }
  for (i in 1:n.iter) {
    beta.i <- Y1 - L%*%Omega
    var.i <- v.L*t(Omega)%*%Omega + v.e
    f.i <- rowSums((beta.i%*%solve(var.i))*beta.i)/Sigma
    q.i <- qvalue::qvalue(pchisq(q = f.i, df = d, lower.tail = F))$qvalues
    ind.remove <- which(q.i<=qvalue.thresh)
    if (length(ind.remove)==0) {return(Omega)}
    Omega <- solve(t(cbind(L[-ind.remove,])/Sigma[-ind.remove])%*%cbind(L[-ind.remove,]) - V,t(cbind(L[-ind.remove,])/Sigma[-ind.remove])%*%Y1[-ind.remove,])
  }
  return(Omega)
}
