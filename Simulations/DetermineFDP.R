###Determine FDP as a function of q-value###
require(cate)
require(qvalue)
require(minfi)
require(sva)
require(dSVA)
require(ruv)

AnalyzeData <- function(Sim, q.points, K.try, n.controls=30, my.method="theory", analyses = c("my", "cate", "sva", "oracle", "dsva", "ruv", "none"), svd.method="fast") {
  out <- list()
  n <- ncol(Sim$Y) + 1
  Q <- qr.Q(qr(cbind(rep(1,n))), complete=T)[,2:n]
  Y <- Sim$Y %*% t(Q)
  L <- Sim$L
  X <- cbind(rep(1,n), c(rep(1,n/2), rep(0,n/2)))
  C <- Sim$C
  B.true <- Sim$B
  Sigma <- Sim$Sigma
  
  ##% of variance on X explained by C##
  X.tilde <- (t(Q) %*% X)[,2]
  Omega.OLS <- as.vector(solve(t(X.tilde)%*%X.tilde,t(X.tilde)%*%C))
  out$PercVar <- sum((Comp.Proj(C)%*%X.tilde)^2)/sum(X.tilde^2)
  
  
  ##No correction##
  if ("none" %in% analyses) {
    tmp <- CalcFDP(B=B.true, q=qvalue(CalcPvalues(X=X, Y=Y)$pvalues)$qvalues, q.points=q.points)
    out$FDP.none <- tmp$FDP
    out$Power.none <- tmp$Power
  }
  
  ##My method##
  if ("my" %in% analyses) {
    out$FDP.my <- matrix(0, nrow=length(K.try), ncol=length(q.points))
    out$Power.my <- matrix(0, nrow=length(K.try), ncol=length(q.points))
    out$PercVar.my <- rep(0,length(K.try))
    for (i in 1:length(K.try)) {
      k <- K.try[i]
      out.my <- Correction(Y = Y, X = X, ind.cov = 2, r.confound = k, method = my.method, min.true.eig = 0, svd.method = svd.method, shrink.Sigma = F)
      out$PercVar.my[i] <- sum((Comp.Proj(t(Q)%*%out.my$C)%*%X.tilde)^2)/sum(X.tilde^2)
      if (i == 1) {
        #out.my.hub <- ChangeMethod(Corr.Object = out.my, X = X, ind.cov = 2, old.method = my.method, new.method = "huber")
        #tmp <- CalcFDP(B=B.true, q=qvalue(out.my.hub$p.values)$qvalues, q.points=q.points)
        #out$FDP.my.hub <- tmp$FDP
        #out$Power.my.hub <- tmp$Power
        #out$FDP.my.naive <- CalcFDP(B=B.true, q=qvalue(ChangeMethod(Corr.Object = out.my, X = X, ind.cov = 2, old.method = "theory", new.method = "naive")$p.values)$qvalues, q.points=q.points)$FDP
      }
      tmp <- CalcFDP(B=B.true, q=qvalue(out.my$p.values)$qvalues, q.points=q.points)
      out$FDP.my[i,] <- tmp$FDP
      out$Power.my[i,] <- tmp$Power
    }
  }
  
  ##CATE##
  if ("cate" %in% analyses) {
    out.cate <- cate.fit(X.primary=cbind(X[,2]), X.nuis=cbind(X[,-2]), Y=t(Y), adj.method="rr", r=K.try[1], fa.method="ml", calibrate=F)
    tmp <- CalcFDP(B=B.true, q=qvalue(CalcPvalues(X=X, Y=Y, C=out.cate$Z)$pvalues)$qvalues, q.points=q.points)
    out$FDP.cate <- tmp$FDP
    out$Power.cate <- tmp$Power
    out$PercVar.cate <- sum((Comp.Proj(t(Q)%*%out.cate$Z)%*%X.tilde)^2)/sum(X.tilde^2)
  }
  
  ##C known##
  if ("oracle" %in% analyses) {
    tmp <- CalcFDP(B=B.true, q=qvalue(CalcPvalues(X=X, Y=Y, C=Q%*%C)$pvalues)$qvalues, q.points=q.points)
    out$FDP.known <- tmp$FDP
    out$Power.known <- tmp$Power
    #X.tilde <- (Q %*% t(Q) %*% X)[,2]
    #out$Reviewer1 <- CalcFDP(B=B.true, q=qvalue(CalcPvalues(X=X, Y=Y %*% t(Q), C=(diag(n)-X.tilde%*%solve(t(X.tilde)%*%X.tilde)%*%t(X.tilde))%*%Q%*%C))$qvalues, q.points=q.points)
    #Reviewr one's idea does NOT work
  }
  
  ##SVA##
  if ("sva" %in% analyses) {
    invisible(capture.output( out.SVA.irw <- sva(dat = Y, mod = X, mod0 = cbind(X[,1]), n.sv = K.try[1], method = "irw") ))
    #out.SVA.2step <- sva(dat = Y, mod = X, mod0 = cbind(X[,1]), n.sv = K.try[1], method = "two-step")
    tmp <- CalcFDP(B=B.true, q=qvalue(CalcPvalues(X=X, Y=Y, C=out.SVA.irw$sv)$pvalues)$qvalues, q.points=q.points)
    out$FDP.sva.irw <- tmp$FDP
    out$Power.sva.irw <- tmp$Power
    out$PercVar.sva <- sum((Comp.Proj(t(Q)%*%out.SVA.irw$sv)%*%X.tilde)^2)/sum(X.tilde^2)
    #out$FDP.sva.2step <- CalcFDP(B=B.true, q=qvalue(CalcPvalues(X=X, Y=Y, C=out.SVA.2step$sv)$pvalues)$qvalues, q.points=q.points)
  }
  
  ##dSVA##
  if ("dsva" %in% analyses) {
    out.dsva <- dSVA(Y = t(Y), X = cbind(X[,2]), ncomp = K.try[1])
    tmp <- CalcFDP(B=B.true, q=qvalue(out.dsva$Pvalue)$qvalues, q.points=q.points)
    out$FDP.dsva <- tmp$FDP
    out$Power.dsva <- tmp$Power
    out$PercVar.dsva <- sum((Comp.Proj(t(Q)%*%out.dsva$Z)%*%X.tilde)^2)/sum(X.tilde^2)
  }
  
  ##RUV##
  #Use n.controls for the number of controls
  if ("ruv" %in% analyses) {
    control.genes <- (1:nrow(Y)) %in% which(abs(B.true) < 1e-8)[sample(x=1:sum(abs(B.true) < 1e-8), size=n.controls, replace = F)]
    out.ruv2 <- ruv::RUV2(Y = t(Y), X = cbind(X[,2]), ctl=control.genes, k = K.try[1])
    out.ruv4 <- ruv::RUV4(Y = t(Y), X = cbind(X[,2]), ctl=control.genes, k = K.try[1])
    tmp.2 <- CalcFDP(B=B.true, q=qvalue(out.ruv2$p)$qvalues, q.points=q.points)
    tmp.4 <- CalcFDP(B=B.true, q=qvalue(out.ruv4$p)$qvalues, q.points=q.points)
    out$FDP.ruv2 <- tmp.2$FDP
    out$Power.ruv2 <- tmp.2$Power
    out$FDP.ruv4 <- tmp.4$FDP
    out$Power.ruv4 <- tmp.4$Power
    out$PercVar.ruv2 <- sum((Comp.Proj(t(Q)%*%out.ruv2$W)%*%X.tilde)^2)/sum(X.tilde^2)
    out$PercVar.ruv4 <- sum((Comp.Proj(t(Q)%*%out.ruv4$W)%*%X.tilde)^2)/sum(X.tilde^2)
  }

  return(out)
}

CalcFDP <- function(B, q, q.points) {  #q are qvalues
  FDP <- rep(0, length(q.points))
  Power <- rep(0, length(q.points))
  ind.B <- abs(B) > 1e-8
  ind.B.null <- !ind.B
  for (i in 1:length(q.points)) {
    ind.i  <- q <= q.points[i]
    FDP[i] <- sum(ind.B.null & ind.i) / sum(ind.i)
    Power[i] <- sum(ind.B & ind.i) / sum(ind.B)
  }
  return(list(FDP=FDP, Power=Power))
}


CalcPvalues <- function(X, Y, C=NULL) {
  if (is.null(C)) {
    Cov <- cbind(X)
  } else {
    Cov <- cbind(X, C)
  }
  df <- nrow(Cov) - ncol(Cov)
  B <- Y %*% Cov %*% solve(t(Cov) %*% Cov)
  Resids <- Y - B %*% t(Cov); Sigma <- rowSums(Resids^2 / df)
  t <- B[,2] / sqrt(Sigma) / sqrt(solve(t(Cov)%*%Cov)[2,2])
  return(list(pvalues=2*pt( -abs(t), df = df ), z=qnorm(pt(t,df=df))))
}

Comp.Proj <- function(X) {
  qr.X <- qr(cbind(X))
  Q <- qr.Q(qr.X)[,1:qr.X$rank]
  Q%*%t(Q)
}

Comp.Frac.Var <- function(X, C, Z=rep(1,length(X))) {
  if (!is.null(Z)) {
    Proj.Z <- Comp.Proj(cbind(Z))
    X <- X - Proj.Z %*% X
    C <- C - Proj.Z %*% C
  }
  sum(as.vector(Comp.Proj(cbind(C))%*%X)^2)/sum(as.vector(X)^2)
}
