# ==============================================================================
# This script implements the estimators studied in this setting. It defines:
#   - MEMC_MAR(): implements the proposed estimator.
#   - oracle(): implements the CCA / oracle estimator.
#   - MI(): implements the multiple imputation estimator.
# ==============================================================================





# ------------------------------------------------------------------------------
# MEMC_MAR(): Implements the proposed SIO estimator (and, optionally, variance estimation).
#
# Inputs:
#   - data.input: a list containing the following elements
#       * X.obs : observed covariates
#       * M     : mediator
#       * Y     : outcome
#       * Z     : shadow variables
#       * TT    : treatment indicator
#       * R     : missing indicator for X
#   - var.cal (logical): if TRUE, the function also performs variance estimation;
#       if FALSE, variance estimation is not computed.
#
# Output:
#   - A list containing all MEMC estimates and (if var.cal = TRUE) their
#     corresponding variance estimates.
# ------------------------------------------------------------------------------
MEMC_MAR <- function(data.input,var.cal=F,...){
  
  R <- data.input$R
  X <- data.input$X.obs
  Y <- data.input$Y
  M <- data.input$M
  Z <- data.input$Z
  TT <- data.input$TT
  ps <- data.input$ps.true
  n=length(TT)
  
  Y_norm = normalize(Y)
  M_norm = normalize(M)
  X_norm = X
  
  rho <- function(x){1+exp(-x)}  
  
  
  # estimation of delta - constrained SMD
  Bas.int.mat = R*polynomial_basis3(M_norm,Y_norm,X_norm[,2:3])
  Bas.ext.mat = polynomial_basis3(M_norm,Y_norm,X_norm[,2:3])
  Bas.MX.mat = R*polynomial_basis2(M_norm,X_norm)
  Bas.X.mat = R*polynomial_basis1(X_norm)
  
  Jn=dim(Bas.int.mat)[2]
  ln=dim(Bas.MX.mat)[2]
  sn=dim(Bas.X.mat)[2]
  
  
  # loss function
  ff0<-function(b){
    b0=b[1:Jn]
    b1=b[(Jn+1):(2*Jn)]
    delta.mat <- TT*rho(Bas.int.mat%*%b1) + (1-TT)*rho(Bas.int.mat%*%b0)

    return(
      sum(     ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% (rowSums(t(Bas.ext.mat*TT))-crossprod(Bas.ext.mat*TT,R*delta.mat)) )^2  +
                 ( tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% (rowSums(t(Bas.ext.mat*(1-TT)))-crossprod(Bas.ext.mat*(1-TT),R*delta.mat)) )^2   )
    )
  }
  
  c.ini = rep(0,2*Jn)
  ui = rbind(-diag(2*Jn),diag(2*Jn))
  ci = c(-rep(50,2*Jn),-rep(50,2*Jn))

  tsls.est <-constrOptim(c.ini, ff0, NULL, ui = ui, ci = ci)$par
  delta.hat <- (1-TT)*rho(Bas.int.mat%*%tsls.est[1:Jn]) + TT*rho(Bas.int.mat%*%tsls.est[(Jn+1):(2*Jn)])
  delta.hat <- as.vector(delta.hat)
  
  
  # estimation of gamma 
  gamma.hat = predict(lm(Y~Bas.MX.mat,weights=TT*R*delta.hat))
 
  
  # estimation of eta 
  eta.hat = predict(lm(gamma.hat~Bas.X.mat,weights=(1-TT)*R*delta.hat))
  theta.hat = mean(R*delta.hat*eta.hat)
  
  
  # estimation of NDE
  Y11.hat = predict(lm(Y~Bas.X.mat,weights=TT*R*delta.hat))
  Y00.hat = predict(lm(Y~Bas.X.mat,weights=(1-TT)*R*delta.hat))
  alpha1.hat = mean(R*delta.hat*Y11.hat)
  alpha0.hat = mean(R*delta.hat*Y00.hat)
  
  nde.hat = theta.hat-alpha0.hat
  nie.hat = alpha1.hat-theta.hat
  ate.hat = alpha1.hat-alpha0.hat
  
  std.y11=NA
  std.y10=NA
  std.y00=NA
  std.nie=NA
  std.nde=NA
  std.ate=NA
  
  # variance estimation 
  if (var.cal==T){
    
    # loss function
    loss_Deta<-function(b){
      D.eta <- 1+exp(Bas.X.mat%*%b)
      return(    sum( R*(predict(lm(delta.hat*(D.eta*(1-TT)-1)~0+Bas.X.mat,weights = R)) )^2)    )
    }
    
    tsls.est.Deta <- optim(rep(0,sn), loss_Deta, control=list(maxit=20000))$par
    D_eta <- 1+exp(Bas.X.mat%*%tsls.est.Deta)
    D_eta <- ifelse(D_eta==1,1.01,D_eta)
    
    loss_Dgamma<-function(b){
      D.gamma <- Bas.MX.mat%*%b
      return(    sum( R*(predict(lm(delta.hat*(D_eta*(1-TT)-D.gamma*TT)~0+Bas.MX.mat,weights = R)) )^2)    )
    }
    
    tsls.est.Dgamma <- optim(rep(0,ln), loss_Dgamma, control=list(maxit=20000))$par
    D_gamma <- Bas.MX.mat%*%tsls.est.Dgamma
    
    
    kappa.hat = eta.hat + D_eta*(1-TT)*(gamma.hat-eta.hat) + D_gamma*TT*(Y-gamma.hat)
    
    
    # estimate representor
    ff1<-function(b){
      b0=b[1:Jn]
      b1=b[(Jn+1):(2*Jn)]
      rho.mat <- TT*rho(Bas.int.mat%*%b1) + (1-TT)*rho(Bas.int.mat%*%b0)
      
      return(
        sum(     ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,R*rho.mat) )^2  +
                   ( tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),R*rho.mat) )^2
                 - 2*R*rho.mat* kappa.hat
        )
      )
    }
    
    tsls.est.rho <- constrOptim(c.ini, ff1, NULL, ui = ui, ci = ci)$par
    rho.hat <- (1-TT)*rho(Bas.int.mat%*%tsls.est.rho[1:Jn]) + TT*rho(Bas.int.mat%*%tsls.est.rho[(Jn+1):(2*Jn)])
    rho.hat <- as.vector(rho.hat)
    
    IF10 = R*delta.hat*kappa.hat - theta.hat + (1-R*delta.hat)*
      ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,R*rho.hat) +
          tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),R*rho.hat)   )
    
    std.y10 = sqrt(mean(IF10^2)/n)
    
    loss_f1x_i<-function(b){
      f1x_i <- 1+exp(Bas.X.mat%*%b)
      return(    sum( R*(predict(lm(delta.hat*(f1x_i*TT-1)~0+Bas.X.mat,weights = R)) )^2)    )
    }
    
    tsls.est.f1x_i <- optim(rep(0,sn), loss_f1x_i, control=list(maxit=20000))$par
    f1x_i <- 1+exp(Bas.X.mat%*%tsls.est.f1x_i) 
    plot(1/f1x_i[R==1],1-1/D_eta[R==1])
    
    chi1.hat <- ifelse(R==1, TT*(Y-Y11.hat)*f1x_i + Y11.hat, 0)
    chi0.hat <- ifelse(R==1, (1-TT)*(Y-Y00.hat)*D_eta + Y00.hat, 0)
    
    # estimate representor
    ff21<-function(b){
      b0=b[1:Jn]
      b1=b[(Jn+1):(2*Jn)]
      sig.mat <- TT*rho(Bas.int.mat%*%b1) + (1-TT)*rho(Bas.int.mat%*%b0)
      
      return(
        sum(     ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,R*sig.mat) )^2  +
                   ( tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),R*sig.mat) )^2
                 - 2*R*sig.mat* chi1.hat))
    }
    
    ff20<-function(b){
      b0=b[1:Jn]
      b1=b[(Jn+1):(2*Jn)]
      sig.mat <-  (1-TT)*rho(Bas.int.mat%*%b0) + TT*rho(Bas.int.mat%*%b1)
      
      return(
        sum(     ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,R*sig.mat) )^2  +
                   ( tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),R*sig.mat) )^2
                 - 2*R*sig.mat* chi0.hat))
    }
    
    tsls.est.sig1 <- constrOptim(c.ini, ff21, NULL, ui = ui, ci = ci)$par
    tsls.est.sig0 <- constrOptim(c.ini, ff20, NULL, ui = ui, ci = ci)$par
    
    sig1.hat <- (1-TT)*rho(Bas.int.mat%*%tsls.est.sig1[1:Jn]) + TT*rho(Bas.int.mat%*%tsls.est.sig1[(Jn+1):(2*Jn)])
    sig1.hat <- as.vector(sig1.hat)
    
    sig0.hat <- (1-TT)*rho(Bas.int.mat%*%tsls.est.sig0[1:Jn]) + TT*rho(Bas.int.mat%*%tsls.est.sig0[(Jn+1):(2*Jn)])
    sig0.hat <- as.vector(sig0.hat)
    
    
    IF11 = R*delta.hat*chi1.hat - alpha1.hat + (1-R*delta.hat)*
      ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,R*sig1.hat) +
          tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),R*sig1.hat)   )
    
    IF00 = R*delta.hat*chi0.hat - alpha0.hat + (1-R*delta.hat)*
      ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,R*sig0.hat) +
          tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),R*sig0.hat)   )
    
    std.y11 = sqrt(mean((IF11)^2)/n)
    std.y00 = sqrt(mean((IF00)^2)/n)
    
    std.nie = sqrt(mean((IF11-IF10)^2)/n)
    std.nde = sqrt(mean((IF10-IF00)^2)/n)
    std.ate = sqrt(mean((IF11-IF00)^2)/n)
  }
  
  
  
  return(list(theta.hat=theta.hat,alpha1.hat=alpha1.hat,alpha0.hat=alpha0.hat,
              nde.hat=nde.hat,nie.hat=nie.hat,ate.hat=ate.hat,
              std.y11=std.y11,std.y10=std.y10,std.y00=std.y00,
              std.nie=std.nie,std.nde=std.nde,std.ate=std.ate
  ))
  
  
}


# ------------------------------------------------------------------------------
# oracle(): Implements the oracle / CCA (complete-case analysis) estimator.
#
# Inputs:
#   - data.input: a list containing the following elements
#       * X     : underlying covariates
#       * X.obs : observed covariates
#       * M     : mediator
#       * Y     : outcome
#       * TT    : treatment indicator
#       * Rb    : missing indicator for covariates
#   - method (character): estimation method to use
#       * "oracle": oracle estimator
#       * "ignore": CCA estimator
#   - var.cal (logical): if TRUE, the function also performs variance estimation;
#       if FALSE, variance estimation is not computed.
#
# Output:
#   - A list containing the oracle / CCA estimates and (if var.cal = TRUE) their
#     corresponding variance estimates.
# ------------------------------------------------------------------------------
oracle <- function(data.input,method="oracle",var.cal=T){  
  
  if (method=="oracle") {
    X <- data.input$X
    Y <- data.input$Y
    M <- data.input$M
    TT <- data.input$TT
  }else if(method=="ignore"){
    R <- data.input$R
    X <- data.input$X.obs[R==1,]
    Y <- data.input$Y[R==1]
    M <- data.input$M[R==1]
    TT <- data.input$TT[R==1]
  }
  n=length(TT)
  M_norm = normalize(M)
  
  Bas.MX.mat = polynomial_basis2(M_norm,X)
  gamma.hat = predict(lm(Y~0+Bas.MX.mat,weights = TT))
  Bas.X.mat = polynomial_basis1(X)
  eta.hat = predict(lm(gamma.hat~0+Bas.X.mat,weights = 1-TT))
  theta.hat = mean(eta.hat)
  
  
  # estimation of nde
  Bas.X.mat = polynomial_basis1(X)
  Y11.hat = predict(lm(Y~0+Bas.X.mat,weights=TT))
  Y00.hat = predict(lm(Y~0+Bas.X.mat,weights=1-TT))
  
  alpha1.hat = mean(Y11.hat) 
  alpha0.hat = mean(Y00.hat) 
  
  nde.hat = theta.hat-alpha0.hat
  nie.hat = alpha1.hat-theta.hat
  ate.hat = alpha1.hat-alpha0.hat
  
  std.y11=NA
  std.y10=NA
  std.y00=NA
  std.nie=NA
  std.nde=NA
  std.ate=NA
  
  if (var.cal==T){
  
  # variance estimation
  f0x = tcrossprod(Bas.X.mat, ginv(crossprod(Bas.X.mat,Bas.X.mat)) ) %*% crossprod(Bas.X.mat,1-TT)
  f0x <- ifelse(f0x<0.01,0.01,f0x)
  f0x <- ifelse(f0x>0.99,0.99,f0x)
  D_eta <- 1/f0x
  
  f0mx <- tcrossprod(Bas.MX.mat, ginv(crossprod(Bas.MX.mat,Bas.MX.mat)) ) %*% crossprod(Bas.MX.mat,1-TT)
  f0mx <- ifelse(f0mx<0.01,0.01,f0mx)
  f0mx <- ifelse(f0mx>0.99,0.99,f0mx)
  D_gamma = f0mx/(1-f0mx)/f0x
  
  
  IF10 = eta.hat + D_eta*(1-TT)*(gamma.hat-eta.hat) + D_gamma*TT*(Y-gamma.hat) - theta.hat
  std.y10 = sqrt(mean(IF10^2)/n)
  
  # IF of Y11/Y00
  IF11 <- TT*(Y-Y11.hat)/(1-f0x) + Y11.hat - alpha1.hat
  IF00 <- (1-TT)*(Y-Y00.hat)/f0x + Y00.hat - alpha0.hat
  
  std.y11 = sqrt(mean((IF11)^2)/n)
  std.y00 = sqrt(mean((IF00)^2)/n)
  
  std.nie = sqrt(mean((IF11-IF10)^2)/n)
  std.nde = sqrt(mean((IF10-IF00)^2)/n)
  std.ate = sqrt(mean((IF11-IF00)^2)/n)
  
  }
  
  return(list(theta.hat=theta.hat,alpha1.hat=alpha1.hat,alpha0.hat=alpha0.hat,
              nde.hat=nde.hat,nie.hat=nie.hat,ate.hat=ate.hat,
              std.y11=std.y11,std.y10=std.y10,std.y00=std.y00,
              std.nie=std.nie,std.nde=std.nde,std.ate=std.ate
  ))
  
}


 


# ------------------------------------------------------------------------------
# MI(): Implements the multiple imputation (MI) estimator.
#
# Inputs:
#   - data.input: a list containing the following elements
#       * X.obs : observed covariates
#       * M     : mediator
#       * Y     : outcome
#       * TT    : treatment indicator
#       * R1    : missing indicator for X1
#       * R2    : missing indicator for X2
#   - bs.rep (integer): number of bootstrap resampling replicates used for
#       inference (e.g., variance / standard error estimation).
#
# Output:
#   - A list containing the MI estimates and their bootstrap-based variance
#     estimates.
# ------------------------------------------------------------------------------
MI <- function(data.input,bs.rep=50){

  R <- data.input$R
  X <- data.input$X.obs
  X.true <- data.input$X.true[,1]
  Y <- data.input$Y
  M <- data.input$M
  TT <- data.input$TT
  n=length(TT)
  
  X[R==0,1] <- NA
  imp <- mice(cbind(TT,M,Y,X),method="rf")
  data.imp <- data.input
  data.imp$X <- as.matrix(complete(imp)[,-1:-3])
  point.est <- oracle(data.imp,var.cal=F)
  
  std.y11=NA
  std.y10=NA
  std.y00=NA
  std.nie=NA
  std.nde=NA
  std.ate=NA
  
  
  # bootstrap
  if(bs.rep>0){
    theta.list <- rep(NA,bs.rep)
    alpha1.list <- rep(NA,bs.rep)
    alpha0.list <- rep(NA,bs.rep)
    for (trial in 1:bs.rep) {
      bs.sample <- sample(n,n,replace=T)
      A <- data.imp$TT[bs.sample]
      M <- data.imp$M[bs.sample]
      Y <- data.imp$Y[bs.sample]
      X <- data.imp$X[bs.sample,]
      
      data.bs <- list(TT=A,M=M,Y=Y,X=X)
      est <- oracle(data.bs,var.cal=F)
      theta.list[trial]=est$theta.hat
      alpha1.list[trial]=est$alpha1.hat
      alpha0.list[trial]=est$alpha0.hat
    }
    
    std.y10=sd(theta.list)
    std.y11=sd(alpha1.list)
    std.y00=sd(alpha0.list)
    std.nie=sd(alpha1.list-theta.list)
    std.nde=sd(theta.list-alpha0.list)
    std.ate=sd(alpha1.list-alpha0.list)
    }
  
  return(list(theta.hat=point.est$theta.hat,alpha1.hat=point.est$alpha1.hat,alpha0.hat=point.est$alpha0.hat,
              nde.hat=point.est$nde.hat,nie.hat=point.est$nie.hat,ate.hat=point.est$ate.hat,
              std.y11=std.y11,std.y10=std.y10,std.y00=std.y00,
              std.nie=std.nie,std.nde=std.nde,std.ate=std.ate))
  
}


