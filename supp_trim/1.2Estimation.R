# ==============================================================================
# This script implements the estimators studied in this setting. It defines:
#   - TRIM(): implements the trimming estimator.
#   - oracle(): implements the CCA / oracle estimator.
#   - MI(): implements the multiple imputation estimator.
# ==============================================================================





# ------------------------------------------------------------------------------
# TRIM(): Implements the trimming estimator.
#
# Inputs:
#   - data.input: a list containing the following elements
#       * X.obs : observed covariates
#       * M     : mediator
#       * Y     : outcome
#       * Z     : shadow variables
#       * TT    : treatment indicator
#       * R1    : missing indicator for X1
#       * R2    : missing indicator for X2
#   - thresholds: thresholds for trimming;
#
# Output:
#   - A list containing all trimming estimates.
# ------------------------------------------------------------------------------
TRIM <- function(data.input,thresholds=c(0.1,0.05),...){

  R1 <- data.input$R1
  R2 <- data.input$R2
  Rb <- R1*R2
  X <- data.input$X.obs
  Y <- data.input$Y
  M <- data.input$M
  Z <- data.input$Z
  TT <- data.input$TT
  n=length(TT)
  
  Y_norm = normalize(Y)
  M_norm = normalize(M)
  X_norm = X
  rho <- function(x){1+exp(-x)}
  
  # estimate delta_r for r = (1,0), (0,1), (0,0)
  Bas.int.mat.all = Rb*polynomial_basis3(M_norm,Y_norm,X_norm)
  Jn.all=ncol(Bas.int.mat.all)
  
  ####### r = (1,0) ########
  Bas.ext.mat.all = R1*polynomial_basis4(M_norm,Y_norm,Z,cbind(X_norm[,c(1,3)]))
  kn.all = ncol(Bas.ext.mat.all)
  
  Jn.set <- 7:Jn.all
  d.knJn = kn.all-Jn.all
  loss.set <- NULL
  delta.hat.set <- NULL
  
  # select tuning parameters
  for (Jn in Jn.set) {
    kn = Jn + d.knJn
    Bas.int.mat = Bas.int.mat.all[,1:Jn]
    Bas.ext.mat = Bas.ext.mat.all[,1:kn]
    
    # loss function
    ff0<-function(b){
      b0=b[1:Jn]
      b1=b[(Jn+1):(2*Jn)]
      delta.mat <- TT*rho(Bas.int.mat%*%b1) + (1-TT)*rho(Bas.int.mat%*%b0)
      
      return(
        sum(  ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% (rowSums(t(Bas.ext.mat*TT))-crossprod(Bas.ext.mat*TT,Rb*delta.mat)) )^2  +
                ( tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% (rowSums(t(Bas.ext.mat*(1-TT)))-crossprod(Bas.ext.mat*(1-TT),Rb*delta.mat)) )^2   )
      )
    }
    
    c.ini = rep(0,2*Jn)
    ui = rbind(-diag(2*Jn),diag(2*Jn))
    ci = c(-rep(50,2*Jn),-rep(50,2*Jn))
    
    tsls.est <-constrOptim(c.ini, ff0, NULL, ui = ui, ci = ci)$par
    delta10.hat <- (1-TT)*rho(Bas.int.mat%*%tsls.est[1:Jn]) + TT*rho(Bas.int.mat%*%tsls.est[(Jn+1):(2*Jn)])
    delta10.hat <- as.vector(delta10.hat)
    
    loss <- sum(  ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% (rowSums(t(Bas.ext.mat*TT))-crossprod(Bas.ext.mat*TT,Rb*delta10.hat)) )^2  +
                    ( tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% (rowSums(t(Bas.ext.mat*(1-TT)))-crossprod(Bas.ext.mat*(1-TT),Rb*delta10.hat)) )^2 )*(1+2*Jn/n)
    loss.set <- c(loss.set,loss)
    delta.hat.set <- cbind(delta.hat.set,delta10.hat)
    
  } 
  sel.10=which.min(loss.set)
  delta10.hat <- delta.hat.set[,sel.10]
  ss.10 = Jn.set[sel.10]
  
  ####### r = (0,1) ########
  Bas.ext.mat.all = R2*polynomial_basis4(M_norm,Y_norm,Z,cbind(X_norm[,c(2,3)]))
  kn.all = ncol(Bas.ext.mat.all)
  
  Jn.set <- 7:Jn.all
  d.knJn = kn.all-Jn.all
  loss.set <- NULL
  delta.hat.set <- NULL
  
  # select tuning parameters
  for (Jn in Jn.set) {
    kn = Jn + d.knJn
    Bas.int.mat = Bas.int.mat.all[,1:Jn]
    Bas.ext.mat = Bas.ext.mat.all[,1:kn]
    
    # loss function
    ff0<-function(b){
      b0=b[1:Jn]
      b1=b[(Jn+1):(2*Jn)]
      delta.mat <- TT*rho(Bas.int.mat%*%b1) + (1-TT)*rho(Bas.int.mat%*%b0)
      
      return(
        sum(  ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% (rowSums(t(Bas.ext.mat*TT))-crossprod(Bas.ext.mat*TT,Rb*delta.mat)) )^2  +
                ( tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)*R2))) %*% (rowSums(t(Bas.ext.mat*(1-TT)))-crossprod(Bas.ext.mat*(1-TT),Rb*delta.mat)) )^2   )
      )
    }
    
    c.ini = rep(0,2*Jn)
    ui = rbind(-diag(2*Jn),diag(2*Jn))
    ci = c(-rep(50,2*Jn),-rep(50,2*Jn))
    
    tsls.est <-constrOptim(c.ini, ff0, NULL, ui = ui, ci = ci)$par
    delta01.hat <- (1-TT)*rho(Bas.int.mat%*%tsls.est[1:Jn]) + TT*rho(Bas.int.mat%*%tsls.est[(Jn+1):(2*Jn)])
    delta01.hat <- as.vector(delta01.hat)
    
    loss <- sum(  ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% (rowSums(t(Bas.ext.mat*TT))-crossprod(Bas.ext.mat*TT,Rb*delta01.hat)) )^2  +
                    ( tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)*R2))) %*% (rowSums(t(Bas.ext.mat*(1-TT)))-crossprod(Bas.ext.mat*(1-TT),Rb*delta01.hat)) )^2)*(1+2*Jn/n)
    loss.set <- c(loss.set,loss)
    delta.hat.set <- cbind(delta.hat.set,delta01.hat)
    
  } 
  sel.01=which.min(loss.set)
  delta01.hat <- delta.hat.set[,sel.01]
  ss.01 = Jn.set[sel.01]
  
  
  ####### r = (0,0) ########
  R0 <- (R1==R2)
  Bas.ext.mat.all = R0*polynomial_basis3(M_norm,Y_norm,cbind(Z,X_norm[,3]))
  kn.all = ncol(Bas.ext.mat.all)
  
  Jn.set <- 7:Jn.all
  d.knJn = kn.all-Jn.all
  loss.set <- NULL
  delta.hat.set <- NULL
  
  # select tuning parameters
  for (Jn in Jn.set) {
    kn = Jn + d.knJn
    Bas.int.mat = Bas.int.mat.all[,1:Jn]
    Bas.ext.mat = Bas.ext.mat.all[,1:kn]
    
    # loss function
    ff0<-function(b){
      b0=b[1:Jn]
      b1=b[(Jn+1):(2*Jn)]
      delta.mat <- TT*rho(Bas.int.mat%*%b1) + (1-TT)*rho(Bas.int.mat%*%b0)
      
      return(
        sum(  ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% (rowSums(t(Bas.ext.mat*TT))-crossprod(Bas.ext.mat*TT,Rb*delta.mat)) )^2  +
                ( tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% (rowSums(t(Bas.ext.mat*(1-TT)))-crossprod(Bas.ext.mat*(1-TT),Rb*delta.mat)) )^2   )
      )
    }
    
    c.ini = rep(0,2*Jn)
    ui = rbind(-diag(2*Jn),diag(2*Jn))
    ci = c(-rep(50,2*Jn),-rep(50,2*Jn))
    
    tsls.est <-constrOptim(c.ini, ff0, NULL, ui = ui, ci = ci)$par
    delta00.hat <- (1-TT)*rho(Bas.int.mat%*%tsls.est[1:Jn]) + TT*rho(Bas.int.mat%*%tsls.est[(Jn+1):(2*Jn)])
    delta00.hat <- as.vector(delta00.hat)
    
    loss <- sum(  ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% (rowSums(t(Bas.ext.mat*TT))-crossprod(Bas.ext.mat*TT,Rb*delta00.hat)) )^2  +
                    ( tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)*R2))) %*% (rowSums(t(Bas.ext.mat*(1-TT)))-crossprod(Bas.ext.mat*(1-TT),Rb*delta00.hat)) )^2)*(1+2*Jn/n)
    loss.set <- c(loss.set,loss)
    delta.hat.set <- cbind(delta.hat.set,delta00.hat)
    
  } 
  sel.00=which.min(loss.set)
  delta00.hat <- delta.hat.set[,sel.00]
  ss.00 = Jn.set[sel.00]
  
  delta.hat = delta10.hat + delta01.hat + delta00.hat - 2
  
  # trimming
  est.trim <- trimming(data.input,delta.hat,thresholds)
  
  
  return(list(est=est.trim$est,n=est.trim$n))
  
  
}


# Intermediate function used in TRIM(). Given `delta.hat`, this function
# applies the trimming rule and returns the resulting estimates.
trimming <- function(data.input,delta.hat,thresholds=c(0.1,0.05)){
  
  est.trim <- NULL
  n.sel <- NULL
  for (th in thresholds) {
    sel <- which(delta.hat < 1/th)
    delta.trim = delta.hat[sel]
    
    R1 <- data.input$R1[sel]
    R2 <- data.input$R2[sel]
    Rb <- R1*R2
    X <- data.input$X.obs[sel,]
    Y <- data.input$Y[sel]
    M <- data.input$M[sel]
    TT <- data.input$TT[sel]
    n=length(TT)
    
    Y_norm = normalize(Y)
    M_norm = normalize(M)
    X_norm = X
    rho <- function(x){1+exp(-x)}
    
    Bas.MX.mat = Rb*polynomial_basis2(M_norm,X_norm)
    Bas.X.mat = Rb*polynomial_basis1(X_norm)
    
    gamma.hat = predict(lm(Y*delta.trim~Bas.MX.mat,weights=TT*Rb))/predict(lm(delta.trim~Bas.MX.mat,weights=TT*Rb))
    eta.hat = predict(lm(gamma.hat*delta.trim~Bas.X.mat,weights=(1-TT)*Rb))/predict(lm(delta.trim~Bas.X.mat,weights=(1-TT)*Rb))
    theta.hat = mean(Rb*delta.trim*eta.hat)
    
    Y11.hat = predict(lm(Y*delta.trim~Bas.X.mat,weights=TT*Rb))/predict(lm(delta.trim~Bas.X.mat,weights=TT*Rb))
    Y00.hat = predict(lm(Y*delta.trim~Bas.X.mat,weights=(1-TT)*Rb))/predict(lm(delta.trim~Bas.X.mat,weights=(1-TT)*Rb))
    alpha1.hat = mean(Rb*delta.trim*Y11.hat)
    alpha0.hat = mean(Rb*delta.trim*Y00.hat)
    
    nde.hat = theta.hat-alpha0.hat
    nie.hat = alpha1.hat-theta.hat
    ate.hat = alpha1.hat-alpha0.hat
    
    est.trim <- c(est.trim,theta.hat,nie.hat,nde.hat,ate.hat)
    n.sel <- c(n.sel,n)
  }
  return(list(est=est.trim,n=n.sel))
}





oracle <- function(data.input,method="oracle",var.cal=T){
  
  # data.input <- datagen(1000,1111)
  
  
  if (method=="oracle") {
    X <- data.input$X
    Y <- data.input$Y
    M <- data.input$M
    TT <- data.input$TT
  }else if(method=="ignore"){
    Rb <- data.input$Rb
    X <- data.input$X.obs[Rb==1,]
    Y <- data.input$Y[Rb==1]
    M <- data.input$M[Rb==1]
    TT <- data.input$TT[Rb==1]
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
    # variance estimates
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
  
  R1 <- data.input$R1
  R2 <- data.input$R2
  Rb <- R1*R2
  X <- data.input$X.obs
  Y <- data.input$Y
  M <- data.input$M
  TT <- data.input$TT
  n=length(TT)
  
  
  X[R1==0,1] <- NA
  X[R2==0,2] <- NA
  imp <- mice(cbind(TT,M,Y,X),maxit=2,m=2,seed=1,printFlag=F)
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


