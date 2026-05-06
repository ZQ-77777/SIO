# ==============================================================================
# This script implements the estimators studied in this setting. It defines:
#   - MEMC_DGP1(): implements the proposed estimator.
#   - oracle(): implements the CCA / oracle estimator.
#   - MI(): implements the multiple imputation estimator.
# ==============================================================================





# ------------------------------------------------------------------------------
# MEMC_DGP1(): Implements the proposed SIO estimator (and, optionally, variance estimation).
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
#   - var.cal (logical): if TRUE, the function also performs variance estimation;
#       if FALSE, variance estimation is not computed.
#
# Output:
#   - A list containing all MEMC estimates and (if var.cal = TRUE) their
#     corresponding variance estimates.
# ------------------------------------------------------------------------------
MEMC_DGP1 <- function(data.input,var.cal=F,...){

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

  
  
  # estimation of gamma 
  Bas.MX.mat.all = Rb*polynomial_basis2(M_norm,X_norm)
  Bas.X.mat.all = Rb*polynomial_basis1(X_norm)
  
  ln.all=dim(Bas.MX.mat.all)[2]
  sn.all=dim(Bas.X.mat.all)[2]
  
  ln.set <- 5:ln.all
  loss.set <- NULL
  gamma.hat.set <- NULL
  
  for (ln in ln.set) {
    Bas.MX.mat = Bas.MX.mat.all[,1:ln]
    
    gamma.lm = lm(Y~0+Bas.MX.mat,weights=TT*Rb*delta.hat)
    gamma.hat = predict(gamma.lm)
    loss = sum(weighted.residuals(gamma.lm)^2) * (1+2*ln/n)
    
    loss.set <- c(loss.set,loss)
    gamma.hat.set <- cbind(gamma.hat.set,gamma.hat)
    
  }
  sel.gamma=which.min(loss.set)
  gamma.hat <- gamma.hat.set[,sel.gamma]
  ss.gamma = ln.set[sel.gamma]

  
  
  
  
  # estimation of eta
  sn.set <- 4:sn.all
  loss.set <- NULL
  eta.hat.set <- NULL
  
  for (sn in sn.set) {
    Bas.X.mat = Bas.X.mat.all[,1:sn]
    
    eta.lm = lm(gamma.hat~0+Bas.X.mat,weights=(1-TT)*Rb*delta.hat)
    eta.hat = predict(eta.lm)
    loss = sum(weighted.residuals(eta.lm)^2) * (1+2*sn/n)
    
    loss.set <- c(loss.set,loss)
    eta.hat.set <- cbind(eta.hat.set,eta.hat)
    
  }
  sel.eta=which.min(loss.set)
  eta.hat <- eta.hat.set[,sel.eta]
  ss.eta = sn.set[sel.eta]
  
  theta.hat = mean(Rb*delta.hat*eta.hat)
  
  
  # estimation of NDE
  # Y11
  loss.set <- NULL
  Y11.hat.set <- NULL
  
  for (sn in sn.set) {
    Bas.X.mat = Bas.X.mat.all[,1:sn]
    
    Y11.lm = lm(Y~Bas.X.mat,weights=TT*Rb*delta.hat)
    Y11.hat = predict(Y11.lm)
    loss = sum(weighted.residuals(Y11.lm)^2) * (1+2*sn/n)
    
    loss.set <- c(loss.set,loss)
    Y11.hat.set <- cbind(Y11.hat.set,Y11.hat)
    
  }
  sel.Y11=which.min(loss.set)
  Y11.hat <- Y11.hat.set[,sel.Y11]
  ss.Y11 <- sn.set[sel.Y11]
  
  # Y00
  loss.set <- NULL
  Y00.hat.set <- NULL
  
  for (sn in sn.set) {
    Bas.X.mat = Bas.X.mat.all[,1:sn]
    
    Y00.lm =lm(Y~Bas.X.mat,weights=(1-TT)*Rb*delta.hat)
    Y00.hat = predict(Y00.lm)
    loss = sum(weighted.residuals(Y00.lm)^2) * (1+2*sn/n)
    
    loss.set <- c(loss.set,loss)
    Y00.hat.set <- cbind(Y00.hat.set,Y00.hat)
    
  }
  sel.Y00=which.min(loss.set)
  Y00.hat <- Y00.hat.set[,sel.Y00]
  ss.Y00 <- sn.set[sel.Y00]

  alpha1.hat = mean(Rb*delta.hat*Y11.hat)
  alpha0.hat = mean(Rb*delta.hat*Y00.hat)
  
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
    Bas.int.mat.all = Rb*polynomial_basis3(M_norm,Y_norm,X_norm)
    Bas.MX.mat = Rb*polynomial_basis2(M_norm,X_norm)
    Bas.X.mat = Rb*polynomial_basis1(X_norm)
    
    ln=dim(Bas.MX.mat)[2]
    sn=dim(Bas.X.mat)[2]
    
    # ======= y10 =================
    # loss function
    Bas.X.mat.qr = Rb*qr_basis1(X_norm)
    snqr <- ncol(Bas.X.mat.qr)
    
    loss_Deta<-function(b){
      D.eta <- 1+exp(Bas.X.mat.qr%*%b)
      return(    sum( Rb*(predict(lm(delta.hat*(D.eta*(1-TT)-1)~0+Bas.X.mat.qr,weights = Rb)) )^2)    )
    }
    
    tsls.est.Deta <- optim(rep(0,snqr), loss_Deta, control=list(maxit=20000))$par
    D_eta <- 1+exp(Bas.X.mat.qr%*%tsls.est.Deta)
    D_eta <- ifelse(D_eta==1,1.01,D_eta)
    
    loss_Dgamma<-function(b){
      D.gamma <- Bas.MX.mat%*%b
      return(    sum( Rb*(predict(lm(delta.hat*(D_eta*(1-TT)-D.gamma*TT)~0+Bas.MX.mat,weights = Rb)) )^2)    )
    }
    
    tsls.est.Dgamma <- optim(rep(0,ln), loss_Dgamma, control=list(maxit=20000))$par
    D_gamma <- Bas.MX.mat%*%tsls.est.Dgamma
    
    kappa.hat = eta.hat + D_eta*(1-TT)*(gamma.hat-eta.hat) + D_gamma*TT*(Y-gamma.hat)
    
    # ---- estimate representers
    ####### r = (1,0) ########
    Bas.ext.mat = R1*polynomial_basis4(M_norm,Y_norm,Z,cbind(X_norm[,c(1,3)]))
    
    ff1<-function(b){
      b0=b[1:Jn]
      b1=b[(Jn+1):(2*Jn)]
      rho.mat <- TT*rho(Bas.int.mat%*%b1) + (1-TT)*rho(Bas.int.mat%*%b0)
      
      return(
        sum(     ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,Rb*rho.mat) )^2  +
                   ( tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),Rb*rho.mat) )^2
                 - 2*Rb*rho.mat* kappa.hat
        )
      )
    }
    
    tsls.est.rho <- optim(c.ini, ff1, NULL, control=list(maxit=20000))$par
    rho10.hat <- (1-TT)*rho(Bas.int.mat%*%tsls.est.rho[1:Jn]) + TT*rho(Bas.int.mat%*%tsls.est.rho[(Jn+1):(2*Jn)])
    rho10.hat <- as.vector(rho10.hat)
    A.rho10   <- tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,Rb*rho10.hat) +
                 tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),Rb*rho10.hat) 

    ####### r = (0,1) ########
    Bas.ext.mat = R2*polynomial_basis4(M_norm,Y_norm,Z,cbind(X_norm[,c(2,3)]))
    
    ff1<-function(b){
      b0=b[1:Jn]
      b1=b[(Jn+1):(2*Jn)]
      rho.mat <- TT*rho(Bas.int.mat%*%b1) + (1-TT)*rho(Bas.int.mat%*%b0)
      
      return(
        sum(     ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,Rb*rho.mat) )^2  +
                   ( tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),Rb*rho.mat) )^2
                 - 2*Rb*rho.mat* kappa.hat
        )
      )
    }
    
    tsls.est.rho <- optim(c.ini, ff1, NULL, control=list(maxit=20000))$par
    rho01.hat <- (1-TT)*rho(Bas.int.mat%*%tsls.est.rho[1:Jn]) + TT*rho(Bas.int.mat%*%tsls.est.rho[(Jn+1):(2*Jn)])
    rho01.hat <- as.vector(rho01.hat)
    A.rho01   <- tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,Rb*rho01.hat) +
      tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),Rb*rho01.hat) 
    
    
    ####### r = (0,0) ########
    R0 <- (R1==R2)
    Bas.ext.mat = R0*polynomial_basis3(M_norm,Y_norm,cbind(Z,X_norm[,3]))
    
    ff1<-function(b){
      b0=b[1:Jn]
      b1=b[(Jn+1):(2*Jn)]
      rho.mat <- TT*rho(Bas.int.mat%*%b1) + (1-TT)*rho(Bas.int.mat%*%b0)
      
      return(
        sum(     ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,Rb*rho.mat) )^2  +
                   ( tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),Rb*rho.mat) )^2
                 - 2*Rb*rho.mat* kappa.hat
        )
      )
    }
    
    tsls.est.rho <- optim(c.ini, ff1, NULL, control=list(maxit=20000))$par
    rho00.hat <- (1-TT)*rho(Bas.int.mat%*%tsls.est.rho[1:Jn]) + TT*rho(Bas.int.mat%*%tsls.est.rho[(Jn+1):(2*Jn)])
    rho00.hat <- as.vector(rho10.hat)
    A.rho00   <- tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,Rb*rho00.hat) +
      tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),Rb*rho00.hat) 
    
    IF10 = Rb*delta.hat*kappa.hat - theta.hat + (1-Rb*delta10.hat)*A.rho10 + (1-Rb*delta01.hat)*A.rho01 + (1-Rb*delta00.hat)*A.rho00
    std.y10 = sqrt(mean(IF10^2)/n)
    
    #  Y11/Y00
    loss_f1x_i<-function(b){
      f1x_i <- 1+exp(Bas.X.mat.qr%*%b)
      return(    sum( Rb*(predict(lm(delta.hat*(f1x_i*TT-1)~0+Bas.X.mat.qr,weights = Rb)) )^2)    )
    }
    
    tsls.est.f1x_i <- optim(rep(0,snqr), loss_f1x_i, control=list(maxit=20000))$par
    f1x_i <- 1+exp(Bas.X.mat.qr%*%tsls.est.f1x_i)
    plot(1/f1x_i[Rb==1],1-1/D_eta[Rb==1])
    
    chi1.hat <- ifelse(Rb==1, TT*(Y-Y11.hat)*f1x_i + Y11.hat, 0)
    chi0.hat <- ifelse(Rb==1, (1-TT)*(Y-Y00.hat)*D_eta + Y00.hat, 0)
    
    # ----- estimate representers
    ####### r = (1,0) ########
    Bas.ext.mat = R1*polynomial_basis4(M_norm,Y_norm,Z,cbind(X_norm[,c(1,3)]))
    
    ff21<-function(b){
      b0=b[1:Jn]
      b1=b[(Jn+1):(2*Jn)]
      sig.mat <- TT*rho(Bas.int.mat%*%b1) + (1-TT)*rho(Bas.int.mat%*%b0)
      
      return(
        sum(     ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,Rb*sig.mat) )^2  +
                   ( tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),Rb*sig.mat) )^2
                 - 2*Rb*sig.mat* chi1.hat))
    }
    
    ff20<-function(b){
      b0=b[1:Jn]
      b1=b[(Jn+1):(2*Jn)]
      sig.mat <-  (1-TT)*rho(Bas.int.mat%*%b0) + TT*rho(Bas.int.mat%*%b1)
      
      return(
        sum(     ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,Rb*sig.mat) )^2  +
                   ( tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),Rb*sig.mat) )^2
                 - 2*Rb*sig.mat* chi0.hat))
    }
    
    
    tsls.est.sig1 <- optim(c.ini, ff21, NULL, control=list(maxit=20000))$par
    tsls.est.sig0 <- optim(c.ini, ff20, NULL, control=list(maxit=20000))$par

    sig1.10.hat <- (1-TT)*rho(Bas.int.mat%*%tsls.est.sig1[1:Jn]) + TT*rho(Bas.int.mat%*%tsls.est.sig1[(Jn+1):(2*Jn)])
    sig1.10.hat <- as.vector(sig1.10.hat)
    A.sig1.10   <- ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,Rb*sig1.10.hat) +
                  tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),Rb*sig1.10.hat)   )
    
    sig0.10.hat <- (1-TT)*rho(Bas.int.mat%*%tsls.est.sig0[1:Jn]) + TT*rho(Bas.int.mat%*%tsls.est.sig0[(Jn+1):(2*Jn)])
    sig0.10.hat <- as.vector(sig0.10.hat)
    A.sig0.10   <- ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,Rb*sig0.10.hat) +
        tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),Rb*sig0.10.hat)   )
    

    ####### r = (0,1) ########
    Bas.ext.mat = R2*polynomial_basis4(M_norm,Y_norm,Z,cbind(X_norm[,c(2,3)]))
    
    ff21<-function(b){
      b0=b[1:Jn]
      b1=b[(Jn+1):(2*Jn)]
      sig.mat <- TT*rho(Bas.int.mat%*%b1) + (1-TT)*rho(Bas.int.mat%*%b0)
      
      return(
        sum(     ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,Rb*sig.mat) )^2  +
                   ( tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),Rb*sig.mat) )^2
                 - 2*Rb*sig.mat* chi1.hat))
    }
    
    ff20<-function(b){
      b0=b[1:Jn]
      b1=b[(Jn+1):(2*Jn)]
      sig.mat <-  (1-TT)*rho(Bas.int.mat%*%b0) + TT*rho(Bas.int.mat%*%b1)
      
      return(
        sum(     ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,Rb*sig.mat) )^2  +
                   ( tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),Rb*sig.mat) )^2
                 - 2*Rb*sig.mat* chi0.hat))
    }
    
    
    tsls.est.sig1 <- optim(c.ini, ff21, NULL, control=list(maxit=20000))$par
    tsls.est.sig0 <- optim(c.ini, ff20, NULL, control=list(maxit=20000))$par
    
    sig1.01.hat <- (1-TT)*rho(Bas.int.mat%*%tsls.est.sig1[1:Jn]) + TT*rho(Bas.int.mat%*%tsls.est.sig1[(Jn+1):(2*Jn)])
    sig1.01.hat <- as.vector(sig1.01.hat)
    A.sig1.01   <- ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,Rb*sig1.01.hat) +
                    tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),Rb*sig1.01.hat)   )
    
    sig0.01.hat <- (1-TT)*rho(Bas.int.mat%*%tsls.est.sig0[1:Jn]) + TT*rho(Bas.int.mat%*%tsls.est.sig0[(Jn+1):(2*Jn)])
    sig0.01.hat <- as.vector(sig0.01.hat)
    A.sig0.01   <- ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,Rb*sig0.01.hat) +
                    tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),Rb*sig0.01.hat)   )
    
    ####### r = (0,0) ########
    Bas.ext.mat = R0*polynomial_basis3(M_norm,Y_norm,cbind(Z,X_norm[,3]))
    
    ff21<-function(b){
      b0=b[1:Jn]
      b1=b[(Jn+1):(2*Jn)]
      sig.mat <- TT*rho(Bas.int.mat%*%b1) + (1-TT)*rho(Bas.int.mat%*%b0)
      
      return(
        sum(     ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,Rb*sig.mat) )^2  +
                   ( tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),Rb*sig.mat) )^2
                 - 2*Rb*sig.mat* chi1.hat))
    }
    
    ff20<-function(b){
      b0=b[1:Jn]
      b1=b[(Jn+1):(2*Jn)]
      sig.mat <-  (1-TT)*rho(Bas.int.mat%*%b0) + TT*rho(Bas.int.mat%*%b1)
      
      return(
        sum(     ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,Rb*sig.mat) )^2  +
                   ( tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),Rb*sig.mat) )^2
                 - 2*Rb*sig.mat* chi0.hat))
    }
    
    
    tsls.est.sig1 <- optim(c.ini, ff21, NULL, control=list(maxit=20000))$par
    tsls.est.sig0 <- optim(c.ini, ff20, NULL, control=list(maxit=20000))$par
    
    sig1.00.hat <- (1-TT)*rho(Bas.int.mat%*%tsls.est.sig1[1:Jn]) + TT*rho(Bas.int.mat%*%tsls.est.sig1[(Jn+1):(2*Jn)])
    sig1.00.hat <- as.vector(sig1.00.hat)
    A.sig1.00   <- ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,Rb*sig1.00.hat) +
                    tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),Rb*sig1.00.hat)   )
    
    sig0.00.hat <- (1-TT)*rho(Bas.int.mat%*%tsls.est.sig0[1:Jn]) + TT*rho(Bas.int.mat%*%tsls.est.sig0[(Jn+1):(2*Jn)])
    sig0.00.hat <- as.vector(sig0.00.hat)
    A.sig0.00   <- ( tcrossprod( TT*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*TT))) %*% crossprod(Bas.ext.mat*TT,Rb*sig0.00.hat) +
                    tcrossprod( (1-TT)*Bas.ext.mat,ginv(crossprod(Bas.ext.mat*(1-TT)))) %*% crossprod(Bas.ext.mat*(1-TT),Rb*sig0.00.hat)   )
    
    
    
    IF11 = Rb*delta.hat*chi1.hat - alpha1.hat + (1-Rb*delta10.hat)*A.sig1.10 + (1-Rb*delta01.hat)*A.sig1.01 + (1-Rb*delta00.hat)*A.sig1.00
    IF00 = Rb*delta.hat*chi0.hat - alpha0.hat +  (1-Rb*delta10.hat)*A.sig0.10 + (1-Rb*delta01.hat)*A.sig0.01 + (1-Rb*delta00.hat)*A.sig0.00
    
    std.y11 = sqrt(mean((IF11)^2)/n)
    std.y00 = sqrt(mean((IF00)^2)/n)

    std.nie = sqrt(mean((IF11-IF10)^2)/n)
    std.nde = sqrt(mean((IF10-IF00)^2)/n)
    std.ate = sqrt(mean((IF11-IF00)^2)/n)
  }

  return(list(theta.hat=theta.hat,alpha1.hat=alpha1.hat,alpha0.hat=alpha0.hat,
              nde.hat=nde.hat,nie.hat=nie.hat,ate.hat=ate.hat,
              std.y11=std.y11,std.y10=std.y10,std.y00=std.y00,
              std.nie=std.nie,std.nde=std.nde,std.ate=std.ate,
              ss.01=ss.01,ss.10=ss.10,ss.00=ss.00,
              ss.gamma=ss.gamma,ss.eta=ss.eta,ss.Y11=ss.Y11,ss.Y00=ss.Y00
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
    
    # estimation of var 
    
    Bas.X.mat.qr = qr_basis1(X)
    snqr <- ncol(Bas.X.mat.qr)
    ln <- ncol(Bas.MX.mat)
    
    loss_Deta<-function(b){
      D.eta <- 1+exp(Bas.X.mat.qr%*%b)
      return(    sum( (predict(lm((D.eta*(1-TT)-1)~0+Bas.X.mat.qr)) )^2)    )
    }
    
    tsls.est.Deta <- optim(rep(0,snqr), loss_Deta, control=list(maxit=20000))$par
    D_eta <- 1+exp(Bas.X.mat.qr%*%tsls.est.Deta)
    D_eta <- ifelse(D_eta==1,1.01,D_eta)
    
    
    loss_Dgamma<-function(b){
      D.gamma <- Bas.MX.mat%*%b
      return(    sum( (predict(lm((D_eta*(1-TT)-D.gamma*TT)~0+Bas.MX.mat)) )^2)    )
    }
    
    tsls.est.Dgamma <- optim(rep(0,ln), loss_Dgamma, control=list(maxit=20000))$par
    D_gamma <- Bas.MX.mat%*%tsls.est.Dgamma
    
    IF10 = eta.hat + D_eta*(1-TT)*(gamma.hat-eta.hat) + D_gamma*TT*(Y-gamma.hat) - theta.hat
    std.y10 = sqrt(mean(IF10^2)/n)
    
    # IF of Y11/Y00
    loss_f1x_i<-function(b){
      f1x_i <- 1+exp(Bas.X.mat.qr%*%b)
      return(  sum( (predict(lm((f1x_i*TT-1)~0+Bas.X.mat.qr)) )^2)    )
    }
    
    tsls.est.f1x_i <- optim(rep(0,snqr), loss_f1x_i, control=list(maxit=20000))$par
    f1x_i <- 1+exp(Bas.X.mat.qr%*%tsls.est.f1x_i)
    
    
    IF11 <- TT*(Y-Y11.hat)*f1x_i + Y11.hat - alpha1.hat
    IF00 <- (1-TT)*(Y-Y00.hat)*D_eta + Y00.hat - alpha0.hat
    
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
  imp <- mice(cbind(TT,M,Y,X),method="norm.predict",maxit=2,m=2,seed=1,printFlag=F)
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


