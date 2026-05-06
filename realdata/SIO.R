rm(list=ls())
source("functions.R")
load("mediation.rdata")
library(MASS)
library(dplyr)


# Current time
current_time <- Sys.time()
print(current_time)


# data preparation
R <- ifelse(is.na(mediation$income),0,1)
mediation$marriage.status <- ifelse(mediation$marriage==1,1,0)
X.com <- mediation[,c("familysize","health","gender","house","age","marriage.status")]
X.com[,c(1,2,4,5)] <- apply(X.com[,c(1,2,4,5)],2,normalize)
X.mis <- mediation$income
X.mis <- normalize(X.mis)
X.mis <- ifelse(is.na(X.mis), -1, X.mis)
Y <- mediation$depression
Y <- Y-mean(Y)
Y_norm <- normalize(Y)
M <- mediation$wellbeing
M <- normalize(M)
Z <- mediation$engel
Z <- normalize(Z)
TT <- ifelse(mediation$jobsatisfication>3,1,0)
n=length(TT)
rho <- function(x){1+exp(-x)}



# ==========================================
# Use engel as shadow variables
# ==========================================
Bas.int.mat = R*polynomial_basis1(X.com,X.mis)
Bas.ext.mat.sio1 = polynomial_basis1(X.com,Z)

Jn=ncol(Bas.int.mat)
kn = ncol(Bas.ext.mat.sio1)

# loss function
ff0<-function(b){
  b0=b[1:Jn]
  b1=b[(Jn+1):(2*Jn)]
  delta.mat <- TT*rho(Bas.int.mat%*%b1) + (1-TT)*rho(Bas.int.mat%*%b0)

  return(
    sum(     ( tcrossprod( TT*Bas.ext.mat.sio1,ginv(crossprod(Bas.ext.mat.sio1*TT))) %*% (rowSums(t(Bas.ext.mat.sio1*TT))-crossprod(Bas.ext.mat.sio1*TT,R*delta.mat)) )^2  +
               ( tcrossprod( (1-TT)*Bas.ext.mat.sio1,ginv(crossprod(Bas.ext.mat.sio1*(1-TT)))) %*% (rowSums(t(Bas.ext.mat.sio1*(1-TT)))-crossprod(Bas.ext.mat.sio1*(1-TT),R*delta.mat)) )^2   )
  )
}

c.ini = rep(0,2*Jn)
ui = rbind(-diag(2*Jn),diag(2*Jn))
ci = c(-rep(50,2*Jn),-rep(50,2*Jn))
print(ci)

tsls.est <-constrOptim(c.ini, ff0, NULL, ui = ui, ci = ci)$par
delta.hat <- (1-TT)*rho(Bas.int.mat%*%tsls.est[1:Jn]) + TT*rho(Bas.int.mat%*%tsls.est[(Jn+1):(2*Jn)])
delta.hat <- as.vector(delta.hat)

# estimation of Y10
Bas.MX.mat = R*polynomial_basis2(M,X.com,X.mis)
Bas.X.mat = R*polynomial_basis1(X.com,X.mis)

gamma.hat = predict(lm(Y~0+Bas.MX.mat,weights=TT*R*delta.hat))
eta.hat = predict(lm(gamma.hat~0+Bas.X.mat,weights=(1-TT)*R*delta.hat))
theta.hat.sio1 = mean(R*delta.hat*eta.hat)

# estimation of Y11/Y00
Y11.hat.sio1 = predict(lm(Y~Bas.X.mat,weights=TT*R*delta.hat))
Y00.hat.sio1 = predict(lm(Y~Bas.X.mat,weights=(1-TT)*R*delta.hat))

alpha1.hat.sio1 = mean(R*delta.hat*Y11.hat.sio1)
alpha0.hat.sio1 = mean(R*delta.hat*Y00.hat.sio1)


# estimation of NIE/NDE
nie.hat.sio1 =  alpha1.hat.sio1 - theta.hat.sio1
print(paste0("NIE.sio1: ",nie.hat.sio1))

nde.hat.sio1 = theta.hat.sio1 - alpha0.hat.sio1
print(paste0("NDE.sio1: ",nde.hat.sio1))

ate.hat.sio1 = alpha1.hat.sio1 - alpha0.hat.sio1
print(paste0("ATE.sio1: ",ate.hat.sio1))




# ========================================================
# MNAR - M,Y seen as shadow variable
# ========================================================
Bas.ext.mat.sio2 = polynomial_basis2(M,X.com,Y_norm)
kn = ncol(Bas.ext.mat.sio2)

ff0<-function(b){
  b0=b[1:Jn]
  b1=b[(Jn+1):(2*Jn)]
  delta.mat <- TT*rho(Bas.int.mat%*%b1) + (1-TT)*rho(Bas.int.mat%*%b0)

  return(
    sum(     ( tcrossprod( TT*Bas.ext.mat.sio2,ginv(crossprod(Bas.ext.mat.sio2*TT))) %*% (rowSums(t(Bas.ext.mat.sio2*TT))-crossprod(Bas.ext.mat.sio2*TT,R*delta.mat)) )^2  +
               ( tcrossprod( (1-TT)*Bas.ext.mat.sio2,ginv(crossprod(Bas.ext.mat.sio2*(1-TT)))) %*% (rowSums(t(Bas.ext.mat.sio2*(1-TT)))-crossprod(Bas.ext.mat.sio2*(1-TT),R*delta.mat)) )^2   )
  )
}

tsls.est <-constrOptim(c.ini, ff0, NULL, ui = ui, ci = ci)$par
delta.hat <- (1-TT)*rho(Bas.int.mat%*%tsls.est[1:Jn]) + TT*rho(Bas.int.mat%*%tsls.est[(Jn+1):(2*Jn)])
delta.hat <- as.vector(delta.hat)


# estimation of Y10
gamma.hat = tcrossprod(Bas.MX.mat, ginv(crossprod(TT*R*delta.hat*Bas.MX.mat,Bas.MX.mat)) ) %*% crossprod(Bas.MX.mat*TT*R*delta.hat,Y)
eta.hat = tcrossprod(Bas.X.mat, ginv(crossprod((1-TT)*R*delta.hat*Bas.X.mat,Bas.X.mat)) ) %*% crossprod(Bas.X.mat*(1-TT)*R*delta.hat,gamma.hat)
theta.hat.sio2 = mean(R*delta.hat*eta.hat)

Y11.hat = tcrossprod(Bas.X.mat, ginv(crossprod(TT*R*delta.hat*Bas.X.mat,Bas.X.mat)) ) %*% crossprod(Bas.X.mat*TT*R*delta.hat,Y)
Y00.hat = tcrossprod(Bas.X.mat, ginv(crossprod((1-TT)*R*delta.hat*Bas.X.mat,Bas.X.mat)) ) %*% crossprod(Bas.X.mat*(1-TT)*R*delta.hat,Y)

alpha1.hat.sio2 = mean(R*delta.hat*Y11.hat)
alpha0.hat.sio2 = mean(R*delta.hat*Y00.hat)

# estimation of NIE/NDE
nie.hat.sio2 =  alpha1.hat.sio2 - theta.hat.sio2
print(paste0("NIE.sio2: ",nie.hat.sio2))

nde.hat.sio2 = theta.hat.sio2 - alpha0.hat.sio2
print(paste0("NDE.sio2: ",nde.hat.sio2))

ate.hat.sio2 = alpha1.hat.sio2 - alpha0.hat.sio2
print(paste0("ATE.sio2: ",ate.hat.sio2))

print(Sys.time())
print(Sys.time()-current_time)

# save data
save(nie.hat.sio1,nie.hat.sio2,nde.hat.sio1,nde.hat.sio2,ate.hat.sio1,ate.hat.sio2,file="SIO.rdata")


