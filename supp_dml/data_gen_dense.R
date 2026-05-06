expit <- function(x){exp(x)/(1+exp(x))}
normalize <- function(x){(x-min(x))/(max(x)-min(x))}


datagen <- function(n,seed=NA,w=0.6,p=10,std=1,...){
  
  if(is.na(seed)){seed=sample(9999999,1)}
  set.seed(seed)
  
  eps1 <- rnorm(n, 0, std)
  eps2 <- rnorm(n, 0, std)
  eps3 <- rnorm(n, 0, std)
  eps4 <- rnorm(n, 0, std)
  eps1p <- rnorm(n, 0, std)
  eps2p <- rnorm(n, 0, std)
  
  X <- matrix(nrow=n,ncol=p)
  X[,1] <- pnorm(eps1,sd=std)
  X[,2] <- pnorm(eps2,sd=std)
  X[,-c(1,2)] <- sapply(3:p, function(j) if(j%%2) runif(n) else rbinom(n,1,0.5) ) 
  
  Z <- cbind(pnorm(w*eps1 + sqrt(1-w^2)*eps1p,sd=std),pnorm(w*eps2 + sqrt(1-w^2)*eps2p,sd=std))
  
  ps.true <- expit(0.5* X[,1]^2 + 0.5 * X[,2] - 0.5 * X[,3]^2 + 0.5* X[,4]) # propensity score of treatment
  summary(ps.true)
  TT <- rbinom(n,1,ps.true)
  
  M.t <- function(t){-1 + 0.3*t + 0.25*X[,1]^2 - 0.5*sin(X[,2]) + 0.2*rowSums(X[,seq(3,p,3)]) - 0.15*rowSums(X[,seq(4,p,3)])}
  Y.t.m <- function(t,m){-1 + 0.4*t + 0.15*t*m - 0.5*m + 0.5*sin(X[,1]) - 0.2*X[,2]^2 + 0.3*rowSums(X[,seq(3,p,3)]) - 0.2*rowSums(X[,seq(4,p,3)])}
  
  M <- M.t(TT) + eps3 
  Y <- Y.t.m(TT,M) + eps4 
  
  M.0.true <- M.t(0)
  Y.1.0.true <- Y.t.m(1, M.0.true) #eta
  Y.1.MX.true <- Y.t.m(1, M) #gamma
  
  # ate
  M.0.true <- M.t(0) 
  M.1.true <- M.t(1)
  Y.0.0.true <- Y.t.m(0,M.0.true) #mu0
  Y.1.1.true <- Y.t.m(1,M.1.true) #mu1
  
  f_M_0x <- dnorm(M, mean=M.t(0),sd=std) # f(m|T=0,X)
  f_M_1x <- dnorm(M, mean=M.t(1),sd=std) # f(m|T=1,X)
  b1.true <- 1/(1-ps.true)
  b2.true <- f_M_0x/f_M_1x/ps.true
  kappa.true <- Y.1.0.true + (1-TT)*b1.true*(Y.1.MX.true-Y.1.0.true) + TT*b2.true*(Y-Y.1.MX.true)
  
  u <- X[,1]*0.5-X[,1]^2+X[,2]+X[,3]^2+X[,4]+2
  mis.p1 <- 1/(1+3*exp(-u)) # pr(1,1)
  mis.p0 <- 1/(3+exp(u)) # pr(1,0), pr(0,1), pr(0,0)
  summary(mis.p1)
  
  R1 <- R2 <- vector(length = n)
  for (i in 1:n) {
    ind = sample(4,1,prob=c(mis.p1[i],mis.p0[i],mis.p0[i],mis.p0[i]))
    R1[i] <- ifelse(ind %in% c(1,2),1,0)
    R2[i] <- ifelse(ind %in% c(1,3),1,0)
  }

  X1.obs <- R1*X[,1] + (1-R1)*(-99)
  X2.obs <- R2*X[,2] + (1-R2)*(-99)
  X.com <- X
  X.obs <- cbind(X1.obs,X2.obs,X[,-c(1,2)])
  
  
  
  return(list(Z=Z,R1=R1,R2=R2,Rb=R1*R2,TT=TT,M=M,Y=Y,X.true=X.com,X=X.com,X.obs=X.obs,
              eta=Y.1.0.true,gamma=Y.1.MX.true,Y.0.0.true=Y.0.0.true,Y.1.1.true=Y.1.1.true,
              ps.true=ps.true, b1.true=b1.true, f_M_0x=f_M_0x, f_M_1x=f_M_1x, b2.true=b2.true,kappa.true=kappa.true,
              mis.p1=mis.p1,mis.p0=mis.p0,ps.true=ps.true))
  
  
}



compute_true <- function(p.l=10,n=200000){
  true.value <- list()
  for (p in p.l) {
    data.comp <- datagen(n=n,seed=123,p=p)
    theta.true <- mean(data.comp$eta)
    alpha0.true <- mean(data.comp$Y.0.0.true)
    alpha1.true <- mean(data.comp$Y.1.1.true)
    nde.true = theta.true-alpha0.true
    nie.true = alpha1.true-theta.true
    ate.true = alpha1.true-alpha0.true
    true.value[[length(true.value)+1]] <- c(theta.true,alpha1.true,alpha0.true,nie.true,nde.true,ate.true)
  }
  true.value
  
}
