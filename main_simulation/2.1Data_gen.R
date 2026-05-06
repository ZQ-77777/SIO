
datagen <- function(n,seed=NA,w=0.6,...){
  
  if(is.na(seed)){seed=sample(9999999,1)}
  set.seed(seed)
  
  eps1 <- rnorm(n)
  eps2 <- rnorm(n)
  eps3 <- rnorm(n)
  eps4 <- rnorm(n)
  eps5 <- rnorm(n)
  eps1p <- rnorm(n)
  eps2p <- rnorm(n)
  
  X1 <- pnorm(eps1)
  X2 <- pnorm(eps2)
  X3 <- rbinom(n,1,0.5)

  
  ps.true <- expit(0.1+X1^2-X2+0.2*X3) # propensity score of treatment
  TT <- rbinom(n,1,ps.true)
  
  M.t <- function(t){ -1 + t - t*X1^2 + 4*X1^2 - sin(X2) + X3}
  Y.t.m <- function(t,m){-1 + t + 3*t*m - 1.5*m + t*X1+ 5*sin(X1) - 2*X2^2 + X3 }
  
  M <- M.t(TT) + eps3 
  Y <- Y.t.m(TT,M) + eps4 
  
  M.0.true <- M.t(0)
  Y.1.0.true <- Y.t.m(1, M.0.true) #eta
  Y.1.MX.true <- Y.t.m(1, M) #gamma
  
  # ate
  M.0.true <- M.t(0) 
  M.1.true <- M.t(1)
  Y.0.0.true <- Y.t.m(0,M.0.true) #mu0
  Y.1.1.true <- Y.t.m(1,M.1.true) #mu0
  
  u <-  1.2+X1 + X2^2 - X3 - 0.5 * TT + TT*X1^2 
  mis.p1 <- 1/(1+3*exp(-u)) # pr(1,1)
  mis.p0 <- 1/(3+exp(u)) # pr(1,0), pr(0,1), pr(0,0)
  
  R1 <- R2 <- vector(length = n)
  for (i in 1:n) {
    ind = sample(4,1,prob=c(mis.p1[i],mis.p0[i],mis.p0[i],mis.p0[i]))
    R1[i] <- ifelse(ind %in% c(1,2),1,0)
    R2[i] <- ifelse(ind %in% c(1,3),1,0)
  }
  

  
  X1.obs <- R1*X1 + (1-R1)*(-99)
  X2.obs <- R2*X2 + (1-R2)*(-99)
  X.com <- cbind(X1,X2,X3)
  X.obs <- cbind(X1.obs,X2.obs,X3)
  
  
  
  return(list(R1=R1,R2=R2,Rb=R1*R2,TT=TT,M=M,Y=Y,X.true=X.com,X=X.com,X.obs=X.obs,
              eta=Y.1.0.true,gamma=Y.1.MX.true,Y.0.0.true=Y.0.0.true,Y.1.1.true=Y.1.1.true,
              mis.p1=mis.p1,mis.p0=mis.p0,ps.true=ps.true))
  
  
}

