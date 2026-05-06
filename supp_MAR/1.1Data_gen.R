
datagen <- function(n,seed=NA,w=0.6,...){
  
  if(is.na(seed)){seed=sample(9999999,1)}
  set.seed(seed)
  
  eps1 <- rnorm(n)
  eps2 <- rnorm(n)
  
  X1 <- runif(n)
  X2 <- runif(n)
  X3 <- rbinom(n,1,0.5)

  
  ps.true <- expit(0.2+X1^2-X2+0.2*X3) # propensity score of treatment
  TT <- rbinom(n,1,ps.true)
  
  eps3 <- rnorm(n)
  eps4 <- rnorm(n)
  
  M.t <- function(t){ -1 + t + 4*X1^2 - sin(X2) + X3}
  Y.t.m <- function(t,m){-1 + t + 3*t*m - 1.5*m + 5*sin(X1) - 2*X2^2 + X3 }
  
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
  
  mis.pr <- expit(1-0.5*TT+2*X2^2-X3+0.5*M-0.2*Y)
  R <- rbinom(n,1,mis.pr)
  
  X1.obs <- R*X1 + (1-R)*(-9999)
  X.com <- cbind(X1,X2,X3)
  X.obs <- cbind(X1.obs,X2,X3)
  
  
  return(list(R=R,TT=TT,M=M,Y=Y,X.true=X.com,X=X.com,X.obs=X.obs,
              eta=Y.1.0.true,gamma=Y.1.MX.true,Y.0.0.true=Y.0.0.true,Y.1.1.true=Y.1.1.true,
              mis.pr=mis.pr,ps.true=ps.true))
  
  
}

#