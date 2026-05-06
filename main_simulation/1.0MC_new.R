# Rscript --vanilla --max-connections=256 1.0MC_new.R > output.log 2>&1 &
rm(list=ls())

source("1.1Data_gen.R")
source("1.2Estimation.R")
source("1.3functions.R")
library(MASS)
library(stats)
library(foreach)
library(doParallel)
library(ggplot2)
library(tidyr)
library(dplyr)
library(mice)


# Current time
current_time <- Sys.time()
print(current_time)


# Monte Carlo setting
seed = sample(9999999,1)
# seed = 7947626 # for replication of the results in paper
numt = c(1000,2000)
J = 500

# true value
data.comp <- datagen(2000000,mis.mech = mis.mech)
theta.true <- mean(data.comp$eta)
alpha0.true <- mean(data.comp$Y.0.0.true)
alpha1.true <- mean(data.comp$Y.1.1.true)
nde.true = theta.true-alpha0.true
nie.true = alpha1.true-theta.true
ate.true = alpha1.true-alpha0.true
print(paste0("obs:",mean(data.comp$Rb)))
rm(data.comp)

# print information
print(numt)
print(paste0("seed=",seed))
print(paste0("J=",J))


# parallel computation settings
ncores = detectCores()
coremax = 120
cl = makeCluster(min(ncores,coremax))
registerDoParallel(cl)


est.final <- foreach(rep=1:J, .combine = rbind, .packages = c("mice","MASS")) %dopar% {
  
  res.df <- NULL
  for (i in seq_along(numt)){
    data.input<-datagen(numt[i],seed+rep)
    or <- oracle(data.input)
    cc <- oracle(data.input,method="ignore")
    mi <- MI(data.input)
    pp <- MEMC_DGP1(data.input,var.cal=T)
    
    res.or  <- data.frame(method="Oracle",rep=rep,n=numt[i],
                          estimand=c("theta","alpha1","alpha0","nie","nde","ate"),
                          true.value=c(theta.true,alpha1.true,alpha0.true,nie.true,nde.true,ate.true),
                          est=c(or$theta.hat,or$alpha1.hat,or$alpha0.hat,or$nie.hat,or$nde.hat,or$ate.hat),
                          std=c(or$std.y10,or$std.y11,or$std.y00,or$std.nie,or$std.nde,or$std.ate))
    res.cc  <- data.frame(method="CCA",rep=rep,n=numt[i],
                          estimand=c("theta","alpha1","alpha0","nie","nde","ate"),
                          true.value=c(theta.true,alpha1.true,alpha0.true,nie.true,nde.true,ate.true),
                          est=c(cc$theta.hat,cc$alpha1.hat,cc$alpha0.hat,cc$nie.hat,cc$nde.hat,cc$ate.hat),
                          std=c(cc$std.y10,cc$std.y11,cc$std.y00,cc$std.nie,cc$std.nde,cc$std.ate))
    res.mi  <- data.frame(method="MI",rep=rep,n=numt[i],
                          estimand=c("theta","alpha1","alpha0","nie","nde","ate"),
                          true.value=c(theta.true,alpha1.true,alpha0.true,nie.true,nde.true,ate.true),
                          est=c(mi$theta.hat,mi$alpha1.hat,mi$alpha0.hat,mi$nie.hat,mi$nde.hat,mi$ate.hat),
                          std=c(mi$std.y10,mi$std.y11,mi$std.y00,mi$std.nie,mi$std.nde,mi$std.ate))
    res.pp  <- data.frame(method="SIO",rep=rep,n=numt[i],
                          estimand=c("theta","alpha1","alpha0","nie","nde","ate"),
                          true.value=c(theta.true,alpha1.true,alpha0.true,nie.true,nde.true,ate.true),
                          est=c(pp$theta.hat,pp$alpha1.hat,pp$alpha0.hat,pp$nie.hat,pp$nde.hat,pp$ate.hat),
                          std=c(pp$std.y10,pp$std.y11,pp$std.y00,pp$std.nie,pp$std.nde,pp$std.ate),
                          ss.01=pp$ss.01,ss.10=pp$ss.10,ss.00=pp$ss.00,
                          ss.gamma=pp$ss.gamma,ss.eta=pp$ss.eta,ss.Y11=pp$ss.Y11,ss.Y00=pp$ss.Y00)
   
    res.df <- dplyr::bind_rows(res.df,res.or,res.cc,res.mi,res.pp)
  }
  res.df
  
} 


stopCluster(cl)
stopImplicitCluster()

# save data
save.image("res_DGP1.rdata")

print(Sys.time())
print(Sys.time()-current_time)

source("1.5table.R")


