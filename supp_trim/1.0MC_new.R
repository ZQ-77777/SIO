# Rscript --vanilla --max-connections=256 1.0MC_new.R > output.log 2>&1 &
rm(list=ls())

source("1.1Data_gen.R")
source("1.2Estimation.R")
source("1.3functions.R")
library(stats)
library(foreach)
library(doParallel)
library(ggplot2)
library(tidyr)
library(dplyr)
library(mice)
library(MASS)



# Current time
current_time <- Sys.time()
print(current_time)
formatted_time <- format(current_time, "%y%m%d_%H%M%S")



# Monte Carlo setting
seed = sample(9999999,1)
# seed = 4727212 # for replication of the results in paper
numt = 1000
thresholds <- c(0.3,0.25,0.2,0.15,0.1,0.05,0)
n.th <- length(thresholds)
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
coremax = 250
cl = makeCluster(min(ncores,coremax))
registerDoParallel(cl)




est.final <- foreach(rep=1:J, .combine = rbind, .packages = c("mice","MASS")) %dopar% {
  
  res.df <- NULL
  data.input<-datagen(numt,seed+rep)
  pp <- TRIM(data.input,thresholds)
  res.pp  <- data.frame(method="trimming",rep=rep,n=numt,
                        estimand=c("theta","nie","nde","ate"),
                        true.value=c(theta.true,nie.true,nde.true,ate.true),
                        est=pp$est,n.sel<-rep(pp$n,each=4),threshold=rep(thresholds,each=4))
 
  res.df <- dplyr::bind_rows(res.df,res.pp)
  res.df
  
} 


stopCluster(cl)
stopImplicitCluster()

# save data
save.image("trimming.rdata")

print(Sys.time())
print(Sys.time()-current_time)

source("1.4summary.R")


