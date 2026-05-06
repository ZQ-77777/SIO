# Rscript --vanilla MC_dml.R > output.log 2>&1 &
rm(list=ls())
library(MASS)
library(stats)
library(foreach)
library(doParallel)
select <- dplyr::select
library(glmnet)
library(dplyr)

# gen.method = "sparse"
gen.method = "dense"
if (gen.method=="dense") {
  source("data_gen_dense.R")
}else{
  source("data_gen_sparse.R")
}
print(gen.method)

# method = "lasso"
method = "nn"
if (method == "nn" & gen.method == "dense") {
  source("Estimation_nn_dense.R")
} else if (method == "nn" & gen.method == "sparse") {
  source("Estimation_nn_sparse.R")
} else {
  source("Estimation_lasso.R")
}
print(method)

# Current time
current_time <- Sys.time()
print(current_time)
formatted_time <- format(current_time, "%y%m%d_%H%M%S")

# Monte Carlo setting
seed = sample(9999999,1)
# seed = 2560340 # for replication of the results of nn,sparse,p=50 in the paper
# seed = 4160527 # for replication of the results of nn,sparse,p=100 in the paper
# seed = 2521789 # for replication of the results of lasso,sparse,p=50 in the paper
# seed = 4250155 # for replication of the results of lasso,sparse,p=100 in the paper
# seed = 2997824 # for replication of the results of nn,dense,p=50 in the paper
# seed = 8227973 # for replication of the results of nn,dense,p=100 in the paper
# seed = 3227207 # for replication of the results of lasso,dense,p=50 in the paper
# seed = 8757510 # for replication of the results of lasso,dense,p=100 in the paper
numt = 3000
# p.l <- seq(10,50,10)
p.l <- 50
J = 500


# true value
true.value <- compute_true(p.l)
print(true.value)

# print information
print(numt)
print(paste0("seed=",seed))
print(paste0("J=",J))
print(paste0("p=",p.l))

# parallel computation settings
cl = makeCluster(40)
registerDoParallel(cl)

est.final <- foreach(rep=1:J, .combine = rbind, .packages = c("torch","caret","MASS","dplyr", "glmnet")) %dopar% {
  
  res.df <- NULL
  for (i in seq_along(p.l)){
    
    res.pp <- tryCatch({
      data.input <- datagen(numt, seed+rep, p=p.l[i])
      pp <- est_dml(data.input, seed+rep)
      
      data.frame(method="DML",rep=rep,n=numt,p=p.l[i],
                 estimand=c("theta","alpha1","alpha0","nie","nde","ate"),
                 true.value=true.value[[i]],
                 est=c(pp$theta,pp$alpha1,pp$alpha0,pp$nie,pp$nde,pp$ate),
                 sd=c(pp$theta_sd,pp$alpha1_sd,pp$alpha0_sd,pp$nie_sd,pp$nde_sd,pp$ate_sd))
      
    }, error = function(e) {
      
      data.frame(method="DML",rep=rep,n=numt,p=p.l[i],
                 estimand=c("theta","alpha1","alpha0","nie","nde","ate"),
                 true.value=true.value[[i]],
                 est=rep(NA,6),
                 sd=rep(NA,6))
    })
    
    res.df <- dplyr::bind_rows(res.df, res.pp)
  }
  res.df
}

summary_res <- est.final %>%
  group_by(estimand, n, p) %>%
  summarise(
    mean_est = mean(est, na.rm = TRUE),
    mean_se  = mean(sd / sqrt(n), na.rm = TRUE),
    true_val = mean(true.value, na.rm = TRUE),
    .groups = "drop"
  )


alpha <- 0.05
z_val <- qnorm(1 - alpha / 2)

est.cover <- est.final %>%
  mutate(
    se = sd / sqrt(n),
    lower = est - z_val * se,
    upper = est + z_val * se,
    cover = (true.value >= lower & true.value <= upper)
  )


cover_res <- est.cover %>%
  group_by(estimand, n, p) %>%
  summarise(
    coverage = mean(cover, na.rm = TRUE),
    .groups = "drop"
  )


final_summary <- left_join(summary_res, cover_res, by = c("estimand", "n", "p"))
print(final_summary)

alpha <- 0.05
z_val <- qnorm(1 - alpha / 2)  

theta_summary_raw <- est.final %>%
  filter(estimand == "theta") %>%
  mutate(
    se = sd / sqrt(n), 
    lower = est - z_val * se,  
    upper = est + z_val * se,  
    cover = as.numeric(true.value >= lower & true.value <= upper)  
  ) %>%
  select(estimand, n, p, est, sd, se, lower, upper, cover)


# save data
save.image(paste0("res_",method,"_",gen.method,"_",p.l[1],"_",formatted_time,".RData"))

print(Sys.time())
print(Sys.time()-current_time)

# 
stopCluster(cl)
stopImplicitCluster()
