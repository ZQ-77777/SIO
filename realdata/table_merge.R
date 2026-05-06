rm(list=ls())
library(dplyr)
load("SIO.rdata")
load("SIO_bs.rdata")
load("MI_CCA.rdata")


sd.merge <- est.final %>%
  mutate(method=factor(method,levels=c("SIO1","SIO2"))) %>% 
  mutate(estimand=factor(estimand,levels=c("NIE","NDE","ATE"))) %>% 
  group_by(estimand,method) %>%
  summarise(se = sd(est),.groups = "drop") %>% pull(se)


final.df <- data.frame(Estimand = rep(c("NIE","NDE","ATE"),each=4),Methods = rep(c("SIO1","SIO2","MI","CCA"),3),
                       PointEstimation=NA, CI=NA)


p.est <- c(nie.hat.sio1,nie.hat.sio2, med.mi$d1,med.cca$d1,nde.hat.sio1,nde.hat.sio2,med.mi$z0,med.cca$z0,ate.hat.sio1,ate.hat.sio2,med.mi$tau.coef,med.cca$tau.coef)
final.df$PointEstimation = round(p.est,3)

final.df$CI[c(1,2,5,6,9,10)] <- paste0("[",round(p.est[c(1,2,5,6,9,10)]-1.96*sd.merge,3),",",
                                              round(p.est[c(1,2,5,6,9,10)]+1.96*sd.merge,3),"]")
final.df[3,4] <- paste0("[",round(med.mi$d1.ci[1],3),",",round(med.mi$d1.ci[2],3),"]")
final.df[7,4] <- paste0("[",round(med.mi$z0.ci[1],3),",",round(med.mi$z0.ci[2],3),"]")
final.df[11,4] <- paste0("[",round(med.mi$tau.ci[1],3),",",round(med.mi$tau.ci[2],3),"]")

final.df[4,4] <- paste0("[",round(med.cca$d1.ci[1],3),",",round(med.cca$d1.ci[2],3),"]")
final.df[8,4] <- paste0("[",round(med.cca$z0.ci[1],3),",",round(med.cca$z0.ci[2],3),"]")
final.df[12,4] <- paste0("[",round(med.cca$tau.ci[1],3),",",round(med.cca$tau.ci[2],3),"]")


print(final.df)


