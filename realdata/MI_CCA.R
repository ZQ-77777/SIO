rm(list=ls())
library(MASS)
library(dplyr)
library(mice)
library(mediation)

# Current time
print(Sys.time())

# data preparation
load("mediation.rdata")
R <- ifelse(is.na(mediation$income),0,1)
mediation$marriage.status <- ifelse(mediation$marriage==1,1,0)
X.com <- mediation[,c("familysize","health","gender","house","age","marriage.status")]
X.mis <- mediation$income
Y <- mediation$depression
M <- mediation$wellbeing
TT <- ifelse(mediation$jobsatisfication>3,1,0)
X.obs <- cbind(X.com,X.mis)
n=length(TT)


# MI
imp <- mice(cbind(TT,M,Y,X.obs),quiet=T)
X.imp <- as.matrix(complete(imp)[,10])
data.df = data.frame(Y,TT,M,X.com,income=X.imp)

med.fit <- lm(M ~ TT + (familysize + familysize^2+ familysize^3 + health + health^2 + health^3 + gender + house + house^2 + house^3 + age + age^2 + age^3 + marriage.status +  income + income^2+ income^3), 
              data = data.df)
out.fit <- lm(Y ~  M+TT + (familysize + familysize^2+ familysize^3 + health + health^2 + health^3 + gender + house + house^2 + house^3 + age + age^2 + age^3 + marriage.status +  income + income^2+ income^3),
              data = data.df)

med.mi <- mediate(med.fit, out.fit, treat = "TT", mediator = "M",
                  robustSE = TRUE, sims = 500)
summary(med.mi)

# CCA
med.fit <- lm(M ~ TT + familysize + familysize^2+ familysize^3 + health + health^2 + health^3 + gender + house + house^2 + house^3 + age + age^2 + age^3 + marriage.status +  income + income^2+ income^3, 
              data = data.df,subset=(R==1))
out.fit <- lm(Y ~ TT + M + familysize + familysize^2+ familysize^3 + health + health^2 + health^3 + gender + house + house^2 + house^3 + age + age^2 + age^3 + marriage.status +  income + income^2+ income^3,
              data = data.df,subset=(R==1))
med.cca <- mediate(med.fit, out.fit, treat = "TT", mediator = "M",
                   robustSE = TRUE, sims = 500)
summary(med.cca)


save(med.mi, med.cca, file = "MI_CCA.rdata")
