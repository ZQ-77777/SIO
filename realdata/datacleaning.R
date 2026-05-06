# data: CFPS 2020 and 2022
rm(list=ls())
library(conflicted)
library(haven)
library(dplyr)
library(reshape)
conflict_prefer("rename","reshape")
removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}

cfps2020person <- read_sas("SAS2020/cfps2020person_202306.sas7bdat") 
cfps2022person <- read_sas("SAS2022/cfps2022person_202410.sas7bdat")
cfps2020famecon <- read_sas("SAS2020/cfps2020famecon_202306.sas7bdat")
cfps2020famconf <- read_sas("SAS2020/cfps2020famconf_202306.sas7bdat") 
cfps2020crossyearid <- read_sas("SAS2020/cfps2020crossyearid_202312.sas7bdat")

cross <- cfps2020crossyearid[cfps2020crossyearid$entrayear<=2020,] 
cross <- cfps2020crossyearid[cfps2020crossyearid$deceased==0,] 
cross2020 <- cross[, c("pid","birthy","gender","marriage_20")]

cfps2020person <- rename(cfps2020person,c('PID' = "pid", 'FID20' = "fid20"))
cfps2022person <- rename(cfps2022person,c('PID' = "pid", 'FID20' = "fid20"))
cfps2020famecon <- cfps2020famecon %>% mutate(engel=FOOD/expense)

data_person <- cfps2020person[,c("pid","QG406","emp_income","QP201")] 
data_person <- data_person %>%
  left_join(cfps2022person %>% dplyr::select(pid, CESD8, QM2016), by = "pid")

data_family <- cfps2020famconf[, c("fid20","pid","familysize20")]
data_eco <- cfps2020famecon[, c("fid20","FR1","engel")]

merge <- merge(data_family, data_person, by = "pid",all=F)
merge <- merge(merge, cross2020, by = "pid",all=F)
merge <- merge(merge, data_eco, by ="fid20", all=F)

merge <- rename(merge,c('familysize20' = "familysize", 'QG406' = "jobsatisfication", 
                        "emp_income"="income","FR1"="house","CESD8"="depression","marriage_20"="marriage",
                        "QM2016"="wellbeing","QP201"="health"))
data_1 <-  merge[merge$familysize > 0
                 & merge$gender %in% c(0, 1)
                 & merge$marriage %in% c(1:5)    
                 & merge$house  %in% c(1, 5)
                 & merge$jobsatisfication %in% c(1:5)  
                 & merge$depression >= 0
                 & merge$wellbeing %in% c(0:10)
                 & merge$health %in% c(1:5)
                 & merge$engel < 1
                 & merge$engel > 0
                 ,]
data_1$income[data_1$income <= 0] <- NA

mediation <- data_1[complete.cases(data_1[-5]),] 
mediation$age <- 2020- mediation$birthy
mediation <- mediation[, -c(1,2,9)]
mediation <- mediation %>% dplyr::filter(age>=35&age<=60) 
mediation$marriage <- ifelse(mediation$marriage==2,1,0)
mediation$house <- ifelse(mediation$house==5,0,1)
mean(!is.na(mediation$income))
save(mediation, file = "mediation.rdata")

