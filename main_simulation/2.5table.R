rm(list=ls())
library(dplyr)
library(tidyr)


## load data
load("res_DGP2.rdata")

est.summary <- est.final %>%
  mutate(estimand = factor(estimand, levels=c("theta","alpha0","alpha1","nie","nde","ate"))) %>% 
  mutate(method = factor(method, levels=c("Oracle","SIO","MI","CCA"))) %>% 
  group_by(method, n, estimand) %>%
  summarise(
    bias = mean(est) - unique(true.value),
    se = sd(est),
    coverage =  mean((est - 1.96 * std <= unique(true.value)) &  (est + 1.96 * std >= unique(true.value))),
    CIlen=mean(std*1.96*2),
    .groups = "drop"
  ) 
print(data.frame(est.summary))


bias_sd_table <- est.summary %>%
  mutate(n_estimand = paste0("n", n, "_", estimand),
         bias = sprintf("%.2f", bias),
         se = sprintf("%.2f", se),
         bias_se = paste0(bias,"(",se,")")) %>% 
  filter(estimand %in% c("nie","nde","ate","theta")) %>% 
  select(method, n_estimand, bias_se) %>%
  pivot_wider(names_from = n_estimand, values_from = bias_se)
print(data.frame(bias_sd_table))


CI_table <- est.summary %>%
  mutate(n_estimand = paste0("n", n, "_", estimand),
         coverage = sub("^0", "", sprintf("%.3f", coverage)),
         CIlen = sprintf("%.2f", CIlen),
         cov_len = paste0(coverage,"(",CIlen,")")) %>%
  filter(estimand %in% c("nie","nde","ate","theta")) %>% 
  select(method, n_estimand, cov_len) %>%
  pivot_wider(names_from = n_estimand, values_from = cov_len)
print(data.frame(CI_table))



library(knitr)
kable(bias_sd_table,format = "latex",escape = FALSE)
kable(CI_table,format = "latex",escape = FALSE)
