# rm(list=ls())

## load data
load("trimming.rdata")
load("SIO_bias_se.rdata")

est.summary <- est.final %>%
  mutate(estimand = factor(estimand, levels=c("theta","nie","nde","ate"),
                           labels = c("MF", "NIE", "NDE", "ATE"))) %>% 
  group_by(threshold, estimand) %>%
  summarise(
    bias = abs(mean(est) - unique(true.value)),
    se = sd(est),
    .groups = "drop"
  ) 
print(data.frame(est.summary))


estimands <- unique(est.summary$estimand)


SIO_1000 <- subset(SIO_summary, n == 1000)
bias_lines <- setNames(abs(SIO_1000$bias), SIO_1000$estimand)   
se_lines   <- setNames(SIO_1000$se, SIO_1000$estimand)



pdf("plt_trim.pdf", height=6 , width=15,pointsize = 17)
# pdf(paste0(formatted_time,".pdf"), height=6 , width=15,pointsize = 17)

par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))   

for (est in estimands) {
  subdf <- subset(est.summary, estimand == est)
  subdf <- subdf[order(-subdf$threshold), ]
  
  plot(subdf$threshold, abs(subdf$bias), type = "b", pch = 19,
       xlab = "Threshold", ylab = "Bias",
       main = paste("Bias of", est),
       xlim = range(est.summary$threshold),
       ylim = range(c(abs(est.summary$bias), bias_lines[est],-0.1), na.rm = TRUE),
       cex.axis = 1.2, cex.lab = 1.6, cex.main = 1.5, mgp = c(2, 0.3, 0))
  
  abline(h = bias_lines[est], col = "blue", lty = 2, lwd = 1.2)  
}

for (est in estimands) {
  subdf <- subset(est.summary, estimand == est)
  subdf <- subdf[order(-subdf$threshold), ]
  
  plot(subdf$threshold, subdf$se, type = "b", pch = 19,
       xlab = "Threshold", ylab = "SE",
       main = paste("SE of", est),
       xlim = range(est.summary$threshold),
       ylim = range(c(est.summary$se, se_lines[est],-0.1), na.rm = TRUE),
       cex.axis = 1.2, cex.lab = 1.6, cex.main = 1.5, mgp = c(2, 0.3, 0))
  
  abline(h = se_lines[est], col = "blue", lty = 2, lwd = 1.2)  
}
dev.off()
