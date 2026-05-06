# rm(list=ls())
options(warn = -1)
library(tidyr)
library(dplyr)
library(latex2exp)

## load data
load("MAR2.rdata")
plot.method = c("Oracle","SIO","MI","CCA")
est.final$method = factor(est.final$method,levels=plot.method)


pdf("plt_MAR.pdf", height=2.8 , width=12,pointsize = 17)
par(mfrow = c(1, 4), mar=c(2.8, 3.2, 2, 2))


# theta
est.final %>% filter(estimand=="theta") %>%
  boxplot(est~n+method,data=.,at=c(1,2,4,5,7,8,10,11),
          ylim=c(theta.true-1.1,theta.true+1.1), cex.axis = 0.8,
          col = c("white","gray"),xaxt = "n",outline=F,
          xlab="",ylab="",las = 2)
abline(h = theta.true, col = "red", lwd = 1.5)
title(main = "(a) MF",adj=0,cex.main=1,line = 1)
title(ylab=TeX("Estimates of $\\theta$"),cex.lab = 1,line = 2)
legend("topleft", fill = c("white","gray"), legend = as.character(numt), horiz = T,cex=0.8)
axis(1, at = c(1.5, 4.5, 7.5, 10.5), tck=-0.04, labels = plot.method, cex.axis = 0.85, line = 0)



# NIE
est.final %>% filter(estimand=="nie") %>%
  boxplot(est~n+method,data=.,at=c(1,2,4,5,7,8,10,11),
          ylim=c(nie.true-1,nie.true+1), cex.axis = 0.8,
          col = c("white","gray"),xaxt = "n",outline=F,
          xlab="",ylab="",las = 2)
abline(h = nie.true, col = "red", lwd = 1.5)
title(main = "(b) NIE",adj=0,cex.main=1,line = 1)
title(ylab=TeX("Estimates of $\\xi^{(1)}$"),cex.lab = 1,line = 2)
legend("topleft", fill = c("white","gray"), legend = as.character(numt), horiz = T,cex=0.8)
axis(1, at = c(1.5, 4.5, 7.5, 10.5), tck=-0.04, labels = plot.method, cex.axis = 0.85, line = 0)



# NDE
est.final %>% filter(estimand=="nde") %>%
  boxplot(est~n+method,data=.,at=c(1,2,4,5,7,8,10,11),
          ylim=c(nde.true-1.5,nde.true+1.5), cex.axis = 0.8,
          col = c("white","gray"),xaxt = "n",outline=F,
          xlab="",ylab="",las = 2)
abline(h = nde.true, col = "red", lwd = 1.5)
title(main = "(c) NDE",adj=0,cex.main=1,line = 1)
title(ylab=TeX("Estimates of $\\zeta^{(0)}$"),cex.lab = 1,line = 2)
legend("topleft", fill = c("white","gray"), legend = as.character(numt), horiz = T,cex=0.8)
axis(1, at = c(1.5, 4.5, 7.5, 10.5), tck=-0.04, labels = plot.method, cex.axis = 0.85, line = 0)



# boxplot of ATE
est.final %>% filter(estimand=="ate") %>%
  boxplot(est~n+method,data=.,at=c(1,2,4,5,7,8,10,11),
          ylim=c(ate.true-1.4,ate.true+1.4), cex.axis = 0.8,
          col = c("white","gray"),xaxt = "n",outline=F,
          xlab="",ylab="",las = 2)
abline(h = ate.true, col = "red", lwd = 1.5)
title(main = "(d) ATE",adj=0,cex.main=1,line = 1)
title(ylab=TeX("Estimates of $\\tau$"),cex.lab = 1,line = 2)
legend("topleft", fill = c("white","gray"), legend = as.character(numt), horiz = T,cex=0.8)
axis(1, at = c(1.5, 4.5, 7.5, 10.5), tck=-0.04, labels = plot.method, cex.axis = 0.85, line = 0)



dev.off()





