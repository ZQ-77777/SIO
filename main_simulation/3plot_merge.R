# rm(list=ls())
options(warn = -1)
library(tidyr)
library(dplyr)
library(latex2exp)



pdf("plt_merge.pdf", height=4 , width=12,pointsize = 17)
par(mfrow = c(2, 4), mar=c(1.5, 3.2, 1.5, 0.5))


# ----------------DGP1------------------
## load data
load("res_DGP1.rdata")
plot.method = c("Oracle","SIO","MI","CCA")
est.final$method = factor(est.final$method,levels=plot.method)
# theta
est.final %>% filter(estimand=="theta") %>%
  boxplot(est~n+method,data=.,at=c(1,2,4,5,7,8,10,11),
          ylim=c(theta.true-1.2/1.5,theta.true+1.2), cex.axis = 0.8,
          col = c("white","gray"),xaxt = "n",outline=F,
          xlab="",ylab="",las = 2)
abline(h = theta.true, col = "red", lwd = 1.5)
title(main = "(a) DGP1, MF",adj=0,cex.main=1,line = 0.7)
title(ylab=TeX("Estimates of $\\theta$"),cex.lab = 0.8,line = 2)
legend("topleft", fill = c("white","gray"), legend = as.character(numt), horiz = T,cex=0.7)
axis(1, at = c(1.5, 4.5, 7.5, 10.5), labels = FALSE, tck = -0.04)
mtext(side = 1, at = c(1.5, 4.5, 7.5, 10.5), text = plot.method,line = 0.3, cex = 0.45)


# NIE
est.final %>% filter(estimand=="nie") %>%
  boxplot(est~n+method,data=.,at=c(1,2,4,5,7,8,10,11),
          ylim=c(nie.true-0.8/1.5,nie.true+0.8), cex.axis = 0.8,
          col = c("white","gray"),xaxt = "n",outline=F,
          xlab="",ylab="",las = 2)
abline(h = nie.true, col = "red", lwd = 1.5)
title(main = "(b) DGP1, NIE",adj=0,cex.main=1,line = 0.7)
title(ylab=TeX("Estimates of $\\xi^{(1)}$"),cex.lab = 0.8,line = 2)
legend("topleft", fill = c("white","gray"), legend = as.character(numt), horiz = T,cex=0.7)
axis(1, at = c(1.5, 4.5, 7.5, 10.5), labels = FALSE, tck = -0.04)
mtext(side = 1, at = c(1.5, 4.5, 7.5, 10.5), text = plot.method,line = 0.3, cex = 0.45)


# NDE
est.final %>% filter(estimand=="nde") %>%
  boxplot(est~n+method,data=.,at=c(1,2,4,5,7,8,10,11),
          ylim=c(nde.true-1.6/1.5,nde.true+1.6), cex.axis = 0.8,
          col = c("white","gray"),xaxt = "n",outline=F,
          xlab="",ylab="",las = 2)
abline(h = nde.true, col = "red", lwd = 1.5)
title(main = "(c) DGP1, NDE",adj=0,cex.main=1,line = 0.7)
title(ylab=TeX("Estimates of $\\zeta^{(0)}$"),cex.lab = 0.8,line = 2)
legend("topleft", fill = c("white","gray"), legend = as.character(numt), horiz = T,cex=0.7)
axis(1, at = c(1.5, 4.5, 7.5, 10.5), labels = FALSE, tck = -0.04)
mtext(side = 1, at = c(1.5, 4.5, 7.5, 10.5), text = plot.method,line = 0.3, cex = 0.45)


# boxplot of ATE
est.final %>% filter(estimand=="ate") %>%
  boxplot(est~n+method,data=.,at=c(1,2,4,5,7,8,10,11),
          ylim=c(ate.true-1.5/1.5,ate.true+1.5), cex.axis = 0.8,
          col = c("white","gray"),xaxt = "n",outline=F,
          xlab="",ylab="",las = 2)
abline(h = ate.true, col = "red", lwd = 1.5)
title(main = "(d) DGP1, ATE",adj=0,cex.main=1,line = 0.7)
title(ylab=TeX("Estimates of $\\tau$"),cex.lab = 0.8,line = 2)
legend("topleft", fill = c("white","gray"), legend = as.character(numt), horiz = T,cex=0.7)
axis(1, at = c(1.5, 4.5, 7.5, 10.5), labels = FALSE, tck = -0.04)
mtext(side = 1, at = c(1.5, 4.5, 7.5, 10.5), text = plot.method,line = 0.3, cex = 0.45)



# -----------------DGP2--------------------------
load("res_DGP2.rdata")
plot.method = c("Oracle","SIO","MI","CCA")
est.final$method = factor(est.final$method,levels=plot.method)
# theta
est.final %>% filter(estimand=="theta") %>%
  boxplot(est~n+method,data=.,at=c(1,2,4,5,7,8,10,11),
          ylim=c(theta.true-1.1/1.2,theta.true+1.1), cex.axis = 0.8,
          col = c("white","gray"),xaxt = "n",outline=F,
          xlab="",ylab="",las = 2)
abline(h = theta.true, col = "red", lwd = 1.5)
title(main = "(e) DGP2, MF",adj=0,cex.main=1,line = 0.7)
title(ylab=TeX("Estimates of $\\theta$"),cex.lab = 0.8,line = 2)
legend("topleft", fill = c("white","gray"), legend = as.character(numt), horiz = T,cex=0.7)
axis(1, at = c(1.5, 4.5, 7.5, 10.5), labels = FALSE, tck = -0.04)
mtext(side = 1, at = c(1.5, 4.5, 7.5, 10.5), text = plot.method,line = 0.3, cex = 0.45)


# NIE
est.final %>% filter(estimand=="nie") %>%
  boxplot(est~n+method,data=.,at=c(1,2,4,5,7,8,10,11),
        ylim=c(nie.true-0.8/1.2,nie.true+0.8), cex.axis = 0.8,
        col = c("white","gray"),xaxt = "n",outline=F,
        xlab="",ylab="",las = 2)
abline(h = nie.true, col = "red", lwd = 1.5)
title(main = "(f) DGP2, NIE",adj=0,cex.main=1,line = 0.7)
title(ylab=TeX("Estimates of $\\xi^{(1)}$"),cex.lab = 0.8,line = 2)
legend("topleft", fill = c("white","gray"), legend = as.character(numt), horiz = T,cex=0.7)
axis(1, at = c(1.5, 4.5, 7.5, 10.5), labels = FALSE, tck = -0.04)
mtext(side = 1, at = c(1.5, 4.5, 7.5, 10.5), text = plot.method,line = 0.3, cex = 0.45)


# NDE
est.final %>% filter(estimand=="nde") %>%
  boxplot(est~n+method,data=.,at=c(1,2,4,5,7,8,10,11),
        ylim=c(nde.true-1.4/1.2,nde.true+1.4), cex.axis = 0.8,
        col = c("white","gray"),xaxt = "n",outline=F,
        xlab="",ylab="",las = 2)
abline(h = nde.true, col = "red", lwd = 1.5)
title(main = "(g) DGP2, NDE",adj=0,cex.main=1,line = 0.7)
title(ylab=TeX("Estimates of $\\zeta^{(0)}$"),cex.lab = 0.8,line = 2)
legend("topleft", fill = c("white","gray"), legend = as.character(numt), horiz = T,cex=0.7)
axis(1, at = c(1.5, 4.5, 7.5, 10.5), labels = FALSE, tck = -0.04)
mtext(side = 1, at = c(1.5, 4.5, 7.5, 10.5), text = plot.method,line = 0.3, cex = 0.45)


# boxplot of ATE
est.final %>% filter(estimand=="ate") %>%
  boxplot(est~n+method,data=.,at=c(1,2,4,5,7,8,10,11),
        ylim=c(ate.true-1.0/1.2,ate.true+1.0), cex.axis = 0.8,
        col = c("white","gray"),xaxt = "n",outline=F,
        xlab="",ylab="",las = 2)
abline(h = ate.true, col = "red", lwd = 1.5)
title(main = "(h) DGP2, ATE",adj=0,cex.main=1,line = 0.7)
title(ylab=TeX("Estimates of $\\tau$"),cex.lab = 0.8,line = 2)
legend("topleft", fill = c("white","gray"), legend = as.character(numt), horiz = T,cex=0.7)
axis(1, at = c(1.5, 4.5, 7.5, 10.5), labels = FALSE, tck = -0.04)
mtext(side = 1, at = c(1.5, 4.5, 7.5, 10.5), text = plot.method,line = 0.3, cex = 0.45)


dev.off()




## plot selected sieve size -------
# -----------------DGP1--------------------------
load("res_DGP1.rdata")
ss.df <- est.final %>% filter(method=="SIO") %>%
  group_by(method,rep,n) %>%
  summarise(ss.01=mean(ss.01),
            ss.10=mean(ss.10),
            ss.00=mean(ss.00),
            ss.gamma=mean(ss.gamma),
            ss.eta=mean(ss.eta),
            ss.Y11=mean(ss.Y11),
            ss.Y00=mean(ss.Y00),.groups = "drop")


ss.vars <- c("ss.10", "ss.01", "ss.00", "ss.gamma", "ss.eta", "ss.Y11", "ss.Y00")

xlab.map <- list(
  ss.01 = expression(J[n]/2),
  ss.10 = expression(J[n]/2),
  ss.00 = expression(J[n]/2),
  ss.gamma = expression(k[n]),
  ss.eta = expression(s[n]),
  ss.Y11 = expression(s[n]),
  ss.Y00 = expression(s[n])
)

title.map <- list(
  ss.10 = "r=(1,0,1)",
  ss.01 = "r=(0,1,1)",
  ss.00 = "r=(0,0,1)",
  ss.gamma = "gamma",
  ss.eta = "eta",
  ss.Y11 = "mu1",
  ss.Y00 = "mu0"
)


pdf("plt_size_DGP1.pdf", height=6, width=15,pointsize = 17)
par(mfrow = c(2, 7), mar = c(5, 3, 2.8, 1)) 

for (n_val in c(1000, 2000)) {
  ss.sub <- subset(ss.df, n == n_val)
  
  for (var in ss.vars) {
    counts <- table(ss.sub[[var]])
    barplot(counts,
            main = paste0(title.map[[var]], "\n", "n = ", n_val),
            xlab = xlab.map[[var]],
            col = "lightblue",
            border = "white",
            ylab = "",
            cex.lab = 1.3)
  }
}
dev.off()


# -----------------DGP2--------------------------
load("res_DGP2.rdata")
ss.df <- est.final %>% filter(method=="SIO") %>%
  group_by(method,rep,n) %>%
  summarise(ss.01=mean(ss.01),
            ss.10=mean(ss.10),
            ss.00=mean(ss.00),
            ss.gamma=mean(ss.gamma),
            ss.eta=mean(ss.eta),
            ss.Y11=mean(ss.Y11),
            ss.Y00=mean(ss.Y00),.groups = "drop")


ss.vars <- c("ss.10", "ss.01", "ss.00", "ss.gamma", "ss.eta", "ss.Y11", "ss.Y00")

xlab.map <- list(
  ss.01 = expression(J[n]/2),
  ss.10 = expression(J[n]/2),
  ss.00 = expression(J[n]/2),
  ss.gamma = expression(k[n]),
  ss.eta = expression(s[n]),
  ss.Y11 = expression(s[n]),
  ss.Y00 = expression(s[n])
)

title.map <- list(
  ss.10 = "r=(1,0,1)",
  ss.01 = "r=(0,1,1)",
  ss.00 = "r=(0,0,1)",
  ss.gamma = "gamma",
  ss.eta = "eta",
  ss.Y11 = "mu1",
  ss.Y00 = "mu0"
)


pdf("plt_size_DGP2.pdf", height=6, width=15,pointsize = 17)
par(mfrow = c(2, 7), mar = c(5, 3, 2.8, 1)) 

for (n_val in c(1000, 2000)) {
  ss.sub <- subset(ss.df, n == n_val)
  
  for (var in ss.vars) {
    counts <- table(ss.sub[[var]])
    barplot(counts,
            main = paste0(title.map[[var]], "\n", "n = ", n_val),
            xlab = xlab.map[[var]],
            col = "lightblue",
            border = "white",
            ylab = "",
            cex.lab = 1.3)
  }
}
dev.off()

