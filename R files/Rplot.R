set.seed(4321)

library(MASS)

library(Matrix)

library(truncnorm)


load("Regen_func.RData")

pdf(file="Rplot.pdf")
#ACF of beta samples from AC
par(mfrow = c(1,3))
acf(ac.out.1e6$beta[,1], ylim = c(0.85,1), main = expression(paste("AC for ", beta(1))))
acf(ac.out.1e6$beta[,2], ylim = c(0.85,1), main = expression(paste("AC for ", beta(2))))
acf(ac.out.1e6$beta[,3], ylim = c(0.85,1), main = expression(paste("AC for ", beta(3))))
dev.off()