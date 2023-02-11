set.seed(4321)

library(MASS)

library(Matrix)

library(truncnorm)


load("Regen_func.RData")

#######################################################################################
####################################################################################### 

#Function of success probabilities for Elliptical D*
etaj <- function(betaj0, zj0, betaj2, zj2){
  ret <- 0
  t.starj <- t(zj0 - z.star)%*%x
  
  if(prod((betaj2 < d & betaj2 > c)) == 1)
  {
    ret <- sum(t.starj*(c* (t.starj > 0) + d * (t.starj < 0))) - 0.5*t(zj0)%*%x%*%inv.xtx%*%t(x)%*%zj0 + 0.5*t(z.star)%*%x%*%inv.xtx%*%t(x)%*%z.star+ log(density_mv(betaj2,inv.xtx.tx%*%z.star,inv.xtx)) + sum(log(density_tn(zj2, betaj2)))
    ret <- ret - log(density_mv(betaj2,inv.xtx.tx%*%zj0,inv.xtx) + density_mv(betaj2,inv.xtx.tx%*%zj0,inv.xtx)*prod(density_tn(zj2,betaj2)) +  density_mv(betaj2,inv.xtx.tx%*%zj2,inv.xtx)*prod(density_tn(zj2,betaj0)) + prod(density_tn(zj2,betaj0)))
    ret <- exp(ret)
  }else{
    ret <- 0
  }
  return(ret)
}


for(j in 1:(nsim-2))
{
  betaj2 <- ac.out.1e6$beta[j+2,]
  zj2 <- ac.out.1e6$z[j+2,]
  betaj0 <- ac.out.1e6$beta[j,]
  zj0 <- ac.out.1e6$z[j,]
  eta[j] <- etaj(betaj0, zj0, betaj2, zj2)
}

# Number of times chain is in the small set
sum(eta > 0)
pos_eta = eta[(eta > 0) == TRUE]
length(pos_eta)

#Number of total regenerations
sum=0
for (i in 1:length(eta)) {
  sum = sum + rbinom(1,1,eta[i])
}
sum


##################################################################################
##################################################################################

#Function of success probabilities for Elliptical D*
etaj <- function(betaj0, zj0, betaj2, zj2){
  ret <- 0
  t.starj <- t(zj0 - z.star)%*%x
  
  if(t(betaj2 - beta.bar)%*%solve(COVARIANCE)%*%(betaj2 - beta.bar) < 10)
  {
    ret <- sum(t.starj*(c* (t.starj > 0) + d * (t.starj < 0))) - 0.5*t(zj0)%*%x%*%inv.xtx%*%t(x)%*%zj0 + 0.5*t(z.star)%*%x%*%inv.xtx%*%t(x)%*%z.star+ log(density_mv(betaj2,inv.xtx.tx%*%z.star,inv.xtx)) + sum(log(density_tn(zj2, betaj2)))
    ret <- ret - log(density_mv(betaj2,inv.xtx.tx%*%zj0,inv.xtx) + density_mv(betaj2,inv.xtx.tx%*%zj0,inv.xtx)*prod(density_tn(zj2,betaj2)) +  density_mv(betaj2,inv.xtx.tx%*%zj2,inv.xtx)*prod(density_tn(zj2,betaj0)) + prod(density_tn(zj2,betaj0)))
    ret <- exp(ret)
  }else{
    ret <- 0
  }
  return(ret)
}


for(j in 1:(nsim-2))
{
  betaj2 <- ac.out.1e6$beta[j+2,]
  zj2 <- ac.out.1e6$z[j+2,]
  betaj0 <- ac.out.1e6$beta[j,]
  zj0 <- ac.out.1e6$z[j,]
  eta[j] <- etaj(betaj0, zj0, betaj2, zj2)
}

# Number of times chain is in the small set
sum(eta > 0)
pos_eta = eta[(eta > 0) == TRUE]
length(pos_eta)

#Number of total regenerations
sum=0
for (i in 1:length(eta)) {
  sum = sum + rbinom(1,1,eta[i])
}
sum




















