#-------Example 1 for estimation----------
# library packages
library(MASS)
library(splines)
library(alabama)
library(fda)

# input functions 
filepathfun = '~/Functions/'
source(paste(filepathfun,'three_step.R',sep=''))
source(paste(filepathfun,'prediction.R',sep=''))
source(paste(filepathfun,'datasim.R',sep=''))


# sample size
n = 600
nsim = 500
q = 3
k = floor(n^(1/(2*q+1))) 
k1 = floor(n^(1/(2*q+1)))
h1 = c(0.3*(2000/n)^(-1/5),0.2*(2000/n)^(-1/5));
h2 = c(0.12*(2000/n)^(-1/5),0.12*(2000/n)^(-1/5))

# result output
mse_est_beta1f = mse_est_beta2f = mse_est_beta3f = mse_est_beta4f = rep(0, nsim)
mse_est_beta1fKer = mse_est_beta2fKer = mse_est_beta3fKer = mse_est_beta4fKer = rep(0, nsim)
mse_est_g1f = mse_est_g2f = rep(0, nsim)
mse_est_g1fKer = mse_est_g2fKer = rep(0, nsim)



for(ns in 1:nsim){
  Dat = vadata1(n)
  Y=Dat$Y;X=Dat$X;Z=Dat$Z;Beta=Dat$Beta;U=Dat$U;g=Dat$g
  #------one iteration
  L = Spest_mult1(q, k=k, k1=k1, y=Y, z=Z, x=X, u=U)
  mse_est_beta1f[ns] = sqrt(mean((Beta[,1] - L$Betanew[,1])^2))
  mse_est_beta2f[ns] = sqrt(mean((Beta[,2] - L$Betanew[,2])^2))
  mse_est_beta3f[ns] = sqrt(mean((Beta[,3] - L$Betanew[,3])^2))
  mse_est_beta4f[ns] = sqrt(mean((Beta[,4] - L$Betanew[,4])^2))
  mse_est_g1f[ns] = sqrt(mean((g[,1] - L$gnew[,1])^2))
  mse_est_g2f[ns] = sqrt(mean((g[,2] - L$gnew[,2])^2))
  LKer = twostep3(L$Betanew, L$gnew, L$gnewdev, h1, h2, y=Y, z=Z, x=X, u=U)
  mse_est_beta1fKer[ns] = sqrt(mean((Beta[,1] - LKer$Betanewlp[,1])^2))
  mse_est_beta2fKer[ns] = sqrt(mean((Beta[,2] - LKer$Betanewlp[,2])^2))
  mse_est_beta3fKer[ns] = sqrt(mean((Beta[,3] - LKer$Betanewlp[,3])^2))
  mse_est_beta4fKer[ns] = sqrt(mean((Beta[,4] - LKer$Betanewlp[,4])^2))
  mse_est_g1fKer[ns] = sqrt(mean((g[,1] - LKer$gnewlp[,1])^2))
  mse_est_g2fKer[ns] = sqrt(mean((g[,2] - LKer$gnewlp[,2])^2))
}

mean(mse_est_beta1f);sd(mse_est_beta1f)
mean(mse_est_beta2f);sd(mse_est_beta2f)
mean(mse_est_beta3f);sd(mse_est_beta3f)
mean(mse_est_g1f);sd(mse_est_g1f)
mean(mse_est_g2f);sd(mse_est_g2f)

mean(mse_est_beta1fKer);sd(mse_est_beta1fKer)
mean(mse_est_beta2fKer);sd(mse_est_beta2fKer)
mean(mse_est_beta3fKer);sd(mse_est_beta3fKer)
mean(mse_est_g1fKer);sd(mse_est_g1fKer)
mean(mse_est_g2fKer);sd(mse_est_g2fKer)

