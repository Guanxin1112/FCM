#-------Example 2----------
# library packages
library(MASS)
library(splines)
library(alabama)
library(fda)

# input functions 
filepathfun = '~/Functions/'
source(paste(filepathfun,'three_step.R',sep=''))
source(paste(filepathfun,'test.R',sep=''))


# sample size
n = 1000
nsim = 500
boots = 400
q = 3
k = floor(n^(1/(2*q+1))) 
k1 = floor(n^(1/(2*q+1))) 
h1 = c(0.3*(2000/n)^(-1/5),0.2*(2000/n)^(-1/5));
h2 = c(0.12*(2000/n)^(-1/5),0.12*(2000/n)^(-1/5))
theta = 0
alpha = 0.05


PC = rep(NA,nsim)
for(ns in 1:nsim){
  Dat = vadata2(n, theta)
  Y=Dat$Y;X=Dat$X;Z=Dat$Z;Beta=Dat$Beta;U=Dat$U;g=Dat$g
  L = Spest_mult1(q, k, k1, y=Y, z=Z, x=X, u=U)
  LKer = twostep31(L$Betanew, L$gnew, L$gnewdev, h1, h2, y=Y, z=Z, x=X, u=U)
  est_Beta11 = LKer$Betanewlp[,1]
  est_Beta12 = LKer$Betanewlp[,2] 
  beta2 = rep(mean(est_Beta12),n)
  Betanewlp = cbind(est_Beta11, beta2)
  #wild bootstrap method can obtain SD on each sample point Ui
  Bootest_Beta12 = matrix(0, nrow=n, ncol=boots/2)
  for(b in 1:(boots/2)){
    LKerB = twostep32(Betanewlp, LKer$gnewlp, LKer$g1newlp, LKer$reslp, h1, h2, z=Z, x=X, u=U)
    Bootest_Beta12[,b] = LKerB$Betanewlp[,2]
  }
  sd_Beta12 = apply(Bootest_Beta12, 1, sd)
  #
  T_Beta12 = rep(0, boots/2)
  for(b in 1:(boots/2)){
    LKerM = twostep32(Betanewlp, LKer$gnewlp, LKer$g1newlp, LKer$reslp, h1, h2, z=Z, x=X, u=U)
    T_Beta12[b] = max(abs(LKerM$Betanewlp[,2] - beta2) * (1 / sd_Beta12))
  }
  
  s12 = sort(T_Beta12)
  C12 = s12[ceiling((boots/2)*(1-alpha))]
  
  PC[ns] = as.numeric(max(abs(est_Beta12 - beta2) * (1 / sd_Beta12)) > C12)
  
}

sum(PC) / nsim 
