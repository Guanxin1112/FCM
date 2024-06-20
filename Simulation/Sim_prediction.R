#-------Example 1 for MSPE----------
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
n0 = 20; m = 600
n = m+n0
q = 3
k = floor(n^(1/(2*q+1))); k1 = floor(n^(1/(2*q+1)))
h1 = c(0.3*(2000/n)^(-1/5),0.2*(2000/n)^(-1/5));
h2 = c(0.12*(2000/n)^(-1/5),0.12*(2000/n)^(-1/5))

nrep = 200
MSPE = rep(NA, nrep)
  
for(i in 1:nrep){
  #print(paste0("Start loop ", i, " :"))
  sample_test = sample((n*0.2):(n*0.8), size=n0, replace = FALSE)
  sample_train = seq(1, n, by=1)[-sample_test]
  
  data <- vadata1(n)
  nonpara = Sp_our(sample_test, sample_train, q, k, k1, y=data$Y, z=data$Z, x=data$X, u=data$U, h1, h2)
  MSPE[i] = nonpara$MSPE_ker
}
  
mean(MSPE)
