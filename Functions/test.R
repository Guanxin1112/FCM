## the identification condition is \|beta(u)\|_{L_2} =1
twostep31<-function(Betaold, gold, golddev, h1, h2, y, z, x, u){
  dz = ncol(z); dx = ncol(x); n = nrow(z)
  #step 1: compute a_j,b_j; j=1:n 
  gnewlp = NULL
  g1newlp = NULL
  Betanewlp = NULL
  wblocknew = NULL
  reslp = NULL
  for(l in 1:dx){
    if(dx > 2){ yl = y - rowSums(gold[,-l] * x[,-l])}
    if(dx == 2){ yl = y - gold[,-l] * x[,-l]}
    if(dx ==1 ){yl = y}
    xl = x[,l]
    etaoldl = rowSums(Betaold[,((l-1)*dz+1):(l*dz)] * z)
    gnewl = rep(0, n)
    g1newl = rep(0, n)
    for(j in 1:n){
      w0 = etaoldl[j]
      wdiff = (etaoldl - w0) / h1[l]
      wnear = dnorm(wdiff) / h1[l]
      z1 = xl
      z2 = wdiff * xl
      fit = lm(yl~z1+z2-1, weights = wnear)
      coeff = fit$coefficients
      gnewl[j] = coeff[1]
      g1newl[j] = coeff[2] / h1[l]
    }
    
    # step 2: compute Betanew
    Betanewl = matrix(0, n, dz)
    gdevcom = as.vector(golddev[,l] * xl)
    ystarl = yl - gold[,l] * xl + gdevcom * etaoldl
    for(j in 1:n){
      u0 = u[j]
      udiff = (u - u0) / h2[l]
      unear = dnorm(udiff) / h2[l]
      zz =  cbind(z, z * udiff)
      xstarl = zz *  gdevcom
      fitb = lm(ystarl~xstarl+0, weights = unear)
      bcoeff = fitb$coefficients
      Betanewl[j,] = bcoeff[1:dz]
    }
    
    # calculate \|beta(u)\|_L_2
    norml =  sqrt(sum(colMeans(Betanewl^2)))
    
    # identification -- divide \|beta(u)\|_L_2
    Betanewl = Betanewl / norml
    etanewl = rowSums(Betanewl * z)
    Betanewlp = cbind(Betanewlp, Betanewl)
    wblocknew = cbind(wblocknew, etaoldl)
    gnewlp = cbind(gnewlp, gnewl)
    g1newlp = cbind(g1newlp, g1newl)
    resl = yl - gnewl * xl
    reslp = cbind(reslp, resl)
  }
  
  
  reslp = reslp - matrix(apply(reslp,2,mean), nrow = n, ncol = dx, byrow = T)
  fit = rowSums(gnewlp * x)
  res = y - fit
  
  list(Betanewlp = Betanewlp, gnewlp = gnewlp, g1newlp = g1newlp, wblocknew = wblocknew, 
       fit = fit, res = res, reslp = reslp)
} 



##for wild bootstrap
twostep32<-function(Betaoldlp, goldlp, golddevlp, reslp, h1, h2, z, x, u){
  dz = ncol(z); dx = ncol(x); n = nrow(z)
  #step 1: compute a_j,b_j; j=1:n 
  gnewlp = NULL
  g1newlp = NULL
  Betanewlp = NULL
  for(l in 1:dx){
    booter = rnorm(n, 0, 1)
    xl = x[,l]
    yl = goldlp[,l] * xl + reslp[,l] * booter
    etaoldl = rowSums(Betaoldlp[,((l-1)*dz+1):(l*dz)] * z)
    # step 2: compute Betanew
    Betanewl = matrix(0, n, dz)
    gdevcom = as.vector(golddevlp[,l] * xl)
    ystarl = yl - goldlp[,l] * xl + gdevcom * etaoldl
    for(j in 1:n){
      u0 = u[j]
      udiff = (u - u0) / h2[l]
      unear = dnorm(udiff) / h2[l]
      zz =  cbind(z, z * udiff)
      xstarl = zz *  gdevcom
      fitb = lm(ystarl~xstarl+0, weights = unear)
      bcoeff = fitb$coefficients
      Betanewl[j,] = bcoeff[1:dz]
    }
    
    # calculate \|beta(u)\|_L_2
    norml =  sqrt(sum(colMeans(Betanewl^2)))
    # identification -- divide \|beta(u)\|_L_2
    Betanewl = Betanewl / norml
    Betanewlp = cbind(Betanewlp, Betanewl)
  }
  
  list(Betanewlp = Betanewlp)
} 



vadata2<-function(n, theta){
  beta11f<-function(u)  0.1*exp(-0.7+3.5*u)
  beta12f<-function(u)  2*sin(2*pi*(u-0.5)^2)
  g0f<-function(u) 5*u
  
  dx = 1; dz = 2
  sZ = matrix(0, 2*n, dz)
  sZ[1,] = rnorm(dz, 0, 1)
  Amat = matrix(0.05, dz, dz)
  diag(Amat) = rep(0.15, dz)
  for(t in 2:(2*n)){
    tvec = as.vector(Amat %*% sZ[(t-1),]) + rnorm((dz+dx-1), 0, 1)
    sZ[t,] = tvec[1:dz]
  }
  Z = apply(sZ[-(1:n),], 2, pnorm)
  X = matrix(rep(1,n), nrow = n, ncol = 1)
  
  U = seq_along(1:n) / n
  
  c = integrate(f = beta12f, 0, 1)$value
  beta12ft<-function(u) c + theta*beta12f(u)
  
  Beta1 = cbind(beta11f(U),beta12ft(U))
  norm1 = sqrt(sum(colMeans(Beta1^2)))
  Beta = Beta1 / norm1
  
  V1 = rowSums(Beta * Z)
  g = g0f(V1) - mean(g0f(V1))

  
  meany = rowSums(g * X)
  eps = rnorm(n, 0, 0.7)
  
  #var(meany)/var(eps)
  Y = meany + eps
  W = V1
  
  list(Y=Y, Z=Z, X=X, U=U, Beta=Beta, g=g, W=W, meany = meany, eps = eps)
}



