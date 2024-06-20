#index function
g0f<-function(u) 5*u
g1f<-function(u) 10*exp(5*u)/(1+exp(5*u))


#varying coefficient function
beta21f<-function(u)  -0.3*(2-u)^3 + 1 
beta22f<-function(u)  (sin(pi*u))^2 + cos(1.5*pi*u)
beta11f<-function(u)  -0.2 + u^2  
beta12f<-function(u)  cos(0.5*pi*u) 


vadata1<-function(n){
  dx = 2
  dz = 2
  sZ = matrix(0, 2*n, dz)
  sX =  matrix(0, 2*n, dx-1)
  sZ[1,] = rnorm(dz, 0, 1)
  sX[1,] = rnorm(dx-1, 0, 1)
  Amat = matrix(0.1, (dz+dx-1), (dz+dx-1))
  diag(Amat) = rep(0.6, (dz+dx-1))
  for(t in 2:(2*n)){
    tvec = as.vector(Amat %*% c(sZ[(t-1),],sX[(t-1),])) + rnorm((dz+dx-1), 0, 1)
    sZ[t,] = tvec[1:dz]
    sX[t,] = tvec[-(1:dz)]
  }
  Z = apply(sZ[-(1:n),], 2, pnorm)
  X = cbind(rep(1,n), sX[-(1:n),1])
  
  U = seq_along(1:n) / n
  Beta1 = cbind(beta11f(U),beta12f(U))
  norm1 = sqrt(sum(colMeans(Beta1^2)))
  Beta1 = Beta1 / norm1
  
  Beta2 = cbind(beta21f(U),beta22f(U))
  norm2 = sqrt(sum(colMeans(Beta2^2)))
  Beta2 = Beta2 / norm2
  
  V1 = rowSums(Beta1 * Z)
  V2 = rowSums(Beta2 * Z)
  g0 = g0f(V1) - mean(g0f(V1))
  g1 = g1f(V2) - mean(g1f(V2))
  g = cbind(g0, g1)
  
  meany = rowSums(g * X)
  eps = rnorm(n, 0, 1)
  
  Y = meany + eps
  Beta = cbind(Beta1, Beta2)
  W = cbind(V1, V2)
  
  list(Y=Y, Z=Z, X=X, U=U, Beta=Beta, g=g, W=W, meany = meany, eps = eps)
}

