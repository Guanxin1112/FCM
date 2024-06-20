bfunini_mult1<-function(q, k, y, z, x, u){
  dz = ncol(z);dx = ncol(x);n = nrow(z)
  basis = create.bspline.basis(rangeval=c(0,1), nbasis =(q+k), norder=q)
  bpu = eval.basis(u, basis)
  ##Initialization step  estimation
  Phi = NULL
  rex = NULL
  for(l in 1:dz){ Phi = cbind(Phi, bpu * z[, l]) }
  for(l in 1:dx){ rex = cbind(rex, Phi * x[, l]) }
  rex = cbind(x, rex)
  atheta = as.vector(ginv(t(rex) %*% rex)%*% t(rex) %*% y) #ridge regression
  
  part_atheta = atheta[-(1:dx)]
  matrix_atheta = matrix(part_atheta, ncol=dx, byrow= F)
  
  #identification -- beta_1 nondecreasing, \|beta_(u)\| = 1
  for(l in 1:dx){
    matrix_atheta[1:(q+k),l] = sort(matrix_atheta[1:(q+k),l])
  }
  
  Betaold = NULL
  for(l in 1:dx){
    btemp = NULL
    for(j in 1:dz){
      ttt = as.vector(bpu %*% matrix_atheta[((j-1)*(q+k)+1):(j*(q+k)),l])
      btemp = cbind(btemp, ttt)
    }
    norml = sqrt(sum(colMeans(btemp^2)))
    btemp = btemp / norml
    matrix_atheta[,l] = matrix_atheta[,l] / norml
    Betaold = cbind(Betaold, btemp)
  }
  
  deltaold = as.vector(matrix_atheta)
  
  return(deltaold)
}



Est_nlfun_mult1<-function(deltaold, q, k, k1, y, z, x, u){
  dz = ncol(z); dx = ncol(x); n = nrow(z)
  basis = create.bspline.basis(rangeval=c(0,1), nbasis =(q+k) , norder=q)
  bpu = eval.basis(u, basis)
  Phi = NULL
  for(l in 1:dz){Phi = cbind(Phi, bpu * z[, l])}
  
  fn<-function(ddd){
    ddd1 = matrix(ddd, ncol = dx, byrow= F)
    D = NULL
    for(l in 1:dx){
      btemp = NULL
      for(j in 1:dz){
        btemp = cbind(btemp,as.vector(bpu %*% ddd1[((j-1)*(q+k)+1):(j*(q+k)),l]))
      }
      wtemp = rowSums(btemp * z)
      wbasis = create.bspline.basis(rangeval=c(min(wtemp),max(wtemp)), nbasis =(q+k1) , norder=q)
      bpwf  = eval.basis(wtemp, wbasis)
      D =  cbind(D, bpwf * x[,l])
    }
    
    lambdaold = as.vector(ginv(t(D) %*% D)%*% t(D) %*% y) 
    res = as.vector(y - D %*% lambdaold)
    # minimize objective function
    L1 = as.numeric(0.5*t(res) %*% (res))
    return(L1)
  }
  
  gn<-function(ddd){
    ddd1 = matrix(ddd, ncol = dx, byrow= F)
    D = NULL
    wblock = NULL
    for(l in 1:dx){
      btemp = NULL
      for(j in 1:dz){
        btemp = cbind(btemp,as.vector(bpu %*% ddd1[((j-1)*(q+k)+1):(j*(q+k)),l]))
      }
      wtemp = rowSums(btemp * z)
      wbasis = create.bspline.basis(rangeval=c(min(wtemp),max(wtemp)), nbasis =(q+k1) , norder=q)
      bpwf  = eval.basis(wtemp, wbasis)
      wblock = cbind(wblock, wtemp)
      D =  cbind(D, bpwf * x[,l])
    }
    lambdaold = as.vector(ginv(t(D) %*% D)%*% t(D) %*% y) 
    res = as.vector(y - D %*% lambdaold)
    
    
    eta = ginv(t(D) %*% D)%*% t(D) %*% Phi
    tildePhi =  Phi- D %*% eta
    lambdaold1 = matrix(lambdaold, ncol=dx, byrow = F)
    ghatdev = NULL
    for(l in 1:dx){
      wbasisl = create.bspline.basis(rangeval=c(min(wblock[,l]),max(wblock[,l])), nbasis =(q+k1) , norder=q)
      bpwdev  = eval.basis(wblock[,l], wbasisl, Lfdobj=1)
      ghatdev0 =  as.vector(bpwdev %*% lambdaold1[,l])
      ghatdevp = ghatdev0 * tildePhi
      ghatdev = cbind(ghatdev, ghatdevp * x[,l])
    }
    dLn1 = 0
    for(i in 1:n){dLn1 = dLn1 + res[i] * ghatdev[i,] }
    dLn1 = - dLn1
    dLn1
  }
  
  # R = optim(par=deltaold, fn=fn, gr = gn, method = "BFGS")
  Amat = matrix(0, nrow = dx*dz*(q+k), ncol = dx*dz*(q+k))
  for(l in 1:dx){
    for(i in (2+(l-1)*dz*(q+k)):(q+k+(l-1)*dz*(q+k))){Amat[i,(i-1)] = -1}
    for(i in (2+(l-1)*dz*(q+k)):(q+k+(l-1)*dz*(q+k))){Amat[i,i] = 1}
  }
  bvec = rep(0, dx*dz*(q+k))  - 1e-10
  # Amat %*% deltaold >= bvec
  R = constrOptim(theta=deltaold, f=fn, grad = gn, ui = Amat, ci = bvec, method = "BFGS")
  
  conver = R$convergence
  deltanew = R$par
  
  # standardized Beta
  deltanew1 = matrix(deltanew, ncol = dx, byrow= F)
  
  #identification -- beta_1 nondecreasing, \|beta_(u)\| = 1 
  nD = NULL
  wblocknew = NULL
  Betanew = NULL
  for(l in 1:dx){
    btemp = NULL
    for(j in 1:dz){
      btemp = cbind(btemp,as.vector(bpu %*% deltanew1[((j-1)*(q+k)+1):(j*(q+k)),l]))
    }
    norml = sqrt(sum(colMeans(btemp^2)))
    btemp = btemp / norml
    wtemp = rowSums(btemp * z)
    wbasis = create.bspline.basis(rangeval=c(min(wtemp),max(wtemp)), nbasis =(q+k1) , norder=q)
    nbpw  = eval.basis(wtemp, wbasis)
    nD =  cbind(nD, nbpw * x[,l])
    wblocknew = cbind(wblocknew, wtemp)
    Betanew = cbind(Betanew, btemp)
    deltanew1[,l] = deltanew1[,l] / norml
  }
  
  deltanew = as.vector(deltanew1)
  lambdanew = as.vector(ginv(t(nD) %*% nD)%*% t(nD) %*% y) 
  lambdanew1 = matrix(lambdanew, ncol=dx, byrow = F)
  gnew = matrix(0, n, dx)
  gnewdev = matrix(0, n, dx)
  for(l in 1:dx){
    wbasisl = create.bspline.basis(rangeval=c(min(wblocknew[,l]),max(wblocknew[,l])), nbasis =(q+k1) , norder=q)
    nbpw  = eval.basis(wblocknew[,l], wbasisl)
    gnew[,l] = as.vector(nbpw %*% lambdanew1[,l]) # - mean(as.vector(nbpw %*% lambdanew1[,l]))
    #
    nbpwdev  = eval.basis(wblocknew[,l], wbasisl, Lfdobj=1)
    gnewdev[,l] =  as.vector(nbpwdev %*% lambdanew1[,l])
  }
  
  list(Betanew = Betanew, deltanew=deltanew, gnew=gnew, lambdanew=lambdanew,
       gnewdev=gnewdev, wblocknew = wblocknew)
}



Spest_mult1<-function(q, k, k1, y, z, x, u){
  deltaini = bfunini_mult1(q, k, y, z, x, u)
  Res_step3 = Est_nlfun_mult1(deltaini, q, k, k1, y, z, x, u)
  Betanew = Res_step3$Betanew 
  deltanew = Res_step3$deltanew
  gnew = Res_step3$gnew
  lambdanew = Res_step3$lambdanew
  gnewdev = Res_step3$gnewdev
  wblocknew = Res_step3$wblocknew
  fit = rowSums(gnew * x)
  res = y - fit
  
  list(Betanew = Betanew,deltanew=deltanew, gnew=gnew,lambdanew=lambdanew,
       gnewdev=gnewdev,wblocknew = wblocknew, fit = fit, res = res)
}



## the identification condition is \|beta(u)\|_{L_2} =1
twostep3<-function(Betaold, gold, golddev, h1, h2, y, z, x, u){
  dz = ncol(z); dx = ncol(x); n = nrow(z)
  #step 1: compute a_j,b_j; j=1:n 
  gnewlp = NULL
  g1newlp = NULL
  Betanewlp = NULL
  wblocknew = NULL
  for(l in 1:dx){
    if(dx > 2){ yl = y - rowSums(gold[,-l] * x[,-l])}
    if(dx == 2){ yl = y - gold[,-l] * x[,-l]}
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
      # Beta1newl[j,] = bcoeff[-(1:dz)] / h2
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
  }
  
  fit = rowSums(gnewlp * x)
  res = y - fit
  
  list(Betanewlp = Betanewlp, gnewlp = gnewlp, g1newlp = g1newlp, wblocknew = wblocknew, 
       fit = fit, res = res)
} 




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
      # Beta1newl[j,] = bcoeff[-(1:dz)] / h2
    }
    
    # calculate \|beta(u)\|_L_2
    norml =  sqrt(sum(colMeans(Betanewl^2)))
    
    # identification -- divide \|beta(u)\|_L_2
    Betanewl = Betanewl / norml
    etanewl = rowSums(Betanewl * z)
    # gnewlpl = rep(0, n)
    # g1newlpl = rep(0, n)
    # for(j in 1:n){
    #   w0 = etanewl[j]
    #   wdiff = (etanewl - w0) / h1
    #   wnear = dnorm(wdiff) / h1
    #   Dl = cbind(xl, xl*wdiff)
    #   fit = lm(yl~0+Dl, weights = wnear)
    #   coeff = fit$coefficients
    #   gnewlpl[j] = coeff[1]
    #   g1newlpl[j] = coeff[2] / h1
    # }
    # Betanewlp = cbind(Betanewlp, Betanewl)
    # wblocknew = cbind(wblocknew, etanewl)
    # gnewlp = cbind(gnewlp, gnewlpl)
    # g1newlp = cbind(g1newlp, g1newlpl)
    
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

