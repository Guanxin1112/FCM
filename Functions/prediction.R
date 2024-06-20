########################## MSPE of our model ##################################
Sp_our<-function(sample_test, sample_train, q, k, k1, y, z, x, u,h1, h2){
  dx = ncol(x); dz = ncol(z)
  basis = create.bspline.basis(rangeval=c(0,1), nbasis =(q+k) , norder=q)
  bpu = eval.basis(u, basis)
  L = Spest_mult1(q, k=k, k1=k1, y=y[-sample_test], z=z[-sample_test,], x=x[-sample_test,], u=u[-sample_test])
  LKer = twostep3(L$Betanew, L$gnew, L$gnewdev, h1, h2, y=y[-sample_test], z=z[-sample_test,], x=x[-sample_test,], u=u[-sample_test])
  delta_train = L$deltanew
  lambda_train = L$lambdanew
  delta_train1 = matrix(delta_train , ncol = dx, byrow= F)
  lambda_train1 = matrix(lambda_train, ncol=dx, byrow = F)
  
  ###spline
  Betapre = NULL
  gnew = matrix(0, length(sample_test), dx)
  for(l in 1:dx){
    btemp = NULL
    for(j in 1:dz){
      btemp = cbind(btemp,as.vector(bpu[sample_test,] %*% delta_train1[((j-1)*(q+k)+1):(j*(q+k)),l]))
    }
    Betapre = cbind(Betapre, btemp)
    wtemp = rowSums(btemp * z[sample_test,])
    wbasis = create.bspline.basis(rangeval=c(min(wtemp),max(wtemp)), nbasis =(q+k1) , norder=q)
    nbpw  = eval.basis(wtemp, wbasis)
    gnew[,l] = as.vector(nbpw %*% lambda_train1[,l]) 
  }

  #####kernel
  gpre_ker = gprel = NULL
  for(l in 1:dx){
    if(dx > 2){ yl = y[-sample_test] - rowSums(L$gnew[,-l] * x[-sample_test,-l])}
    if(dx == 2){ yl = y[-sample_test] - L$gnew[,-l] * x[-sample_test,-l]}
    if(dx == 1){ yl = y[-sample_test]}
    xl = x[-sample_test,l]
    etaoldl =  rowSums(L$Betanew[,((l-1)*dz+1):(l*dz)] * z[-sample_test,])
    wnew = rowSums(Betapre[,((l-1)*dz+1):(l*dz)] * z[sample_test,])
    for(j in 1:length(sample_test)){
      wdiff = (etaoldl - wnew[j]) / h1[l]
      Cmat = cbind(xl, xl*wdiff)
      Mmat = diag(dnorm(wdiff)/ h1[l],length(sample_train),length(sample_train))
      gprel[j] = c(1,0) %*% ginv( t(Cmat) %*% Mmat %*% Cmat ) %*% t(Cmat) %*% Mmat %*% yl
    }
    gpre_ker = cbind(gpre_ker, gprel)
  }
  ypre_ker= rowSums(gpre_ker * x[sample_test,])
  
  y_test <-  y[sample_test]
  MSPE_ker = mean((y_test - ypre_ker)^2) 
  
  list(ypre_ker = ypre_ker, MSPE_ker = MSPE_ker )
  
}

