boot_conf_irf2 <- function(YY, AA, cons, lag_irf, ures, niter, ord){
  
  # get confidence intervals for the irf (impulse response function) using
  # bootstrapping, for VECM model only, cointegration rank 3 used
  #
  # INPUT
  # YY: original Data
  # AA: list of estimated coefficients of reduced VAR
  #     Ahat[[1]] = A_1, ... Ahat[[p]] = A_p
  # cons: constant term from reduced VAR estimation
  # lag_irf: number of lags included in the calculation of the irf
  # ures: residuals from reduced form VAR
  # niter: number of bootstrap runs
  # ord: causal order of variables in original data set
  # 
  # OUTPUT
  # a list with two elements:
  #   $[[1]] lower bound of confidence interval
  #   $[[2]] upper bound of confidence interval
  #     each of which contains a k^2 x lag_irf matrix, with k=#variables, where
  #     each row contains the upper/lower bound of the confidence interval of
  #     the response of one variable on another varialbe at all time lags
  #     (0,...,lag_irf):
  #     rows 1...k contain the bounds of the 1st var to all other vars
  #     rows k+1...2k contain the bounds of the 2nd var to all other vars
  #     etc.
  
  t <- nrow(ures) # number of samples
  k <- ncol(ures) # number of variables
  p <- length(AA) # number of lags
  
  MM <- array(0, c(k*k, lag_irf+1, niter))
  realord = list()
  
  # start bootstrap
  it=1
  while (it <= niter) {
    
    if (it%%10==0) {
      cat("bootstrap run", it, "out of", niter, "\n")
    }
    
    ll <- 100
    lb <- t%%ll  # 7
    mc <- t%/%ll # (t-lb)/ll
    pp <- sample(mc,mc, replace=TRUE)*ll
    ind <- 1:(t-lb)
    
    for(w in 1:mc) {
      ind[(w*ll - ll +1):(w*ll)] <- (pp[w]-ll+1): pp[w]
    }
    if (lb!=0) {
      s1 <- 1:lb
      ind <- c(s1, ind + lb)
    }
    
    # generate new data
    unew <- ures[ind,]
    Ynew <- matrix(0,nrow=(t+p),ncol=k)
    Ynew[1:p,] <- as.matrix(YY[1:p,])
    for(i in (p+1):(t+p)) {
      for(j in 1:p) {
        Y <- AA[[j]]%*%t(YY[i-j,])
        Ynew[i,] <- Ynew[i,]+Y 
      }
      Ynew[i,] <- cons + Ynew[i,] + unew[i-p,]
    }
    Ynew <- as.data.frame(Ynew)
    
    res_new <- VARLiNGAM(tsdata2canonicalform(Ynew,nlags), est_meth = "vecm", ntests = FALSE, regstats=FALSE, corank=3)
    AA_new <- res_new$Mhat 
    unewh <- res_new$resid
    realord[[it]] = res_new$var_order
    B0_new <- choleski(unewh,ord)
    MM[,,it] <- irf(AA_new, B0_new, lag_irf, unewh) # irf on bootstrap data
    if (sum(realord[[it]]==ord)<3){it <- it-1}
    it = it+1
  }
  
  # save bootstrap results
  DOWN <- matrix(nrow=k*k, ncol=lag_irf+1) # lower bound
  UP <- matrix(nrow=k*k, ncol=lag_irf+1) # upper bound
  for (m in 1:(lag_irf+1)) {
    for (n in 1:(k*k)) {
      mg<-MM[n,m,]
      DOWN[n,m] <- quantile(mg,0.005) 
      UP[n,m] <- quantile(mg,0.995)
    }
  }
  ris <- as.list(1:2)
  ris[[1]] <- DOWN
  ris[[2]] <- UP
  ris
}
