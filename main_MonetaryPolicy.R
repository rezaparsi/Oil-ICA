# main_MonetaryPolicy <- function(path, bootstrapping = FALSE, boot_irf = FALSE,
#                        subsamples = FALSE) {

  # This code reproduces the results of Chapter 5 of the OBES paper
  #    'Causal Inference by Independent Component Analysis:
  #                Theory and Applications'
  # by A. Moneta, D. Entner, P.O. Hoyer, and A. Coad

  # path: path to the data file including file name, f.ex. "myPath/theFile.txt"
  # bootstrapping: boolean, if TRUE bootstrap coefficients of VAR and SVAR to
  #    determine significant coefficients (Tables 3 and 4)
  # boot_irf: boolean, if TRUE bootstrap impulse response functions to get
  #    confidence bands (Figure 7)
  # subsamples: boolean, if TRUE do analysis on subsamples (Figure 8)
  require(tseries)
  require(urca)
  
  set.seed(1)

  # ------ Load Data ------

  # The data are available at Ilian Mihov's homepage:
  # http://www.insead.edu/facultyresearch/faculty/personal/imihov/documents/mmp.zip
  # preprocess data as explainend in the paper
  # columns: GDP, PGDP, NBR, BR, FFR, PSCCOM
  path = paste(getwd(),"/AggregateDataH.txt",sep = "")
  Data <- read.table(path, header=FALSE)
  
  # replace Total reserves TR with borrowed reserves BR
  #Data[,4] <- Data[,4] - Data[,3] # TR - NBR
  colnames(Data)=c("Oil Prod.%","Log Oil Price","K Ind.")

  
  # ------ Lag Selection with different criteria------

  print('Lag Selection')
  SC <- matrix(nrow=15,ncol=3)
  for (i in 3:15){
    # put Data in canonical form with i time lags
    temp <- tsdata2canonicalform(Data,nlags=i)
    res <- VAR_estim(temp,"vecm",regstats=FALSE,corank=3)
    SC[i,] <- lag_selection(res$residuals,i)
  }
  lag_ak <- which.min(SC[,1])
  lag_hq <- which.min(SC[,2])
  lag_sc <- which.min(SC[,3])
  cat('Akaike: ', lag_ak, '\n')         # 7 - AKAIKE criterion
  cat('Hannan-Quinn: ', lag_hq, '\n')   # 4 - Hannan-Quinn criterion
  cat('Schwarz: ', lag_sc, '\n')        # 3 - Schwarz criterion

  nlags <- 24
  Data_can <- tsdata2canonicalform(Data,nlags)


  # ------ Estimate VAR-LiNGAM ------
  result <- VARLiNGAM(Data_can,"ols", pruning=FALSE, corank=0, ntests = FALSE,regstats=FALSE)

  nodelabels <- list()
  nodelabels[[1]] <- "Oil Prod.%"
  nodelabels[[2]] <- "LR Oil Price"
  nodelabels[[3]] <- "K Ind."

  B0hat <- result$Bhat[[1]]
  ord <- result$var_order # causal order of variables

  print("causal variable order: ")
  print(nodelabels[ord])

  ## reduced-form residuals (from VAR)
  u_res<-result$resid
  cor(u_res)

  # just another way (fct) to get the instantaneous effect matrix, by using the
  # causal order from LiNGAM; B0ch = B0hat if pruning==FALSE in VARLiNGAM
  B0ch <- choleski(u_res,ord)
  
  
  ################
  #B0hat = B0ch

  ###########3

  ## structural residuals (from SVAR)
  Gamma0 <- diag(ncol(u_res))-B0hat # = P^(-1)
  v_res <- t(Gamma0 %*% t(u_res))

  print('Analysis of structural residuals: histograms, qq-plots')
  Gauss_Stats(v_res) # qq-plot, histogram, kurtosis
  print(Gauss_Tests(v_res)) # some normality tests

  print("correlation matrix of the structural resiuduals: ")
  print(cor(v_res)) 
  # this should be about diagonal, necessary condition for independence
  # i.e. structural residuals uncorrelated


  # ------ Bootstrap ------

  if (!bootstrapping) {

    # print coefficients
    cat('\ncoefficients of A1\n')
    print(round(result$Mhat[[1]],4))
    cat('\ncoefficients of A2\n')
    print(round(result$Mhat[[2]],4))
    cat('\ncoefficients of B0\n')
    print(round(result$Bhat[[1]],4))
    cat('\ncoefficients of B1\n')
    print(round(result$Bhat[[2]],4))
    cat('\ncoefficients of B2\n')
    print(round(result$Bhat[[3]],4))

  }

  if (bootstrapping) {

    ## bootstrap to prune coefficients
    print("Start bootstrap to get standard errors of coefficients.")

    # notation in paper: B = B_0, Gamma_i = B_i, A_i = M_i, i>0
    # sample from reduced-form residuals u_res and create new data

    cons<-result$const  # constant
    AA <- result$Mhat # reduced VAR estimates: A_1,...,A_p
    Bhat<-result$Bhat # structural VAR estimates: B_0,...,B_p
    nboot <- 100 # number of bootstrap samples

    # this prints coefficients and standard deviations
    boot_sd2(Data, cons, AA, Bhat, u_res, ord, nlags, nboot)

  } # end bootstrapping


  # ------ Impulse response functions with bootstrap ------

  ## calculate Impulse Response Functions
  print("Caluculate Impulse Response Functions.")  
  AA <- result$Mhat
  lag_irf <- 16
  IRF <- irf(AA, B0hat, lag_irf, u_res)

  # new window for plot of IRFs
  x11()
  par(mfrow=c(3,3))
  nd <- colnames(Data) # variable names
  n <- 0

  ## if not bootstrap confidence intervals of IRFs then plot IRFs without CI
  if (!boot_irf) {    

    for (i in 1:6){
      for (j in 1:6){

	n <- n+1   

	if (is.element(i,c(1,2,5)) && is.element(j,c(3,4,5)) ) {
	  di <- abs(max(IRF[n,]) - min(IRF[n,]))
	  mini <- min(IRF[n,]) - di/2 
	  maxi <- max(IRF[n,]) + di/2

	  matplot(0:25, IRF[n,], t="l", ylim=c(mini,maxi), xlab="", ylab="")
	  abline(h=0, col="grey", lwd=0.5)
	  title(paste("Responses of", nd[i], "to", nd[j]) )
	} # end if

      }
    }

  } # end !boot_irf

  if (boot_irf) {

    ## bootstrap confidence bands of IRFs and plot IRFs with CI
    print("Start bootstrap to get confidence bands of IRFs.")

    nit <- 100 #1000 # number of iteration in the bootstrap
    CI <- boot_conf_irf2(Data, AA, result$const, lag_irf, u_res, nit, ord)

    for (i in 1:3){
      for (j in 1:3){

	n <- n+1   

	if (is.element(i,c(1,2,3,4)) && is.element(j,c(1,2,3,4)) ) {
	  cirf <- c(IRF[n,],CI[[1]][n,],CI[[2]][n,])
	  di <- abs(max(cirf) - min(cirf))
	  mini <- min(cirf) - di/2 
	  maxi <- max(cirf) + di/2

	  matplot(0:16, IRF[n,], t="l", ylim=c(mini,maxi), xlab="", ylab="",xlim=c(0.6,16))
	  matplot(0:16, CI[[1]][n,], t="l", lty=2, add=TRUE,xlim=c(0.6,16))
	  matplot(0:16, CI[[2]][n,], t="l", lty=2, add=TRUE,xlim=c(0.6,16))
	  abline(h=0, col="grey", lwd=0.5)
	  title(paste("Responses of", nd[i], "to", nd[j]) )
	} # end if

      }
    }

 } # end boot_irf


  # ------ Robustness Analysis ------

#   if (subsamples) {
# 
#     set.seed(1)
#     macro_subsamples(Data, "vecm")
# #     macro_subsamples(Data, "ols")
# #     macro_subsamples(Data, "lad")
# #     macro_subsamples(Data, "lad", TRUE) # fm-LAD
# 
# 
#   }


# }
