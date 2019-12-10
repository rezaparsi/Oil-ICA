#  *****************************************************************
#  Copyright (c) Erik G. Learned-Miller, 2004.
#  *****************************************************************
#  RADICAL   Solve the ICA problem in arbitrary dimension.
#  
#     R translation: B. M. Kelm, 2006
# 
#     Version 1.1. Major bug fix. Faster entropy estimator.
#  
#     Apr.1, 2004. Major bug fix. Whitening matrix was wrong. Thanks
#       to Sergey Astakhov for spotting this one.
# 
#     Mar.28, 2004. Sped up inner loop by about 20% with a better
#       entropy estimator.
# 
#     Version 1.0. First release.
#   
#     RADICAL(X, n.comp=dim(X)[2], K=150, AUG_FLAG=1, reps=30, sweeps = n.comp-1 , stdev=0.175)
# 
#     X:        the matrix of mixed components, with one component per
#               column.
#
#     n.comp:   The number of independent components to construct.
#  
#     K:        The number of angles at which to evaluate the contrast
#               function. The ICA contrast function will be evaluated
#               at K evenly spaced rotations from -Pi/4 to Pi/4. For
#               small data sets (less than a few hundred points), the
#               default value of 150 should work well. For larger data
#               sets, very small benefits may be realized by
#               increasing the value of K, but keep in mind that the
#               algorithm is linear in K.
# 
#     AUG_FLAG: This flag is set to 1 by default, which indicates
#               that the data set will be "augmented" as discussed
#               in the paper. If this flag is set to 0, then the
#               data will not be augmented, and the next two 
#               arguments are meaningless. For large data
#               sets with more than 10000 points, this flag should
#               usually be set to 0, as there is usually no need to
#               augment the data in these cases.
# 
#     reps:     This is the number of replicated points for each  
#               original point. The default value is 30. The larger
#               the number of points in the data set, the smaller
#               this value can be. For data sets of 10,000 points or
#               more, point replication should be de-activated by setting
#               AUG_FLAG to 0 (see above).
# 
#     stdev:    This is the standard deviation of the replicated points. I
#               can't give too much guidance as to how to set this
#               parameter, but values much larger or smaller than
#               the default don't seem to work very well in the
#               experiments I've done. 
#
#
#     OUTPUT:
#
#     Wopt:     the "unmixing matrix" Wopt. 
#
#     Yopt:     Wopt applied to the mixed components X produces the 
#               approximately independent components Yopt, with one 
#               component per column.
#     pc:       the principal components (result from "princomp")
#
#     center:   the data mean
#

radical <- function(X, n.comp=dim(X)[2], K=150, AUG_FLAG=1, reps=30, sweeps = n.comp-1 , stdev=0.175) {
  
  vasicekm <- function(v,m) {
    len  <- length(v);
    vals <- sort(v);
    
    # Note that the intervals overlap for this estimator.
    intvals <- vals[(m+1):len]-vals[1:(len-m)];
    hvec    <- log(intvals);
    
    return(sum(hvec))	
  }
  
  radicalOptTheta <- function(x,stdev,m,reps,K,range) {
    # m is the number of intervals in an m-spacing
    # reps is the number of points used in smoothing
    # K is the number of angles theta to check for each Jacobi rotation.
    d <- dim(x)[2]
    N <- dim(x)[1]
    
    # This routine assumes that it gets whitened data.
    # First, we augment the points with reps near copies of each point.
    if (reps==1) {
      xAug <- x;
    } else {#     If the input data X is 1000x5, for example, then Yopt should
      #     also be 1000x5, and Wopt will be 5x5.
      #     ************************************************************* 
      
      # xAug = randn(d,N*reps)*stdev+repmat(x,[1,reps]);
      xAug <- t(matrix(rnorm(d*N*reps,sd=stdev),d,N*reps)+as.vector(t(x)))
    }
    
    # Then we rotate this data to various angles, evaluate the sum of 
    # the marginals, and take the min.
    perc <- range/(pi/2);
    numberK <- perc*K;
    start <- floor(K/2-numberK/2)+1;
    endPt <- ceiling(K/2+numberK/2);
    
    ent <- rep(0, K)
    for (i in 1:K) {
      # Map theta from -pi/4 to pi/4 instead of 0 to pi/2.
      # This will allow us to use Amari-distance for test of
      # convergence.
      theta <- (i-1)/(K-1)*pi/2-pi/4;
      # rot=[cos(theta) -sin(theta); sin(theta) cos(theta)];
      rot <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
      rotPts <- xAug %*% t(rot);
      
      marginalAtTheta <- rep(0,d)
      for (j in 1:d) {
        marginalAtTheta[j]=vasicekm(rotPts[,j],m);
      }
      ent[i] <- sum(marginalAtTheta);
    }
    
    entSort   <- sort(ent, index.return = T);
    thetaStar <- (entSort$ix[1]-1)/(K-1)*pi/2-pi/4;
    # rotStar=[cos(thetaStar) -sin(thetaStar); sin(thetaStar) cos(thetaStar)];
    rotStar   <- matrix(c(cos(thetaStar), sin(thetaStar), -sin(thetaStar), cos(thetaStar)), 2, 2);
    
    cat(sprintf("rotated %5.2f degrees.\n",thetaStar/(2*pi)*360));
    
    return(list(thetaStar=thetaStar,rotStar=rotStar))
  }
  
  X   <- as.matrix(X)
  N   <- dim(X)[1]            # note that X is N x dim (different from the Matlab version!)
  dim <- n.comp;
  m   <- floor(sqrt(N));      # m for use in m-spacing estimator.
  
  if (AUG_FLAG==0) {reps=1}
  sweeps <- as.integer(sweeps)
  
  # ****************
  # Whiten the data. Store the whitening operation to combine with
  # rotation matrix for total solution.
  pcX     <- princomp(X)
  W       <- sweep(pcX$loadings[,1:n.comp], 2, pcX$sdev[1:n.comp], FUN="/"); # Whitening_mat
  Xwhite  <- sweep(X, 2, pcX$center) %*% W
  
  oldTotalRot  <- diag(rep(1, dim))
  totalRot     <- oldTotalRot
  sweepIter    <- 0                      # Current sweep number.
  xcur         <- Xwhite 
  
  # K represents the number of rotations to examine on the FINAL
  # sweep. To optimize performance, we start with a smaller number of
  # rotations to examine. Then, we increase the
  # number of angles to get better resolution as we get closer to the
  # solution. For the first half of the sweeps, we use a constant
  # number for K. Then we increase it exponentially toward the finish.
  finalK <- K;
  startKfloat <- (finalK/1.3^(ceiling(sweeps/2)));
  newKfloat <- startKfloat;
  
  for (sweepNum in 1:sweeps) {
    cat(sprintf("Sweep # %d of %d.\n",sweepNum,sweeps));
    range <- pi/2;
    # % Compute number of angle samples for this sweep.
    
    if (sweepNum>(sweeps/2)) {
      newKfloat <- newKfloat*1.3;
      newK <- floor(newKfloat);
    } else {
      newKfloat <- startKfloat;
      newK <- max(30,floor(newKfloat)); 
    }
    # *********************************************************
    # Iterate over all possible Jacobi rotations.
    # *********************************************************
    for (i in 1:(dim-1)) {
      for (j in (i+1):dim) {
        cat(sprintf("Unmixing dimensions %02d and %02d ...",i,j));
        # **********************************************
        # Extract dimensions (i,j) from the current data.
        # **********************************************
        curSubSpace <- cbind(xcur[,i], xcur[,j]);
        
        # ***************************************************
        # Find the best angle theta for this Jacobi rotation.
        # ***************************************************
        ROT <- radicalOptTheta(curSubSpace,stdev,m,reps,newK,range);
        # ROT: thetaStar,rotStar
        
        # *****************************************
        # Incorporate Jacobi rotation into solution.
        # *****************************************
        newRotComponent      <- diag(rep(1, dim));
        newRotComponent[i,i] <- cos(ROT$thetaStar);
        newRotComponent[i,j] <- -sin(ROT$thetaStar);
        newRotComponent[j,i] <- sin(ROT$thetaStar);
        newRotComponent[j,j] <- cos(ROT$thetaStar);
        totalRot <- newRotComponent %*% totalRot;
        xcur <- Xwhite %*% t(totalRot);
      }
    }
    
    oldTotalRot <- totalRot;
  }
  
  Wopt <- W %*% t(totalRot)
  
  return( list( Wopt = Wopt,
                Yopt = X %*% Wopt, 
                pc = pcX, 
                center = pcX$center ))
}

#if (F) {
## DEMO 2d

cat("Loading 2D example...\n")

path <- "./examples/";
data_2d_ind    <- t(as.matrix(read.table(file=paste(path, "data_2d_ind", sep=""))));
data_2d_mixed  <- t(as.matrix(read.table(file=paste(path, "data_2d_mixed", sep=""))));
A_2d           <- t(as.matrix(read.table(file=paste(path, "A_2d", sep=""))));

cat("Starting RADICAL...\n")
res <- radical(data_2d_mixed, AUG_FLAG=0)
data_2d_unmixed <- res$Yopt
colnames(data_2d_unmixed) <- colnames(data_2d_ind)

X11(height=4, width=12);
par(mfcol=c(1, 3), mai=c(.5,.5,.5,.5));
plot(data_2d_ind, main="Original independent components.")
plot(data_2d_mixed, main="Mixed components.")
plot(data_2d_unmixed, main="Unmixed components.");
#}