require(fastICA)

result = list()
resultchol = list()
order = list()
sim = 100
pinup = rep(0,sim)
set.seed(1)

for (counter in 1:sim){
  
  n = 100
  x = rep(0,n)
  y = rep(0,n) 
  z = rep(0,n)
  p = rep(0,n)


  
 name = "norm2"
 Random = cbind(rlnorm(n-1,0,0.5), rlnorm(n-1,0,0.5),rlnorm(n-1,0,0.5),rlnorm(n-1,0,0.5))


# Random[,2] = rep(0,n-1)
# 
# 
# 
# o = 10
# size = n*o/100
# outliar = sample.int(n-1,size)
# Random[outliar[1:size/2],2]= ((n-1)/size)^(0.5)
# Random[outliar[size/2 + 1:size],2]= -((n-1)/size)^(0.5)

# o = 10
# size = n*o/100
# outliar = sample.int(n-1,size)
# Random[outliar[1:size/2],3]= ((n-1)/size)^(0.5)
# Random[outliar[size/2 + 1:size],3]= -((n-1)/size)^(0.5)
# 
# o = 10
# size = n*o/100
# outliar = sample.int(n-1,size)
# Random[outliar[1:size/2],4]= ((n-1)/size)^(0.5)
# Random[outliar[size/2 + 1:size],4]= -((n-1)/size)^(0.5)

#Random = Random + cbind( 0* outliar * Random[,1],rep(0,n-1),rep(0,n-1),rep(0,n-1)) 

# lamda = 1
# a = -1
# b = 2
# x0 = 2* (((b-a)^2)/3 +((a+b)^2)/lamda)^(-0.5)
# y0 = -((a+b) * x0)/(2*lamda) 
# Random = cbind(3^0.5 *runif(n-1,-1,1),3^0.5 *runif(n-1,-1,1),3^0.5 *runif(n-1,-1,1),
#                x0 * runif(n-1,a,b) + y0 * rpois(n-1,lamda))


  
  for (i in 2:n){
    x[i] = 0.90 *x[i-1] + Random[i-1,1]
    y[i] = 0.3 *x[i] + 0.90 *y[i-1] + Random[i-1,2]
    z[i] = 0.3 *x[i] + 0.3 *y[i] + 0.90 *z[i-1] + Random[i-1,3]
    p[i] = 0.3 *x[i] + 0.3 *y[i] + 0.3 *z[i] + 0.90 *p[i-1] + Random[i-1,4]
  }
  
  
  Data = as.data.frame(cbind(x,y,z,p))
  Data = Data[-c(1:10),]
  Data_can <- tsdata2canonicalform(Data,nlags)
  
  VARres <- VAR_estim(Data_can, "ols", regstats = FALSE, corank = 0, fmlad = FALSE)
  
  u = VARres$residuals
  reslg <- estimate(t(u))
  
  OM <- cov(u[,reslg$k]) # covariance matrix of reduced form residuals 
  P <- t(chol(OM)) # choleski triangular matrix
  D <- matrix(0,4,4) # scaling matrix to get 1's on diag of W
  diag(D) <- diag(P)
  W <- D %*% solve(P) # rotation matrix
  B0n <- (diag(4) - W)
  v_res=t(solve(W)%*%t(u))
  
  resultchol[[counter]] = B0n
  order[[counter]] = reslg$k
  pinup[counter] = reslg$pinup
  
}

reza1 = rep(0,sim)
reza2 = rep(0,sim)
reza3 = rep(0,sim)
reza4 = rep(0,sim)
reza5 = rep(0,sim)
reza6 = rep(0,sim)

for (counter in 1:sim){
  
  reza1[counter] = resultchol[[counter]][2,1] 
  reza2[counter] = resultchol[[counter]][3,1]
  reza3[counter] = resultchol[[counter]][4,1] 
  reza4[counter] = resultchol[[counter]][3,2] 
  reza5[counter] = resultchol[[counter]][4,2] 
  reza6[counter] = resultchol[[counter]][4,3]
}


orders = matrix(0,sim,5)

for (i in 1:sim){
  orders[i,1] = order[[i]][[1]]
  orders[i,2] = order[[i]][[2]]
  orders[i,3] = order[[i]][[3]]
  orders[i,4] = order[[i]][[4]]
  orders[i,5] = pinup[i]
  }

ss = logical(sim)
for (i in 1:sim){
  if(orders[i,1] == 1 && orders[i,2] == 2 && orders[i,3] == 3 && orders[i,4] == 4){
    ss[i] = TRUE
  }
}

simresult[[name]] = matrix(c(mean(reza1),var(reza1),mean(reza2),var(reza2),mean(reza3),var(reza3)
                                 ,mean(reza4),var(reza4),mean(reza5),var(reza5),mean(reza6),var(reza6),sum(ss==TRUE),sim),7,2,byrow = TRUE)

ordersum[[name]] = matrix(c(sum(orders[,1]==1),sum(orders[,2]==1),sum(orders[,3]==1),sum(orders[,4]==1),
                            sum(orders[,1]==2),sum(orders[,2]==2),sum(orders[,3]==2),sum(orders[,4]==2),
                            sum(orders[,1]==3),sum(orders[,2]==3),sum(orders[,3]==3),sum(orders[,4]==3),
                            sum(orders[,1]==4),sum(orders[,2]==4),sum(orders[,3]==4),sum(orders[,4]==4)),4,4,byrow = TRUE)

print(ordersum)
print(simresult)






