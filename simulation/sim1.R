require(fastICA)

result = list()
resultchol = list()
order = list()
sim = 100
for (counter in 1:sim){

n = 300
x = rep(0,n)
y = rep(0,n) 
z = rep(0,n)
p = rep(0,n)

Random = cbind(runif(n-1,-0.5,0.5),runif(n-1,-0.5,0.5),runif(n-1,-0.5,0.5),
               runif(n-1,-0.5,0.5))
#Random = cbind(rlnorm(n-1,0,0.5),rlnorm(n-1,0,0.8),rlnorm(n-1,0,0.5),
               #rlnorm(n-1,0,0.5))

for (i in 2:n){
  x[i] = 0.90 *x[i-1] + Random[i-1,1]
  y[i] = 0.3 *x[i] + 0.90 *y[i-1] + Random[i-1,2]
  z[i] = 0.3 *x[i] + 0.3 *y[i] + 0.90 *z[i-1] + Random[i-1,3]
  p[i] = 0.3 *x[i] + 0.3 *y[i] + 0.3 *z[i] + 0.90 *p[i-1] + Random[i-1,4]
}


Data = as.data.frame(cbind(x,y,z,p))
Data = Data[-c(1:10),]
#Data_can <- tsdata2canonicalform(Data,nlags)


# ------ Estimate VAR-LiNGAM ------
#result <- VARLiNGAM(Data_can,"ols", pruning=FALSE, corank=0,ntests = FALSE)
u = resid(VAR(Data,10))
est = estimate(t(u))

OM <- cov(u) # covariance matrix of reduced form residuals 
P <- t(chol(OM)) # choleski triangular matrix
D <- matrix(0,4,4) # scaling matrix to get 1's on diag of W
diag(D) <- diag(P)
W <- D %*% solve(P) # rotation matrix
B0n <- (diag(4) - W)

resultchol[[counter]] = B0n
result[[counter]] = est$B
order[[counter]] = est$k

}

reza1 = rep(0,sim)
reza11 = rep(0,sim)
reza2 = rep(0,sim)
reza22 = rep(0,sim)
reza3 = rep(0,sim)
reza33 = rep(0,sim)
reza4 = rep(0,sim)
reza44 = rep(0,sim)

for (counter in 1:sim){
  
  reza1[counter] = result[[counter]][2,1] 
  reza11[counter] = resultchol[[counter]][2,1]
  reza2[counter] = result[[counter]][3,1]
  reza22[counter] = resultchol[[counter]][3,1]
  reza3[counter] = result[[counter]][4,1] 
  reza33[counter] = resultchol[[counter]][4,1]
  reza4[counter] = result[[counter]][3,2] 
  reza44[counter] = resultchol[[counter]][3,2]
  reza5[counter] = result[[counter]][4,2] 
  reza55[counter] = resultchol[[counter]][4,2]
  reza6[counter] = result[[counter]][4,3]
  reza66[counter] = resultchol[[counter]][4,3]
  
  rezaord[counter] = order[[counter]][[1]]
}


var(reza1)-var(reza11)
var(reza2)-var(reza22)
var(reza3)-var(reza33)
var(reza4)-var(reza44)
var(reza5)-var(reza55)
var(reza6)-var(reza66)
sum(rezaord == 1)
sum(rezaord == 2)
sum(rezaord == 3)
sum(rezaord == 4)

orders = matrix(0,sim,4)

for (i in 1:sim){
orders[i,1] = order[[i]][[1]]
orders[i,2] = order[[i]][[2]]
orders[i,3] = order[[i]][[3]]
orders[i,4] = order[[i]][[4]]}

ss = logical(sim)
for (i in 1:sim){
  if(orders[i,1] == 1 && orders[i,2] == 2 && orders[i,3] == 3 && orders[i,4] == 4){
    ss[i] = TRUE
  }
}

#result[["Bhat"]][[1]]

