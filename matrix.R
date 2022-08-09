# Setup---------
# library(coda)
set.seed(407)
iterations = 500000
m = 20
n = 20
p = 0.25
data = matrix(rbinom(n*m,1,p), nrow=m, ncol=n)
r = rowSums(data)
c = colSums(data)
w = matrix(runif(n*m), nrow=m, ncol=n)
w2 = matrix(rbeta(n*m,1,3), nrow=m, ncol=n)
# Functions-------
swapprob = function(i1,i2,j1,j2,A,w) {
  pswap = 0
  if (A[i1,j1]==1 & A[i2,j2]==1) {
    pswap = w[i1,j2]*w[i2,j1]/(w[i1,j1]*w[i2,j2]+w[i1,j2]*w[i2,j1])
  }
  else if (A[i1,j1]==0 & A[i2,j2]==0) {
    pswap = w[i1,j1]*w[i2,j2]/(w[i1,j1]*w[i2,j2]+w[i1,j2]*w[i2,j1])
  }
  return(pswap)
}
swap = function(i1,i2,j1,j2,A) {
  temp = A[c(i1,i2),j2]
  A[c(i1,i2),j2] = A[c(i1,i2),j1]
  A[c(i1,i2),j1] = temp
  return(A)
}
diagdiv = function(A){
  ones = which(A==1, arr.ind=TRUE)
  div = sum(abs(ones[,1]-ones[,2]))
  div = div/(n*sum(A))
  return(div)
}
cscore = function(A){
  c = 0
  m = nrow(A)
  r = rowSums(A)
  for (i in 2:m){
    for (j in 1:i){
      s = sum(A[i,]==1 & A[j,]==1)
      c = c + (r[i]-s)*(r[j]-s)
    }
  }
  return(2*c/(n*(n-1)))
}
s2bar = function(A){
  s = 0
  m = nrow(A)
  S = A %*% t(A)
  for (i in 2:m){
    for (j in (1:(i-1))){
      s = s + S[i,j]^2
    }
  }
  return(2*s/(m*(m-1)))
}

# Curveball (function taken directly from Fout[2020])
curveball_trade <- function(A, w){
  # Weighted curveball swapping
  # A: binary matrix
  # W: weight matrix (positive elements, same size as A)
  
  # sample two rows at random
  rows <- sample(1:nrow(A))
  r1 <- A[rows[1],]
  r2 <- A[rows[2],]
  w1 <- w[rows[1],]
  w2 <- w[rows[2],]
  # get columns shared and not shared
  A12 <- which(r1+r2==2)  # Shared between r1 and r2
  A1m2 <- which(r1-r2>0 & w2>0)  # in r1 but not r2 (and not a struct zero in 2)
  A2m1 <- which(r2-r1>0 & w1>0)  # in r2 but not r1 (and not a struct zero in 1)
  
  A1m2sz <- which(r1-r2>0 & w2==0)  # in r1 but not r2 (and a struct zero in 2)
  A2m1sz <- which(r2-r1>0 & w1==0)  # in r2 but not r1 (and a struct zero in 1)
  
  # Check if can trade
  can_trade <- !(length(A1m2)==0 | length(A2m1)==0)
  if(!can_trade){
    return(A)
  }
  # shuffle the tradable indices
  B <- sample(c(A1m2, A2m1))
  B1m2 <- B[1:length(A1m2)]
  B2m1 <- B[(length(A1m2)+1):length(B)]
  # compute trade probability
  r1_01 <- B1m2[is.na(match(B1m2, A1m2))]  # cols that went 0 -> 1 in row 1
  r2_01 <- B2m1[is.na(match(B2m1, A2m1))]  # cols that went 0 -> 1 in row 2
  p_trade <- prod(w[rows[1],r1_01])*prod(w[rows[2],r2_01])/(prod(w[rows[1],r1_01])*prod(w[rows[2],r2_01]) + prod(w[rows[2],r1_01])*prod(w[rows[1],r2_01]))
  if(runif(1) < p_trade){
    # construct traded rows
    r1B <- rep(0, length(r1))
    r1B[c(A12, A1m2sz, B1m2)] <- 1  # shared 1's, sz 1's, and traded 1's
    r2B <- rep(0, length(r2))
    r2B[c(A12, A2m1sz, B2m1)] <- 1  # shared 1's, sz 1's, and traded 1's
    
    A[rows[1],] <- r1B
    A[rows[2],] <- r2B
  }
  return(A)
}
# MCMC-------
stat = vector(length=iterations)
stat2 = vector(length=iterations)
swapcount=0
A = data
B = data
for (z in 1:iterations) {
# Rectangle Loop
i1 = sample.int(m, size=1)
j1 = sample.int(n, size=1)
if (A[i1,j1] == 1) {
  j2 = sample(which(A[i1,]==0),1)
  i2 = sample(which(A[,j2]==1),1)
}
else if (A[i1,j1] == 0) {
  i2 = sample(which(A[,j1]==1),1)
  j2 = sample(which(A[i2,]==0),1)
}
# RL Swap
if (A[i2,j2]==A[i1,j1] & A[i1,j2]==A[i2,j1]) {
  u = runif(1)
  if (u < swapprob(i1,i2,j1,j2,A,w)) { 
    A=swap(i1,i2,j1,j2,A)
    swapcount = swapcount+1
  }
}
stat[z] = cscore(A)
B = curveball_trade(B,w)
stat2[z] = cscore(B)
}
# Analysis-----
mcmcstat1 = mcmc(stat)
effectiveSize(mcmcstat1)
mcmcstat2 = mcmc(stat2)
effectiveSize(mcmcstat2)