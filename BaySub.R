#####################################################################
##### BaySub algorithm 
##### Author Yetian Fan
##### Data 2021-07-20
##### The pairs of M-value in Data.Rdata
#####################################################################
BaySub <- function(nstep,k,initialtimes){
library(gtools)

data = load("Data.Rdata")

Variance = 1
seednum = ceiling(runif(1,0,1)*10^8)
set.seed(seednum)
meantruncate = -50
sdtruncate = 0
lengthtruncate =100

n = dim(X)[1]
m = dim(X)[2]

##### initial the parameters
### Z[i,j], W[i,j], alpha[s], ita[s], mu[s], sigmaval[s], gammaval, delta, Prolambda[k], lambda[n], R[j]
Z = round(matrix(runif(n*m,0,1),nrow=n))
W = round(matrix(runif(k*m,0,1),nrow=k))
alpha = matrix(NA,nrow=2,ncol=nstep)
alpha[,1] = c(-1,1)*initialtimes
ita = matrix(NA,nrow=2,ncol=nstep)
ita[,1] = c(1,1)*initialtimes
mu = matrix(NA,nrow=2,ncol=nstep)
mu[,1] = c(-1,1)*initialtimes
sigmaval = matrix(NA,nrow=2,ncol=nstep)
sigmaval[,1] = c(1,1)*initialtimes
gammaval = matrix(NA,nrow=1,ncol=nstep)
gammaval[1] = 0.5
delta = matrix(NA,nrow=1,ncol=nstep)
delta[1] = 1-rbeta(1,1,100)
betaval = matrix(NA,nrow=1,ncol=nstep)
betaval[1] = 0.5
Prolambda = matrix(NA,nrow=k,ncol=nstep)
Prolambda[,1] = rep(1/k,1,k)
lambda = matrix(NA,nrow=n,ncol=nstep)
lambda[,1] = ceiling(matrix(runif(n,0,1)*k,nrow=n,ncol=1))
R = matrix(NA,nrow=m,ncol=nstep)
R[,1] = round(runif(m,0,1))
P = matrix(NA,nrow=7, ncol=nstep-1)

for (s in 2:nstep){
  ##### update gammaval #####
  gammaval[s] = rbeta(1,sum(R[,s-1]==1)+1,sum(R[,s-1]==0)+1)
  
  ##### update betaval #####
  betaval[s] = rbeta(1,sum(W==1)+1,sum(W==0)+1)
  
  ##### update delta #####
  value1 = 0
  for (j in 1:m){
    value1 = value1 + sum(Z[,j]==R[j,s-1])
  }
  delta[s] = rbeta(1,value1+1,n*m-value1+1)
  
  ##### update R #####
  value1 = 0   
  for (j in 1:m){
    Prob = c(0,0)
    value1 = sum(Z[,j]==0)
    Prob[1] = log(delta[s])*value1+log(1-delta[s])*(n-value1)+log(1-gammaval[s])
    Prob[2] = log(delta[s])*(n-value1)+log(1-delta[s])*value1+log(gammaval[s])
    Pro = 1/(1+exp(Prob[2]-Prob[1]))
    u=runif(1,0,1)
    R[j,s] = (u>Pro)+0
  }
  
  ##### update Prolambda #####
  value1 = matrix(NA,nrow=k,ncol=1)   
  for(l in 1:k){
    value1[l] = sum(lambda[,s-1]==l)+1
  }
  Prolambda[,s] = rdirichlet(1,value1)
  
  ##### calculate the Norm matrix of  X and Y
  logNX = array(NA,dim=c(n,m,2))
  logNY = array(NA,dim=c(n,m,2))
  for(i in 1:n){
    for(j in 1:m){
      for(l in 1:2){
        logNY[i,j,l] = -(Y[i,j]-alpha[l,s-1])^2/(2*ita[l,s-1]) - 1/2*log(2*pi*ita[l,s-1])
        logNX[i,j,l] = -(X[i,j]-mu[l,s-1])^2/(2*sigmaval[l,s-1]) - 1/2*log(2*pi*sigmaval[l,s-1])
      }
    }
  }
  
  ##### update lambda ##### 
  for(i in 1:n){
    Prob = matrix(0,nrow=1,ncol=k)
    for(l in 1:k){
      for(j in 1:m){
        Boolean = ((Z[i,j]!=W[l,j])+1)
        Prob[l] = Prob[l]+logNY[i,j,Boolean]
      }
      Prob[l] = Prob[l]+ log(Prolambda[l,s])
    }
    lambda.try=ceiling(runif(1,0,1)*k)
    r=min(Prob[lambda.try]-Prob[lambda[i,s-1]],0)
    u=log(runif(1,0,1))
    if(u<r)
      lambda[i,s]=lambda.try
    else
      lambda[i,s]=lambda[i,s-1]
  }
  
  
  ##### update W[k,j] ####
  PW = matrix(0,nrow=k,ncol=m)
  PW.try = matrix(0,nrow=k,ncol=m)
  W.try = matrix(NA,nrow=k,ncol=m)
  for(i in 1:k){
    for(j in 1:m){
      PW[i,j] = PW[i,j] + W[i,j]*log(betaval[s])+(1-W[i,j])*log(1-betaval[s])
      
      W.try[i,j] = floor(runif(1,0,1)*2)
      PW.try[i,j] = PW.try[i,j] + W.try[i,j]*log(betaval[s])+(1-W.try[i,j])*log(1-betaval[s])
    }
  }
  for(i in 1:n){
    id = lambda[i,s]
    for(j in 1:m){
      Boolean = ((Z[i,j]!=W[id,j])+1)
      PW[id,j] = PW[id,j] + logNY[i,j,Boolean]
      
      Boolean = ((Z[i,j]!=W.try[id,j])+1)
      PW.try[id,j] = PW.try[id,j] + logNY[i,j,Boolean]
    }
  }
  for(i in 1:k){
    for(j in 1:m){
      r=min(PW.try[i,j]-PW[i,j],0)
      u=log(runif(1,0,1))
      if(u<r)
        W[i,j] = W.try[i,j]
    }
  }
  
  ##### update Z[i,j] #### 
  PZ = array(0,dim=c(n,m,2))
  for(i in 1:n){
    for(j in 1:m){
      PZ[i,j,1] = logNY[i,j,W[lambda[i,s],j]+1] +  logNX[i,j,1] 
      PZ[i,j,2] = logNY[i,j,2-W[lambda[i,s],j]] +  logNX[i,j,2]
      if(R[j,s]==0){
        PZ[i,j,1] = PZ[i,j,1] + log(delta[s]) 
        PZ[i,j,2] = PZ[i,j,2] + log(1-delta[s])
      }else{
        PZ[i,j,1] = PZ[i,j,1] + log(1-delta[s])
        PZ[i,j,2] = PZ[i,j,2] + log(delta[s])
      }
    }
  }
  for(i in 1:n){
    for(j in 1:m){
      Prob = 1/(1+exp(PZ[i,j,2]-PZ[i,j,1]))
      u=runif(1,0,1)
      Z[i,j] = (u>Prob)+0
    }
  }
  
  ##### update mu and sigmaval #####
  sumX = c(0,0) 
  sumX2 = c(0,0) 
  countNum = c(0,0) 
  for(i in 1:n){
    for(j in 1:m){
      sumX[Z[i,j]+1] = sumX[Z[i,j]+1] + X[i,j]
      sumX2[Z[i,j]+1] = sumX2[Z[i,j]+1] + X[i,j]^2
      countNum[Z[i,j]+1] = countNum[Z[i,j]+1] +1
    }
  }
  for(i in 1:2){
    mu[i,s] = rnorm(1,mean=(sumX[i]/countNum[i]),sd=(sqrt(sigmaval[i,s-1]/countNum[i])))
    while((mu[i,s] < meantruncate)|(mu[i,s] > meantruncate+lengthtruncate)){
      mu[i,s] = rnorm(1,mean=(sumX[i]/countNum[i]),sd=(sqrt(sigmaval[i,s-1]/countNum[i])))
    }
    Gamalpha = (countNum[i]/2)+1
    Gambeta = 0.5*sumX2[i]-mu[i,s]*sumX[i]+0.5*countNum[i]*mu[i,s]^2
    sigmaval[i,s] = (1/rgamma(1,Gamalpha,Gambeta))
    while((sigmaval[i,s] < sdtruncate)|(sigmaval[i,s] > sdtruncate+lengthtruncate)){
      sigmaval[i,s] = (1/rgamma(1,Gamalpha,Gambeta))
    }
  }
  
  ##### calculate phi[i,j]  #####
  phi = matrix(NA,nrow=n,ncol=m)
  for(i in 1:n){
    for(j in 1:m){
      phi[i,j] = (Z[i,j]!=W[lambda[i,s],j])
    }
  } 
  
  ##### update alpha and ita #####
  alpha.try = c(0,0)
  alpha.try[1] = rnorm(1,mean=alpha[1,s-1],sd=Variance)
  alpha.try[2] = rnorm(1,mean=alpha[2,s-1],sd=Variance) 
  while((alpha.try[1] < meantruncate)|(alpha.try[1] > meantruncate+lengthtruncate)){
    alpha.try[1] = rnorm(1,mean=alpha[1,s-1],sd=Variance)
  }
  while((alpha.try[2] < meantruncate)|(alpha.try[2] > meantruncate+lengthtruncate)){
    alpha.try[2] = rnorm(1,mean=alpha[2,s-1],sd=Variance)
  }
  Palpha = c(0,0)
  Palpha.try = c(0,0)
  for (i in 1:n){
    for (j in 1:m){
      Palpha[phi[i,j]+1] = Palpha[phi[i,j]+1] - 0.5*(Y[i,j]-alpha[phi[i,j]+1,s-1])^2/ita[phi[i,j]+1,s-1]
      Palpha.try[phi[i,j]+1] = Palpha.try[phi[i,j]+1] - 0.5*(Y[i,j]-alpha.try[phi[i,j]+1])^2/ita[phi[i,j]+1,s-1]
    }
  }
  for(i in 1:2){
    r = min(Palpha.try[i]-Palpha[i],0)
    u=log(runif(1,0,1))
    if(u<r){
      alpha[i,s] = alpha.try[i]
    }else{
      alpha[i,s] = alpha[i,s-1]
    }
  }
  
  sumY = c(0,0)
  sumY2 = c(0,0)
  countNum = c(0,0)
  for(i in 1:n){
    for(j in 1:m){
      sumY[phi[i,j]+1] = sumY[phi[i,j]+1] + Y[i,j]
      sumY2[phi[i,j]+1] = sumY2[phi[i,j]+1] + Y[i,j]^2
      countNum[phi[i,j]+1] = countNum[phi[i,j]+1] +1
    }
  }
  
  for(i in 1:2){
    Gamalpha = (countNum[i]/2)+1
    Gambeta = 0.5*sumY2[i]-alpha[i,s]*sumY[i]+0.5*countNum[i]*alpha[i,s]^2
    ita[i,s] = (1/rgamma(1,Gamalpha,Gambeta))
    while((ita[i,s] < sdtruncate)|(ita[i,s] > sdtruncate+lengthtruncate)){
      ita[i,s] = (1/rgamma(1,Gamalpha,Gambeta))
    }
  }
  
  ###### calculate the likelihood #######
  count1 = sum(phi==0)
  count2 = n*m - count1
  P1 = -0.5*(m*n)*log(2*pi) - 0.5*count1*log(ita[1,s]) - 0.5*count2*log(ita[2,s])
  for (i in 1:n){
    for (j in 1:m){
      P1 = P1 - 0.5*(Y[i,j]-alpha[phi[i,j]+1,s])^2/ita[phi[i,j]+1,s]
    }
  }
  
  count1 = sum(Z==0)
  count2 = n*m - count1
  P2 = -0.5*(m*n)*log(2*pi) - 0.5*count1*log(sigmaval[1,s]) - 0.5*count2*log(sigmaval[2,s])
  for (i in 1:n){
    for (j in 1:m){
      P2 = P2 - 0.5*(X[i,j]-mu[Z[i,j]+1,s])^2/sigmaval[Z[i,j]+1,s]
    }
  }
  
  count1 = matrix(NA,nrow=1,ncol=k)
  for(i in 1:k){
    count1[i]=sum(lambda[,s]==i)
  }
  P3 = sum(log(Prolambda[,s])*count1)
  
  count1 = 0
  for (j in 1:m){
    count1 = count1 + sum(Z[,j]==R[j,s])
  }
  P4 = log(delta[s])*count1+(m*n-count1)*log(1-delta[s])
  
  count1 = sum(R[,s]==1)
  P5 = log(gammaval[s])*count1+(m-count1)*log(1-gammaval[s])
  
  count1 = sum(W==1)
  P6 = log(betaval[s])*count1 + (k*m - count1)*log(1-betaval[s])
  
  P[1,s-1] = P1
  P[2,s-1] = P2 
  P[3,s-1] = P3
  P[4,s-1] = P4
  P[5,s-1] = P5 
  P[6,s-1] = P6
  P[7,s-1] = P1 + P2 + P3 + P4 + P5 + P6
}

####### calculate the accuracy  #######
listNum = permutations(k,k)
NoCompare = dim(listNum)[1]
TestAccuracy = matrix(NA,nrow=1,ncol=nstep)
IDMax = matrix(NA,nrow=1,ncol=nstep)
for(s in 1:nstep){
  Accuracy = matrix(NA,nrow=1,ncol=NoCompare)
  for (i in 1:NoCompare){
    CompareNum = lambda[,s]
    for (j in 1:k){
      id = which(lambda[,s]==j)
      CompareNum[id] = listNum[i,j]
    }
    Accuracy[i] = sum(CompareNum==Target)/n
  }
  TestAccuracy[s] = max(Accuracy)
  IDMax[s] = which.max(Accuracy)
}

return(list(TestAccuracy=TestAccuracy,lambda=lambda[,nstep],W=W))
}

################################################################################
################################################################################

ptm <- proc.time()

k = 3                     # different paths from the corresponding normal tissue
nstep = 200               # number of iterations
initialtimes = 1          # initial value for variables, which can be selected from the set {1,1.2,1.4,1.6,1.8,2,...}.

Result = BaySub(nstep,k,initialtimes)
Accuracy = Result$TestAccuracy
lambda = Result$lambda
W = Result$W

print(proc.time() - ptm)