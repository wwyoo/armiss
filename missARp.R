#function to compute AR(p) covariance matrix
covmat<-function(phi,sigma2){
 if(length(phi)==0){  #to avoid error if no AR order selected
  Sigma<-sigma2*diag(rep(1,T))
 }
 else{
  Corr<-toeplitz(ARMAacf(ar=phi,lag.max=T-1))
  gamma0<-sigma2/(1-crossprod(phi,Corr[1,2:(1+length(phi))]))
  Sigma<-as.numeric(gamma0)*Corr
 }
return(Sigma)
}

#compute AR coefficients from pacf using Durbin-Levinson 
pacf2acf<-function(x){
p<-length(x)
B<-matrix(0,p,p)
diag(B)<-x
if(p==1){
  phi<-B[, p]
  }
else{
  for(i in 2:p){
    B[1:(i-1),i]<-B[1:(i-1),(i-1)]-B[i,i]*rev(B[1:(i-1),(i-1)])
    }
  phi<-B[, p]
  }
return(phi)
}

#negative log-likelihood of the observed data
loglike<-function(theta,PO,XO){
muO<-theta[1]*rep(1,length(XO))
tpacf<-theta[3:length(theta)]
acf<-pacf2acf(tanh(tpacf))
SigmaOO<-PO%*%covmat(phi=acf,sigma2=exp(theta[2]))%*%t(PO)
browser()
f<--dmvnorm(XO,mean=muO,sigma=SigmaOO,log=TRUE)
return(f)
}

aic<-function(k,PO,XO){
#initial estimates from Yule-Walker
init<-ar.yw(XO,aic=FALSE,order.max=k,demean=TRUE)
mu0<-init$x.mean
tpacf0<-atanh(runif(k,-1,1))
tsigma20<-log(init$var.pred)
theta0<-c(mu0,tsigma20,tpacf0)

#maximize likelihood of only the observed data (on transformed space)
result<-optim(theta0, loglike, PO=PO, XO=XO)
aick<-2*result$value+2*(k+2)
estimatek<-result$par
codek<-result$convergence
theta<-c(aic=aick,code=codek,estimate=estimatek)
return(theta)
}

#----------------------------------------------------------------------------
#function to do imputation 
ar.miss2<-function(data, order=NA, max.k=12){
  T <- length(data)
  oidx <- which(is.na(data) == FALSE)  #observed data index
  permu <- c(oidx, (1:T)[-oidx])
  tobs <- length(oidx)

  P <- elem(permu)  #permutation matrix
  PO <- P[1:tobs, ]
  XO <- data[oidx]  #observed data

if(is.na(order)==TRUE){
  for(i in 1:max.k){
    res[i,1:(4+i)]<-aic(k=i, PO=PO, XO=XO) 
    }
  ind<-which(res[,1]==min(res[,1]))
  theta<-res[ind, 3:(4+ind)]
}
else{
   res<-aic(k=order, PO=PO, XO=XO)
   theta<-res[3:(4+order)]
}

thetaest<-c(mu=theta[1],phi=pacf2acf(tanh(theta[3:length(theta)])),sigma2=exp(theta[2]))
return(thetaest)
}

