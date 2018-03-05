#file: missAR1.R
#brief: Estimate AR(1) process with missing data
#author: William Weimin Yoo

#fn 

library(mvtnorm)  #package to sample multivariate normal

#Below are functions written to calculate different statistics
#-----------------------------------------------------------------------------
#compute permutation matrix
elem<-function(permu){
P<-matrix(0,T,T)
permu<-cbind(1:T,permu)
P[permu]<-1
return(P)
}

#-------------------------------------------------------------------------------
#calculate initial estimate of rho
rho<-function(X){
Xbarnolast<-mean(X[-length(X)])
Xbarnofirst<-mean(X[-1])
cross<-sum((X[-1]-Xbarnofirst)*(X[-length(X)]-Xbarnolast))  #Kronecker product
var<-sum((X[-length(X)]-Xbarnolast)^2)
rho<-cross/var
return(rho)
}

#-------------------------------------------------------------------------------
#compute initial estimate of sigma2
sigma2<-function(X,mu,rhohat){
sum((X[-1]-mu-rhohat*(X[-length(X)]-mu))^2)/(length(X)-3)
}

#------------------------------------------------------------------------------
#compute mu, rho and sigma2 estimates
obj<-function(theta,PO,XO){
muO<-theta[1]*PO%*%rep(1,T)
rh<-tanh(theta[2])
Sigma<-(exp(theta[3])/(1-rh^2))*rh^abs(row(Sigma0)-col(Sigma0))
SigmaOO<-PO%*%Sigma%*%t(PO)
f<--dmvnorm(XO,mean=muO,sigma=SigmaOO,log=TRUE)
return(f)
}

#----------------------------------------------------------------------------
#function to do imputation 
ar.miss<-function(data){
T<-length(data)
Sigma0<-matrix(0,T,T)
permu<-c(which(data!=88888),which(data==88888))
missing<-length(which(data==88888))

P<-elem(permu)  #permutation matrix
PO<-P[1:(T-missing),]
XO<-data[data!=88888]  #observed data

#initial estimates
mu0<-mean(XO)
trho0<-atanh(rho(XO))
tsigma20<-log(sigma2(XO,mu=mu0,rho=rho(XO)))

#maximize likelihood of only the observed data (on transformed space)
theta<-nlm(obj, c(mu0,trho0,tsigma20),PO=PO,XO=XO)$estimate

thetaest<-c(mu=theta[1],rho=tanh(theta[2]),sigma2=exp(theta[3]))
return(thetaest)
}














