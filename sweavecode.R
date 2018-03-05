###########################################################################
T=365  #number of data points
m<-7   #moving average window length
B<-1000  #number of parametric bootstrap iteration
n<-c(50,100,200,250)  #number of observed data
rho.true<-c(0,0.5,0.75,0.95)  #true AR(1) coefficient
com<-length(n)*length(rho.true)  #total combination
mu.true<-0  #true mean
sigma2.true<-1  #true innovation variance
Sigma0<-matrix(0,T,T)

#store intermediate results
resulttrue<-matrix(0,length(rho.true),2)
result<-matrix(0,com,2)

set.seed(1234)  #set seed
 for(i in 1:4){
  for(j in 1:4){
   index<-sample(1:T,(T-n[i]),replace=FALSE)  #assume missing at random      
   ar1<-suppressWarnings(arima.sim(n=T,list(ar=c(rho.true[j])),sd=sqrt(sigma2.true))) #sample from AR(1)
   ar1.m<-ar1
   ar1.m[index]<-88888   #represents missing values
   resulttrue[j,]<-boot.est(mu=mu.true,rho=rho.true[j],sigma2=sigma2.true) #alpha0 and beta0 in article
   u<-j+length(rho.true)*(i-1)
   result[u,]<-method(ar1.m)
  }
}












