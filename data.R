T<-365
n<-200
ar1<-arima.sim(n=T,list(ar=c(0.5,0.2), sd=1)
index<-sample(1:T,(T-n),replace=FALSE) 
ar1[index]<-NA
ar1sim<-as.ts(ar1)
save(ar1sim, file="ar1sim.RData", compress=TRUE)










data(rec)
year<-rep(1:38, each=12)[1:length(rec)]
newfish <- cbind(as.vector(rec), year)
index<-strata(newfish, "year", rep(3, 38), "srswor")$ID_unit
newfish[index]<-NA
newfish<-as.ts(newfish[,1])


trec<-diff(log(rec))
f1<-1
W<-cbind(cos(2*pi*f1*time(trec)), sin(2*pi*f1*time(trec)) )
rec.p<-lm(trec~W)$resid

year<-rep(1:38, each=12)[1:length(rec.p)]
newfish <- cbind(as.vector(rec.p), year)
index<-strata(newfish, "year", rep(3, 38), "srswor")$ID_unit
newfish[index]<-NA
newfish<-as.ts(newfish[,1])

data<-lsfit$resid
ar.mle(data)

source('missARp.R')

ar.missf(data)
ar.miss2(data, order=1)
Rprof()
system.time( car.miss2(data, order=1) )
Rprof(NULL)

summaryRprof("Rprof.out")

ar.missf<-cmpfun(ar.miss)
set.seed(1234)
ar.miss(ar1.sim)
