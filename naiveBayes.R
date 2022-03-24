library(rootSolve)
library(fitdistrplus)
# Naive Bayes

naive_thresh = function(sampled,c1,c2,bw=5,betadens=FALSE,makeplot=FALSE,...){
 dat1 = sampled$fibre_het[sampled$clust==1]
 dat2 = sampled$fibre_het[sampled$clust==2]
 datall = sampled$fibre_het

 if(length(dat1)<2) return(100.0)
 if(length(dat2)<2) return(0.0)

 if(betadens){
  b1 = fitdist(dat1/100,"beta")
  b2 = fitdist(dat2/100,"beta")
  ball = fitdist(datall/100,"beta")

  f1 = function(x) dbeta(x/100, b1$estimate[["shape1"]],b1$estimate[["shape2"]])/100 
  f2 = function(x) dbeta(x/100, b2$estimate[["shape1"]],b2$estimate[["shape2"]])/100
  fall = function(x) dbeta(x/100, ball$estimate[["shape1"]],ball$estimate[["shape2"]])/100
 }else{
  d1 = density(dat1,n=1001,from=0,to=100,bw=bw)
  d2 = density(dat2,n=1001,from=0,to=100,bw=bw)
  dall = density(datall,n=1001,from=0,to=100,bw=bw)

  d1$y = d1$y/(sum(d1$y)*(d1$x[2]-d1$x[1]))
  d2$y = d2$y/(sum(d2$y)*(d2$x[2]-d2$x[1]))
  dall$y = dall$y/(sum(dall$y)*(dall$x[2]-dall$x[1]))

  f1 = approxfun(d1$x,d1$y)
  f2 = approxfun(d2$x,d2$y)
  fall = approxfun(dall$x,dall$y)
 }

 f1adj = function(x) f1(x)*c1
 f2adj = function(x) f2(x)*c2

 # Probability cluster 2 given mutation load
 P1 = function(x) c1*f1(x)/(c1*f1(x)+c2*f2(x))
 P2 = function(x) c2*f2(x)/(c1*f1(x)+c2*f2(x))

 diffP = function(x) P2(x)-P1(x)
 candidates = uniroot.all(diffP,c(0,100))
 
 if(length(candidates)==0){ # If no threshold can be found by intersection
  if(diffP(50.0)>0.0){ # Everything is cluster 2 (deficient)
    tres = 0.0
  }else{
    tres = 100.0 # Everything is cluster 1 (normal)
  } 
 }else{
   tres = max(candidates)
 }

 if(makeplot){
  plist = list()
  plist[["kde"]] = function(){
   lwd = 3
   #op = par(mfrow=c(1,2))
   hist(datall,50,freq=FALSE,xlab="Mutation load (%)",cex.axis=1.55,cex.lab=2.25,main="Kernel density estimates of PDF",cex.main=1.55,col="lightgray",border="darkgrey",xlim=c(0,100))
   curve(fall,from=0,to=100,add=TRUE,col="darkgrey",lwd=lwd)
   curve(f1adj,from=0,to=100,add=TRUE,col="blue",lwd=lwd)
   curve(f2adj,from=0,to=100,add=TRUE,col="red",lwd=lwd)
   abline(v=tres,col="grey",lwd=2,lty=2)
   legend("topleft",lwd=lwd,col=c("blue","red","darkgrey"),legend=c("Normal fibres","Deficient fibres","All fibres"),cex=1.55,bty = "n")
  }
  
  plist[["naive"]] = function(){
   lwd = 3
   curve(P1,from=0,to=100,ylim=c(0,1),ylab="Probability of belonging to cluster",xlab="Mutation load (%)",lwd=lwd,cex.axis=1.55,cex.lab=2.25,main="Naive Bayes classifier\nthreshold estimate",cex.main=1.55,col="blue")
   curve(P2,from=0,to=100,add=TRUE,col="red",lwd=lwd)
   abline(h=0.5,col="grey",lwd=2,lty=2)
   abline(v=tres,col="grey",lwd=2,lty=2)
   #legend("topleft",lwd=lwd,col=c("blue","red"),legend=c("Normal fibres","Deficient fibres"),bg="white")
   #par(op)
  }
  return(plist)
 }else{
  return(tres)
 }
}






