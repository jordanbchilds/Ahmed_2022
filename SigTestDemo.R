# Make some fake data
P01 = data.frame(val = c(rnorm(50,50,1),rnorm(50,70,5)), lab = rep(c("A","B"),each=50))
P02 = data.frame(val = c(rnorm(40,41,8),rnorm(50,70,7.5)), lab = c(rep("A",40),rep("B",50)))

# Fit kernel density estimates to data with each label and find intersection
estimateSwitch = function(subj){
     dA = density(subj$val[subj$lab=="A"],n=1001,from=0,to=100,bw=5)
     dB = density(subj$val[subj$lab=="B"],n=1001,from=0,to=100,bw=5)

     crosses = abs(diff(sign(dB$y-dA$y)))>0
     threshind = which(crosses)[which.max(dA$y[which(crosses)+1])]
     thresh = dA$x[threshind]
     threshdens = dA$y[threshind]
     return(c(thresh,threshdens))
}

# Bootrap above to get uncertainty about intersection
bootstrapSwitch = function(subj){
   inds = sample(1:dim(subj)[1],replace=TRUE)
   return(estimateSwitch(subj[inds,])[1])
}

# Point estimates for two fake subjects
sw01 = estimateSwitch(P01)
sw02 = estimateSwitch(P02)

plot(density(P01$val),xlim=c(0,100))
abline(v=sw01[1],lty=2,lwd=3)
points(density(P02$val),type="l",col="red")
abline(v=sw02[1],lty=2,lwd=3,col="red")

# Bootstrap measure of uncertainty about location of intersections
bs01 = replicate(1000,bootstrapSwitch(P01))
bs02 = replicate(1000,bootstrapSwitch(P02))

plot(density(bs01),xlim=range(c(bs01,bs02)),main="",xlab="switch")
points(density(bs02),type="l",col="red")

# Want to test whether observed differences are significant
# Construct test statistic based on boostrap estimates above
teststat = abs((bs01 - sw01[1]) - (bs02 - sw02[1]))
swdiff = abs(sw01[1]-sw02[1])
pval = sum(teststat > swdiff)/length(teststat)
plot(density(teststat),main=paste("p-value:",signif(pval,3)))
abline(v=swdiff,lwd=2,lty=2)
