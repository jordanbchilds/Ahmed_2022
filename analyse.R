#library(devtools)
#install_github('gastonstat/colortools')

library(colortools)
library(data.table)
library(fitdistrplus)
library(gmp)
library(mclust)
library(Hotelling)
library(rootSolve)
library(beanplot)
source("naiveBayes.R")

threshFromScratch = function(dtf,predgmm,c1=0.5,c2=0.5,makeplot=FALSE){
  miniA = makemini(dtf)

  probclass = predgmm(miniA)$z
  probs = probclass[,1]

  dtf$clust = 0
  pcutoff = 0.6
  dtf$clust[probs>pcutoff] = 1
  dtf$clust[probs<(1-pcutoff)] = 2
  dtf$probs = probs

  dtf$biochem[dtf$clust==0] = "Unknown"
  dtf$biochem[dtf$clust==1] = "Normal"
  dtf$biochem[dtf$clust==2] = "Deficient"

  res = naive_thresh(dtf,c1,c2,makeplot=makeplot)
  return(res[1])
}

set.seed(42)

usecols = c("adj_CI","adj_CIV","adj_porin")

# https://towardsdatascience.com/unsupervised-machine-learning-clustering-analysis-d40f2b34ae7e

dat = fread("rawdat_a.csv",sep=",",stringsAsFactors=FALSE)
dat = dat[nchar(dat$caseno)==3,]

#dat$caseno[dat$caseno=="P08"]=paste(dat$caseno[dat$caseno=="P08"],dat$Batch[dat$caseno=="P08"],sep="_")
dat = dat[!((dat$caseno=="P08")&(dat$Batch==3)),]

probexact = function(ndef,nrand,allfibs,deficient) chooseZ(deficient,ndef)*chooseZ(allfibs-deficient,nrand-ndef)/chooseZ(allfibs,nrand)
probatleast = function(ndef,nrand,allfibs,deficient) sum(probexact(ndef:deficient,nrand,allfibs,deficient))
samplesneeded = function(ndef,allfibs,deficient){
  nrand = ndef-1
  patleast = 0
  while((patleast<0.5)&(nrand<allfibs)){
    nrand = nrand + 1
    patleast = probatleast(sampdeficient,nrand,allfibs,deficient)
  }
  return(nrand)
}

dat$pc[dat$controls=="control"]="c"
dat=dat[!is.na(dat$raw_porin),]

quant = 0.025
deltas = c(quantile(dat$raw_porin,quant,na.rm=TRUE),quantile(dat$raw_CI,quant,na.rm=TRUE),quantile(dat$raw_CIV,quant,na.rm=TRUE))
names(deltas) = c("raw_porin","raw_CI","raw_CIV")

lower = -9999999
dat$adj_porin = pmax(lower,dat$raw_porin-deltas["raw_porin"])
dat$adj_CI = pmax(lower,dat$raw_CI-deltas["raw_CI"])
dat$adj_CIV = pmax(lower,dat$raw_CIV-deltas["raw_CIV"])
dat$adj_CI = dat$raw_CI-deltas["raw_CI"]
dat$adj_CIV = dat$raw_CIV-deltas["raw_CIV"]
dat$theta_CI = 360*atan(dat$adj_CI/dat$adj_porin)/(2*pi)
dat$theta_CIV = 360*atan(dat$adj_CIV/dat$adj_porin)/(2*pi)
dat$rat_CI = dat$raw_CI/dat$raw_porin
dat$rat_CIV = dat$raw_CIV/dat$raw_porin

dat$pd_CI = NA
dat$pd_CIV = NA
# Calculate control regression coefficients
rcoef = list()
dt = dat[dat$controls=="control",]
for(comp in c("CI","CIV")){
   x = dt$raw_porin
   y = dt[[paste("raw",comp,sep="_")]]

    pc = prcomp( ~ x + y)
    slope = pc$rotation[2,1] / pc$rotation[1,1]
    intercept = pc$center[2] - slope*pc$center[1]
    #dat[[paste("pd",comp,sep="_")]][dat$caseno==case] = -(intercept+slope*x-y)/sqrt(1+slope^2)
    rcoef[[comp]] = list(slope=slope,intercept=intercept)
}
for(case in sort(unique(dat$caseno))){
  x = dat$raw_porin[dat$caseno==case]
  for (comp in c("CI","CIV")){
    y = dat[[paste("raw",comp,sep="_")]][dat$caseno==case]
    dat[[paste("pd",comp,sep="_")]][dat$caseno==case] = -(rcoef[[comp]]$intercept+rcoef[[comp]]$slope*x-y)/sqrt(1+rcoef[[comp]]$slope^2)
  }
}

dat$sno = dat$caseno

dat$cellid = paste(dat$sno,sprintf("%04d",dat$Fibre),dat$Batch,sep="_")
# Discard repeated occurences of cellid
dat = dat[match(unique(dat$cellid),dat$cellid),]

axrange = c(-15,5)
linesat = c(-6,-4.5,-3,0,3,4.5,6)
ticksat = c(-20,-15,-10,-6,-4.5,-3,0,3,4.5,6,10,15,20)
thetarange = c(-30,90)
tickstheta = seq(thetarange[1],thetarange[2],10)
pdrange = c(-1.2,0.4)
tickspd = round(seq(pdrange[1],pdrange[2],0.2),2)
ratrange = range(c(dat$rat_CI,dat$rat_CIV),na.rm=TRUE,finite=TRUE)
ratrange = c(0.5,1.25)
ticksrat = round(seq(ratrange[1],ratrange[2],0.2),2)
spotcol = rgb(0,0,0,0.35)
ctcol = rgb(1,0,0,0.05)
mutcols = colorRamp(c(rgb(24/255,134/255,151/255),rgb(1,1,1),rgb(145/255,67/255,24/255)),space="Lab")
gcols=tetradic("purple",plot=FALSE)

op=par(mfrow=c(1,2))
for(ycol in c("adj_CI","adj_CIV")){
  plot(dat$adj_porin,dat[[ycol]],col=spotcol,pch=16,cex=0.2,xlab="raw_porin",ylab=ycol)
}
par(op)

tasnimcols = c("red","green","blue","purple","orange")
dat$col="yellow"


ctrl = dat[dat$pc=="c",]

NREPS = 1000
threshes=data.frame(row=1:NREPS)
threshdenses=data.frame(row=1:NREPS)

hiclust = FALSE

#looknums = sprintf("P%02d",1:17)
#looknums = looknums[looknums%in%dat$sno]

looknums = sort(unique(dat$caseno))
looknums = looknums[!grepl("C",looknums)]

unique(dat$caseno)

# Discard zeros
dat = dat[(dat$raw_CIV!=0)&(dat$raw_porin!=0)&(dat$raw_CI!=0),]
dat = dat[!is.na(dat$sno),]

mat = log(dat[,c("raw_CIV","raw_porin","raw_CI")])
colnames(mat) = gsub("adj","raw",colnames(mat))
rownames(mat) = dat$cellid

snos = sort(unique(dat$sno))
snocols = rainbow(length(snos),alpha=0.1)
names(snocols)=snos

write.table(dat,"IFdata.txt",sep="\t",quote=FALSE,row.names=FALSE)

dat$fibre_het[dat$sno=="P06"] = NA
dat = dat[!dat$sno%in%c("D11","D19"),]

gmms = list()

makemini = function(dt,usecols = c("adj_CI","adj_CIV","adj_porin")){
  mat = dt[,c("adj_CIV","adj_porin","adj_CI")]
  minimat = data.frame(mat)[,usecols]
  #colnames(minimat) = gsub("adj","raw",colnames(minimat))
  return(minimat)
}

cramp = colorRamp(c("red","white","blue"),space="Lab",interpolate="linear")
cfunc = function(x, alpha=1.0) {arr = cramp(x)/255; return(rgb(arr[1],arr[2],arr[3],alpha))}

samplerows = function(df,N){
  return(df[sample(nrow(df),N,replace=TRUE),])
}

diffdfs = function(df1,df2){ 
  return(as.numeric(sqrt(rowSums((df1-df2)^2))))
}

dts=list()
normcl=list()
pdf("Report3.pdf",width=18,height=18,pointsize=20)
patplots=list()
for(snum in looknums){
 patplots[[snum]]=list()
 s = as.character(snum)
 dt = dat[(dat$sno==s),]
 dt$biochem = "Normal"

 ct = data.frame(dat)[(dat$controls=="control"),]

 minimat = makemini(dt)
 mb = Mclust(minimat,1:2)
 if(max(mb$BIC[1,],na.rm=TRUE)>0.95*max(mb$BIC[2,],na.rm=TRUE)){ # If two clusters are not much better than one cluster, just use one:
   mb = Mclust(minimat,1)
 }

 ctvals = ct[,usecols]
 patvals = minimat[mb$z[,1]>0.5,]
 colnames(ctvals) = colnames(patvals)

 Nsamps = 15000
 ctsA = samplerows(ctvals,Nsamps)
 ctsB = samplerows(ctvals,Nsamps)
 z1 = samplerows(patvals,Nsamps)
 normclust = 1
 if(mb$G==2){
   #print((sum(mb$z[,2]>0.5)/length(mb$z[,2])))
   z2 = samplerows(minimat[mb$z[,2]>0.5,],Nsamps)
   if((sum(diffdfs(ctsA,z1)>diffdfs(ctsA,z2))/Nsamps)>0.5){normclust = 2} # Cluster 2 is closest to controls, therefore normal
   if((sum(mb$z[,1]>0.5)/length(mb$z[,1]))<0.1){normclust = 2} # Cluster 1 is less than 10% of total fibres, therefore deficient
   if((sum(mb$z[,2]>0.5)/length(mb$z[,2]))<0.1){normclust = 1} # Cluster 2 is less than 10% of total fibres, therefore deficient
   ctrlctrl = diffdfs(ctsA,ctsB)
   #res = hotelling.test(ctvals,patvals,perm=TRUE)
   #print(res)
  }
 if(mb$G==1){
   ctrlpat = diffdfs(ctsA,z1)
   if(mean(ctrlpat)>0.5){normclust = -1}
  }

 normcl[[snum]] = normclust

 if(normclust == 1) {
    predgmm = function(x) predict(mb,x)
 }
 if(normclust == 2) {
   predgmm = function(x){
    orig = predict(mb,x)
    lup = c(1,2)
    names(lup) = c(2,1)
    orig$classification = lup[orig$classification]
    orig$z[,1:2] = orig$z[,2:1]
    return(orig)   
   }
 }
 if(normclust == -1) {
   predgmm = function(x){
    orig = predict(mb,x)
    orig$z = 1.0 - orig$z
    return(orig)
   }
 }
 
  gmms[[snum]] = predgmm
  probclass = predgmm(minimat)$z
  probs = probclass[,1]

  dt$colours = sapply(probs,cfunc,alpha=0.5)
  dt$clust = 0
  pcutoff = 0.6
  dt$clust[probs>pcutoff] = 1
  dt$clust[probs<(1-pcutoff)] = 2
  dt$probs = probs

  dt$colours[dt$clust==0] = rgb(75/255,0/255,130/255,0.5)
  dt$biochem[dt$clust==0] = "Unknown"
  dt$biochem[dt$clust==1] = "Normal"
  dt$biochem[dt$clust==2] = "Deficient"

  if(mb$G==1) {
   dt$colours = rgb(75/255,0/255,130/255,0.5)
   dt$biochem = "Unknown"
   dt$clust = 0
  }

  dt$biochem = factor(dt$biochem,levels=c("Unknown","Normal","Deficient"))

 sampled = dt[!is.na(fibre_het),]
 unsampled = dt[is.na(fibre_het),]

 nsamp = length(sampled$fibre_het)
 homogenate = mean(dt$hom_het,na.rm=TRUE)

 deficient = sum(dt$biochem=="Deficient")
 sampdeficient = sum(dt$biochem=="Deficient"&!is.na(dt$fibre_het))
 samp = dim(sampled)[1]
 allfibs = dim(dt)[1]
 randest = samplesneeded(sampdeficient,allfibs,deficient)
 print(paste("From",snum,"Tasnim observed",allfibs,"fibres,",deficient,"of them were deficient.  Tasnim sampled",samp,"fibres for genotyping, of which",sampdeficient,"were deficient.  If sampling were random, she would have had to sample",randest,"to be likely to get that many (or more) with probability > 0.5."))
 op=par(mfrow=c(2,2))

pil = list()
 for(ycol in c("adj_CI","adj_CIV")){
   comp = gsub("adj_","",ycol)
   x = dt$adj_porin
   y = dt[[paste("adj",comp,sep="_")]]
   cx = ct$adj_porin
   cy = ct[[paste("adj",comp,sep="_")]]
  rx = c(-0.75,1.75)
  ry = rx

  ccols = c(rgb(0,0,1,0.55),rgb(1,0,0,0.55),rgb(75/255,0/255,130/255,0.55))

  FWT = signif(100*sum(dt$biochem=="Normal")/length(dt$biochem),3)
  FD = signif(100*sum(dt$biochem=="Deficient")/length(dt$biochem),3)
  mlab1=paste(snum,"Sampled:",nsamp,"Homogenate mutation level:",homogenate,"(%)")
  mlab2=paste("Total:",length(dt$biochem),"Normal:",FWT,"(%) Deficient:",FD,"(%)")
  mfull = paste(mlab1,mlab2,sep="\n")

  mod = lm(cy~cx)
  synx = seq(rx[1],rx[2],length.out=50)
  pred = predict(mod,newdata=data.frame(cx=synx),se.fit=TRUE, interval="prediction",na.action=na.omit,level=0.95)$fit
  mid = pred[,1]
  up = pred[,3]
  low = pred[,2]
  pil[[comp]]=list(low=low,mid=mid,up=up)
  
  patplots[[snum]][[ycol]] = function(){
   plot(NULL,xlab="log(Porin)",ylab=paste("log(",substring(ycol,5,nchar(ycol)),")",sep=""),main=mfull,xlim=rx,ylim=ry,cex.axis=1.55,cex.lab=1.55)
   points(cx,cy,pch=16,cex=0.5,col=rgb(0,0,0,0.15))
   points(dt$adj_porin,dt[[ycol]],pch=16,col=dt$colours,border=NULL)
   points(sampled$adj_porin,sampled[[ycol]],col=sampled$col,pch=15+as.numeric(sampled$fibre_type),cex=0.5)
   #abline(a=rcoef[[comp]]$intercept,b=rcoef[[comp]]$slope,lwd=2,col="red",lty=2)
   #abline(mod,lwd=2,col="grey")
   #points(synx,up,type="l",lty=2,col="grey")
   #points(synx,low,type="l",lty=2,col="grey")
  }
  patplots[[snum]][[ycol]]()
 }

   comp = "CIV"
   x = dt$adj_CI
   y = dt[[paste("adj",comp,sep="_")]]
   cx = ct$adj_CI
   cy = ct[[paste("adj",comp,sep="_")]]
  rx = c(-0.75,1.75)
  ry = rx

  ccols = c(rgb(0,0,1,0.55),rgb(1,0,0,0.55),rgb(75/255,0/255,130/255,0.55))

  FWT = signif(100*sum(dt$biochem=="Normal")/length(dt$biochem),3)
  FD = signif(100*sum(dt$biochem=="Deficient")/length(dt$biochem),3)
  mlab1=paste(snum,"Sampled:",nsamp,"Homogenate mutation level:",homogenate,"(%)")
  mlab2=paste("Total:",length(dt$biochem),"Normal:",FWT,"(%) Deficient:",FD,"(%)")
  mfull = paste(mlab1,mlab2,sep="\n")
  mod = lm(cy~cx)
  synx = seq(rx[1],rx[2],length.out=50)
  pred = predict(mod,newdata=data.frame(cx=synx),se.fit=TRUE, interval="prediction",na.action=na.omit,level=0.95)$fit
  mid = pred[,1]
  up = pred[,3]
  low = pred[,2]   
  pil[[comp]]=list(low=low,mid=mid,up=up)

  patplots[[snum]][["CI_CIV"]] = function(){
   plot(NULL,xlab="log(CI)",ylab=paste("log(",substring(ycol,5,nchar(ycol)),")",sep=""),main=mfull,xlim=rx,ylim=ry,cex.axis=1.55,cex.lab=1.55)
   points(cx,cy,pch=16,cex=0.5,col=rgb(0,0,0,0.15))
   points(dt$adj_CI,dt[[ycol]],pch=16,col=dt$colours,border=NULL)
   points(sampled$adj_CI,sampled[[ycol]],col=sampled$col,pch=15+as.numeric(sampled$fibre_type),cex=0.5)
   #abline(a=rcoef[[comp]]$intercept,b=rcoef[[comp]]$slope,lwd=2,col="red",lty=2)
   #abline(mod,lwd=2,col="grey")
   #points(synx,up,type="l",lty=2,col="grey")
   #points(synx,low,type="l",lty=2,col="grey")
  }
  patplots[[snum]][["CI_CIV"]]()

   getThreshes = function(sampled,c1=0.5,c2=0.5){
     if(length(unique(sampled$clust))==1){
        return(c(NA,NA))
     }
     
     cl1 = bootsamp(sampled$fibre_het[sampled$clust==1])
     cl2 = bootsamp(sampled$fibre_het[sampled$clust==2])

     d2 = density(cl2,n=1001,from=0,to=100,bw=1)
     d1 = density(cl1,n=1001,from=0,to=100,bw=1)

     if(length(unique(cl2))>1){
       b2 = fitdist(cl2/100, "beta")
       b2dens = function(x) dbeta(x,b2$estimate[["shape1"]],b2$estimate[["shape2"]])
     }else{
       b2dens = function(x) return(0)
     }

     if(length(unique(cl1))>1){
       b1 = fitdist(cl1/100, "beta")
       b1dens = function(x) dbeta(x,b1$estimate[["shape1"]],b1$estimate[["shape2"]])
     }else{
       b1dens = function(x) return(0)
     }

     props = seq(0.3,0.99,length.out=1001)
     diffs = b1dens(props) - b2dens(props)
     
     fA = approxfun(d1$x,d1$y/(sum(d1$y)*(d1$x[2]-d1$x[1]))*c1)
     fB = approxfun(d2$x,d2$y/(sum(d2$y)*(d2$x[2]-d2$x[1]))*c2)

     diff = function(x) fA(x) - fB(x)
     roots = uniroot.all(diff,c(0,100))
     denses = fA(roots)
     thresh = roots[which.max(denses)]
     threshdens = denses[which.max(denses)]
     if(length(roots)==0){
        thresh = NA
        threshdens = NA
     }
     return(c(thresh,threshdens))
   }

  bootsamp = function(x) sample(x,length(x),replace=TRUE)
  bootsampled = function(sampled){
   inds = 1:dim(sampled)[1]
   bsinds = sample(inds,length(inds),replace=TRUE)
   return(sampled[bsinds,])
  }

  getThreshes = function(sampled,c1,c2){
     res=c(Inf,Inf)
     while(!is.finite(res[1])){
      samp = bootsampled(sampled)
      res = naive_thresh(samp,c1,c2)
     }
     return(c(res,res))
  }

  
  # counts
  count1 = length(dt$clust[(dt$clust==1)])
  count2 = length(dt$clust[(dt$clust==2)])
  countall = length(dt$clust)

 if(length(sampled$caseno)>0 ){
   if((sum(sampled$biochem=="Deficient")>10)&(length(unique(sampled$biochem))>1)){
      threshmat = replicate(NREPS,getThreshes(sampled,count1/countall,count2/countall))
     }else{
      threshmat = matrix(NA,2,NREPS)
     }

   threshes[snum]=threshmat[1,]
   threshdenses[snum]=threshmat[2,]
   patplots[[snum]][["mutationlevel"]] = function(){
    boxplot(fibre_het~biochem,data=sampled,notch=FALSE,outline=FALSE,pars=list(boxwex=0.4),ylab="Mutation level (%)",ylim=c(0,100),xlab="RC status",cex.axis=1.55,cex.lab=1.55)
    x2 = 3+runif(sum(sampled$clust==2),-0.15,0.15)
    x1 = 2+runif(sum(sampled$clust==1),-0.15,0.15)
    x0 = 1+runif(sum(sampled$clust==0),-0.15,0.15)
    points(x2,sampled$fibre_het[sampled$clust==2],pch=16,col=ifelse(length(unique(sampled$clust))>1,ccols[2],ccols[3]))
    points(x1,sampled$fibre_het[sampled$clust==1],pch=16,col=ifelse(length(unique(sampled$clust))>1,ccols[1],ccols[3]))
    points(x0,sampled$fibre_het[sampled$clust==0],pch=16,col=ccols[3])
    points(x2,sampled$fibre_het[sampled$clust==2],pch=15+as.numeric(sampled$fibre_type[sampled$clust==2]),col=sampled$col[sampled$clust==2],cex=0.45)
    points(x1,sampled$fibre_het[sampled$clust==1],pch=15+as.numeric(sampled$fibre_type[sampled$clust==1]),col=sampled$col[sampled$clust==1],cex=0.45)
    points(x0,sampled$fibre_het[sampled$clust==0],pch=15+as.numeric(sampled$fibre_type[sampled$clust==0]),col=sampled$col[sampled$clust==0],cex=0.45)
    abline(h= median(threshmat[1,],na.rm=TRUE),col="grey",lwd=2,lty=2)
    abline(h=homogenate,col="grey",lwd=2,lty=3)
   }
 }else{
  patplots[[snum]][["mutationlevel"]] = function(){plot.new()}
 }
 patplots[[snum]][["mutationlevel"]]()

par(op)

layout(matrix(c(1,1,2,2,3,3,
                4,4,4,5,5,5),nrow=2,byrow=TRUE))

for (colname in c("adj_CI","adj_CIV","adj_porin")){
   x = dt[[colname]]
  cx = ct[[colname]]

  densx = density(x)
  denscx = density(cx)
  rangexax = range(c(densx$x,denscx$x))
  rangeyax = range(c(densx$y,denscx$y))
  patplots[[snum]][[paste("dens",colname)]] = function(){
  plot(denscx,xlim=rangexax,ylim=rangeyax,lwd=2,main=snum,lty=2,xlab=paste("log(",substring(colname,5,nchar(colname)),")",sep=""),cex.axis=1.55,cex.lab=1.55,cex.main=1.55)
  points(densx,type="l",lwd=2)
  points(dt[[colname]],rep(0,length(dt[[colname]])),pch=3,col=dt$colours)
  legend("topleft",c("Control","Patient"),lwd=2,lty=c(2,1))
  }
  patplots[[snum]][[paste("dens",colname)]]()
}

 if(length(sampled$caseno)>0&sum(sampled$biochem=="Deficient")>10&length(unique(as.character(sampled$biochem))>1)){

   nt = naive_thresh(sampled,count1/countall,count2/countall,makeplot=TRUE)
   patplots[[snum]][["kde"]] = nt$kde
   patplots[[snum]][["naive"]] = nt$naive
 }else{
   patplots[[snum]][["kde"]] = function() {plot.new()}
   patplots[[snum]][["naive"]] = function() {plot.new()}
 }
 patplots[[snum]][["kde"]]() 
 patplots[[snum]][["naive"]]() 
 par(mfrow=c(1,1))

 mat = log(dt[,c("raw_CIV","raw_porin","raw_CI")])
 colnames(mat) = gsub("adj","raw",colnames(mat))
 rownames(mat) = dt$cellid

 dts[[s]] = dt


}

dts = do.call(rbind,dts)

threshes$row=NULL
threshdenses$row=NULL

threshes=threshes[,sort(colnames(threshes))]
threshdenses=threshdenses[,sort(colnames(threshdenses))]

par(op)

drops=list(c(),c("P02", "P05", "P10", "P14"),c("P02", "P05", "P06", "P10", "P14", "P18"))
d = drops[3]
threshsm = threshes[,!colnames(threshes)%in%d[[1]]]
threshsm[["All"]] = apply(threshsm,1,median,na.rm=TRUE)
threshplot = function(){
 op = par(mai=c(1.5,1.75,0.5,0.5))
 stripchart(threshsm,vertical=TRUE,method="jitter",jitter=0.2,pch=16,col=rgb(1,0,0,0.05),ylab="Mutation level threshold (%)",ylim=c(0,100),main="",cex.axis=1.75,cex.lab=2.55)
 beanplot(threshsm,add=TRUE,outline=FALSE,pars=list(boxwex=0.5),cex.axis=1.25,cex.lab=2.55,col=rgb(0,0,0,0),ll=0,wd=1.5,beanlinewd=0.5,beanlines="median",overallline="median",maxwidth=1.0,axes=FALSE)
 par(op)
}
threshplot()

dev.off()
print(summary(threshsm))
apply(threshsm,2,IQR,na.rm=TRUE)

png("OverallThreshold.png",width=1000,height=1000,pointsize=18)
 threshplot()
dev.off()

bthresh = function(threshes,A,B){
 N = dim(threshes)[1]
 x = c(threshes[[A]],threshes[[B]])
 x = sample(x,length(x),replace=TRUE)
 a = x[1:N]; b = x[(N+1):(2*N)]
 return(mean(a-b,na.rm=TRUE))
}

bprop = function(threshes,A,B,N){
  a = sample(threshes[[A]],N,replace=TRUE)
  b = sample(threshes[[B]],N,replace=TRUE)
  return(sum(a>b,na.rm=TRUE)/N)
}

permutep = function(dat,A,B){
  countA1 = length(dat$clust[(dat$sno==A)&(dat$clust==1)])
  countA2 = length(dat$clust[(dat$sno==A)&(dat$clust==2)])
  countA = length(dat$clust[(dat$sno==A)])
  countB1 = length(dat$clust[(dat$sno==B)&(dat$clust==1)])
  countB2 = length(dat$clust[(dat$sno==B)&(dat$clust==2)])
  countB = length(dat$clust[(dat$sno==B)])
  # Select only the rows from the two patients A and B
  datAB = dat[(dat$sno%in%c(A,B))&(!is.na(dat$fibre_het)),]
  # Shuffle the patient labels & split the data by shuffled labels
  datAB$sno = sample(datAB$sno,replace=FALSE)
  dtAt = datAB[datAB$sno==A,]
  dtBt = datAB[datAB$sno==B,]
  
  At = threshFromScratch(dtAt,gmms[[A]],countA1/(countA),countA2/(countA))
  Bt = threshFromScratch(dtBt,gmms[[B]],countB1/(countB),countB2/(countB))
  return(Bt-At)
}

getSig3 = function(dat,A,B,N=5000,makeplot=FALSE){
  countA1 = length(dat$clust[(dat$sno==A)&(dat$clust==1)])
  countA2 = length(dat$clust[(dat$sno==A)&(dat$clust==2)])
  countA = length(dat$clust[(dat$sno==A)])
  countB1 = length(dat$clust[(dat$sno==B)&(dat$clust==1)])
  countB2 = length(dat$clust[(dat$sno==B)&(dat$clust==2)]) 
  countB = length(dat$clust[(dat$sno==B)])
  resA = threshFromScratch(dat[(dat$sno==A)&(!is.na(dat$fibre_het)),],gmms[[A]],countA1/(countA),countA2/(countA))
  resB = threshFromScratch(dat[(dat$sno==B)&(!is.na(dat$fibre_het)),],gmms[[B]],countB1/(countB),countB2/(countB))
  est = resB-resA
  teststat = replicate(N,permutep(dat,A,B))
  teststat = teststat[!is.na(teststat)]
  pval2 = sum(abs(teststat)>abs(est))/length(teststat)
  if(makeplot){
    plot(density(abs(teststat)),main=paste("p-value:",formatC(pval2,5)),xlab=paste("|",B,"-",A,"|"),lwd=3,xlim=c(0,100),cex.lab=2.55,cex.axis=1.55,cex.main=2)
    abline(v=abs(est),col="green",lwd=2,lty=2)
    legend("topright",legend=c("H0 (no difference)","Observed"),col=c("black","green"),lwd=c(3,2),lty=c(1,2),cex=1.5)
  }
  return(pval2)
}

pdf("TestStatistic.pdf")
#png("TestStat.png",width=1000,height=1000,pointsize=18)
op = par(mai=c(1.1,1.1,0.9,0.1))
combs = combn(colnames(threshsm)[colnames(threshsm)!="All"],2)
ncombs = dim(combs)[2]
sigs = replicate(ncombs,NA)
diffs = replicate(ncombs,NA)
probs_diff = replicate(ncombs,NA)
for(i in 1:ncombs){
  print(paste(i,"of",ncombs,combs[1,i],combs[2,i]))
  sigs[i] = tryCatch(getSig3(dts,combs[1,i],combs[2,i],N=10000,makeplot=TRUE),error=function(err) NaN)
  diffs[i] = median(threshes[[combs[2,i]]]-threshes[[combs[1,i]]],na.rm=TRUE)
  probs_diff[i] = ifelse(diffs[i]>0.0,bprop(threshes,combs[1,i],combs[2,i],1000),bprop(threshes,combs[2,i],combs[1,i],1000))
}
par(op)
dev.off()

sigs2 = sigs
sigs2[is.na(diffs)] = NaN

sigsadj = p.adjust(sigs2,method="hochberg")

pdf("Mutloads.pdf",width=18,height=18,pointsize=20)
dml = dat[!is.na(dat$fibre_het),]
 op = par(mai=c(1.5,1.75,0.5,0.5))
 stripchart(dml$fibre_het~dml$caseno,vertical=TRUE,method="jitter",jitter=0.15,pch=16,col=rgb(1,0,0,0.25),ylab="Mutation level (%)",ylim=c(0,100),main="",cex.axis=1.75,cex.lab=2.55)
 beanplot(dml$fibre_het~dml$caseno,add=TRUE,outline=FALSE,pars=list(boxwex=1.5),cex.axis=1.25,cex.lab=2.55,col=rgb(0,0,0,0),ll=0,wd=3.5,beanlinewd=0.5,beanlines="median",overallline="median",maxwidth=1.0,axes=FALSE)
 par(op)
dev.off()

defs = aggregate(as.character(dts$biochem)=="Deficient",by=list(dts$caseno),sum)
fibs = aggregate(as.character(dts$biochem)=="Deficient",by=list(dts$caseno),length)
deffracs = data.frame(caseno=defs[,1],deficient=100*defs$x/fibs$x)
write.table(deffracs,"PropDef.txt",sep="\t",quote=FALSE,row.names=FALSE)

dat$fibre_type[dat$fibre_type=="2,3"]=2.5

outdat = data.frame(
 sampleID=dat$caseno,
 fibreID=dat$Fibre,
 sbj_type=ifelse(grepl("C",dat$caseno, fixed=TRUE),"control","patient"),
 Porin=dat$adj_porin,
 CI=dat$adj_CI,
 CIV=dat$adj_CIV,
 CopyNumber=as.numeric(dat$fibre_total_cn),
 FibreType=as.numeric(dat$fibre_type),
 FibreType2=as.numeric(dat$fibre_type_2),
 MutationLoad=as.numeric(dat$fibre_het)
)

longdat = melt(setDT(outdat),id.vars=c("sampleID","fibreID","sbj_type"),variable.name="channel")
write.table(longdat,"DataChildsFormat_adj.txt",sep="\t",quote=FALSE,row.names=FALSE)

outdat = data.frame(
 sampleID=dat$caseno,
 fibreID=dat$Fibre,
 sbj_type=ifelse(grepl("C",dat$caseno, fixed=TRUE),"control","patient"),
 Porin=dat$raw_porin,
 CI=dat$raw_CI,
 CIV=dat$raw_CIV,
 CopyNumber=as.numeric(dat$fibre_total_cn),
 FibreType=as.numeric(dat$fibre_type),
 FibreType2=as.numeric(dat$fibre_type_2),
 MutationLoad=as.numeric(dat$fibre_het)
)

longdat = melt(setDT(outdat),id.vars=c("sampleID","fibreID","sbj_type"),variable.name="channel")
write.table(longdat,"DataChildsFormat_raw.txt",sep="\t",quote=FALSE,row.names=FALSE)

