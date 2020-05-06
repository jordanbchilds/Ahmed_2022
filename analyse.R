library(colortools)
library(data.table)
library(fitdistrplus)
library(gmp)
library(mclust)
library(umap)
library(Hotelling)

set.seed(42)

usecols = c("adj_CI","adj_CIV","adj_porin")

# https://towardsdatascience.com/unsupervised-machine-learning-clustering-analysis-d40f2b34ae7e

dat = fread("rawdat_a.csv",sep=",",stringsAsFactors=FALSE)
dat = dat[nchar(dat$caseno)==3,]

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

dat$cellid = paste(dat$sno,sprintf("%04d",dat$Fibre),sep="_")

ctrl = dat[dat$pc=="c",]

threshes=data.frame(row=1:1000)
threshdenses=data.frame(row=1:1000)

hiclust = FALSE

nums = unique(dat$caseno)
looknums = sort(nums[grepl("P",nums)])
looknums = looknums[looknums%in%dat$sno]

sort(unique(dat$caseno))

# Discard repeated occurences of controls
dat = dat[match(unique(dat$cellid),dat$cellid),]

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
dat = dat[!dat$sno%in%c("P11","P19"),]
looknums = looknums[!looknums%in%c("P11","P19")]

gmms = list()

makemini = function(dt,usecols = c("adj_CI","adj_CIV","adj_porin")){
  mat = dt[,c("adj_CIV","adj_porin","adj_CI")]
  minimat = data.frame(mat)[,usecols]
  colnames(minimat) = gsub("adj","raw",colnames(minimat))
  return(minimat)
}

pdf("Report3.pdf",width=18,height=18,pointsize=20)
for(snum in looknums){
 s = as.character(snum)
 dt = dat[(dat$sno==s),]
 dt$biochem = "Normal"
 ct = data.frame(dat)[(dat$controls=="control")&(dat$Batch%in%unique(dt$Batch)),]

 minimat = makemini(dt)
 mb = Mclust(minimat,1:2)

	 if(max(mb$BIC[1,],na.rm=TRUE)>0.95*max(mb$BIC[2,],na.rm=TRUE)){ # If two clusters are not much better than one cluster, just use one:
	   mb = Mclust(minimat,1)
	 }

      gmms[[snum]] = mb

    cramp = colorRamp(c("red","white","blue"),space="Lab",interpolate="linear")
	cfunc = function(x, alpha=1.0) {arr = cramp(x)/255; return(rgb(arr[1],arr[2],arr[3],alpha))}

	samplerows = function(df,N){
         return(df[sample(nrow(df),N,replace=TRUE),])
	}

      diffdfs = function(df1,df2){
         return(as.numeric(sqrt(rowSums((df1-df2)^2))))
	}
      ctvals = ct[,usecols]
      patvals = minimat[mb$z[,1]>0.5,]
      colnames(ctvals) = colnames(patvals)

      Nsamps = 150000
	ctsA = samplerows(ctvals,Nsamps)
      ctsB = samplerows(ctvals,Nsamps)
      z1 = samplerows(patvals,Nsamps)
      
	 if(mb$G==2){
         z2 = samplerows(minimat[mb$z[,2]>0.5,],Nsamps)
         if((sum(diffdfs(ctsA,z1)>diffdfs(ctsA,z2))/Nsamps)>0.5){probs=mb$z[,2]}else{probs=mb$z[,1]}
       }else{
         ctrlctrl = diffdfs(ctsA,ctsB)
         ctrlpat = diffdfs(ctsA,z1)
         
         print(snum)
	   print(mean(ctrlpat))
         if(mean(ctrlpat)>0.5){probs = 1.0-mb$z[,1]}else{probs = mb$z[,1]}

         res = hotelling.test(ctvals,patvals,perm=TRUE)
         print(res)
	}
	 dt$colours = sapply(probs,cfunc,alpha=0.5)
       if(mb$G==1) dt$colours = rgb(75/255,0/255,130/255,0.5)
	 dt$clust = 0
	 pcutoff = 0.6
	 dt$clust[probs>pcutoff] = 1
	 dt$clust[probs<(1-pcutoff)] = 2
	 dt$probs = probs

 dt$biochem[dt$clust==0] = "Unknown"
 dt$biochem[dt$clust==1] = "Normal"
 dt$biochem[dt$clust==2] = "Deficient"

 if(mb$G==1) {
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

 op=par(mfrow=c(2,2))

pil = list()
 for(ycol in c("adj_CI","adj_CIV")){
   comp = gsub("adj_","",ycol)
   x = dt$adj_porin
   y = dt[[paste("adj",comp,sep="_")]]
  cx = dat$adj_porin[(dat$controls=="control")&(dat$Batch%in%unique(dt$Batch))]
  cy = dat[[paste("adj",comp,sep="_")]][(dat$controls=="control")&(dat$Batch%in%unique(dt$Batch))]
  rx = c(-0.75,1.75)
  ry = rx

  ccols = c(rgb(0,0,1,0.55),rgb(1,0,0,0.55),rgb(75/255,0/255,130/255,0.55))

  FWT = signif(100*sum(dt$biochem=="Normal")/length(dt$biochem),3)
  FD = signif(100*sum(dt$biochem=="Deficient")/length(dt$biochem),3)
  mlab1=paste(snum,"Sampled:",nsamp,"Homogenate mutation load:",homogenate,"(%)")
  mlab2=paste("Total:",length(dt$biochem),"Normal:",FWT,"(%) Deficient:",FD,"(%)")
  mfull = paste(mlab1,mlab2,sep="\n")

  plot(NULL,xlab="log(Porin)",ylab=paste("log(",substring(ycol,5,nchar(ycol)),")",sep=""),main=mfull,xlim=rx,ylim=ry,cex.axis=1.55,cex.lab=1.55)
  points(cx,cy,pch=16,cex=0.5,col=rgb(0,0,0,0.15))
  points(dt$adj_porin,dt[[ycol]],pch=16,col=dt$colours,border=NULL)
  points(sampled$adj_porin,sampled[[ycol]],col=sampled$col,pch=16,cex=0.5)

  mod = lm(cy~cx)
  synx = seq(rx[1],rx[2],length.out=50)
  pred = predict(mod,newdata=data.frame(cx=synx),se.fit=TRUE, interval="prediction",na.action=na.omit,level=0.95)$fit
  mid = pred[,1]
  up = pred[,3]
  low = pred[,2]   
  pil[[comp]]=list(low=low,mid=mid,up=up)
 }

   comp = "CIV"
   x = dt$adj_CI
   y = dt[[paste("adj",comp,sep="_")]]
  cx = dat$adj_CI[(dat$controls=="control")&(dat$Batch%in%unique(dt$Batch))]
  cy = dat[[paste("adj",comp,sep="_")]][(dat$controls=="control")&(dat$Batch%in%unique(dt$Batch))]
  rx = c(-0.75,1.75)
  ry = rx

  ccols = c(rgb(0,0,1,0.55),rgb(1,0,0,0.55),rgb(75/255,0/255,130/255,0.55))

  FWT = signif(100*sum(dt$biochem=="Normal")/length(dt$biochem),3)
  FD = signif(100*sum(dt$biochem=="Deficient")/length(dt$biochem),3)
  mlab1=paste(snum,"Sampled:",nsamp,"Homogenate mutation load:",homogenate,"(%)")
  mlab2=paste("Total:",length(dt$biochem),"Normal:",FWT,"(%) Deficient:",FD,"(%)")
  mfull = paste(mlab1,mlab2,sep="\n")

  plot(NULL,xlab="log(CI)",ylab=paste("log(",substring(ycol,5,nchar(ycol)),")",sep=""),main=mfull,xlim=rx,ylim=ry,cex.axis=1.55,cex.lab=1.55)
  points(cx,cy,pch=16,cex=0.5,col=rgb(0,0,0,0.15))
  points(dt$adj_CI,dt[[ycol]],pch=16,col=dt$colours,border=NULL)
  points(sampled$adj_CI,sampled[[ycol]],col=sampled$col,pch=16,cex=0.5)
  mod = lm(cy~cx)
  synx = seq(rx[1],rx[2],length.out=50)
  pred = predict(mod,newdata=data.frame(cx=synx),se.fit=TRUE, interval="prediction",na.action=na.omit,level=0.95)$fit
  mid = pred[,1]
  up = pred[,3]
  low = pred[,2]   
  pil[[comp]]=list(low=low,mid=mid,up=up)

   getThreshes = function(sampled){
     if(length(unique(sampled$clust))==1){
        return(c(NA,NA))
     }
     
     cl1 = bootsamp(sampled$fibre_het[sampled$clust==1])
     cl2 = bootsamp(sampled$fibre_het[sampled$clust==2])

     d2 = density(cl2,n=1001,from=0,to=100,bw=5)
     d1 = density(cl1,n=1001,from=0,to=100,bw=5)

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


     thresh = d1$x[max(which(d1$y>d2$y))]
     threshdens = props[which.min(abs(diffs))]
     return(c(thresh,threshdens))
   }

  bootsamp = function(x) sample(x,length(x),replace=TRUE)

 if(length(sampled$caseno)>0 ){
   if((sum(sampled$biochem=="Deficient")>10)&(length(unique(sampled$biochem)>1))){
      threshmat = replicate(1000,getThreshes(sampled))
     }else{
      threshmat = matrix(NA,2,1000)
     }

   threshes[snum]=threshmat[1,]
   threshdenses[snum]=threshmat[2,]
   
   boxplot(fibre_het~biochem,data=sampled,notch=FALSE,outline=FALSE,pars=list(boxwex=0.4),ylab="Mutation load (%)",ylim=c(0,100),xlab="RC status",cex.axis=1.55,cex.lab=1.55)
   x2 = 3+runif(sum(sampled$clust==2),-0.15,0.15)
   x1 = 2+runif(sum(sampled$clust==1),-0.15,0.15)
   x0 = 1+runif(sum(sampled$clust==0),-0.15,0.15)
   points(x2,sampled$fibre_het[sampled$clust==2],pch=16,col=ifelse(length(unique(sampled$clust))>1,ccols[2],ccols[3]))
   points(x1,sampled$fibre_het[sampled$clust==1],pch=16,col=ifelse(length(unique(sampled$clust))>1,ccols[1],ccols[3]))
   points(x0,sampled$fibre_het[sampled$clust==0],pch=16,col=ccols[3])
   points(x2,sampled$fibre_het[sampled$clust==2],pch=16,col=sampled$col[sampled$clust==2],cex=0.35)
   points(x1,sampled$fibre_het[sampled$clust==1],pch=16,col=sampled$col[sampled$clust==1],cex=0.35)
   points(x0,sampled$fibre_het[sampled$clust==0],pch=16,col=sampled$col[sampled$clust==0],cex=0.35)
   abline(h= median(threshmat[1,]),col="grey",lwd=2,lty=2)
   abline(h=homogenate,col="grey",lwd=2,lty=3)
 }else{
  plot.new()
 }

for (colname in c("adj_CI","adj_CIV","adj_porin")){
   x = dt[[colname]]
  cx = dat[[colname]][(dat$controls=="control")&(dat$Batch%in%unique(dt$Batch))]

  densx = density(x)
  denscx = density(cx)
  rangexax = range(c(densx$x,denscx$x))
  rangeyax = range(c(densx$y,denscx$y))
  plot(denscx,xlim=rangexax,ylim=rangeyax,lwd=2,main=snum,lty=2,xlab=paste("log(",substring(colname,5,nchar(colname)),")",sep=""),cex.axis=1.55,cex.lab=1.55)
  points(densx,type="l",lwd=2)
  points(dt[[colname]],rep(0,length(dt[[colname]])),pch=3,col=dt$colours)
  legend("topleft",c("Control","Patient"),lwd=2,lty=c(2,1))
}

 if(length(sampled$caseno)>0){
   hist(sampled$fibre_het,seq(-2.5,102.5,5),xlim=c(0,100),xlab = "Mutation load (%)",main=snum,cex.axis=1.55,cex.lab=1.55)
   points(sampled$fibre_het,rep(0,length(sampled$fibre_het)),col=sampled$colours,pch=3,lwd=3)
   abline(v=mean(sampled$hom_het),col="grey",lty=3,lwd=3)
 }else{
  plot.new()
 }

 mat = log(dt[,c("raw_CIV","raw_porin","raw_CI")])
 colnames(mat) = gsub("adj","raw",colnames(mat))
 rownames(mat) = dt$cellid

}


threshes$row=NULL
threshdenses$row=NULL

threshes=threshes[,sort(colnames(threshes))]
threshdenses=threshdenses[,sort(colnames(threshdenses))]

drops=list(c(),c("P02", "P05", "P10", "P15"),c("P02", "P05", "P06", "P10", "P15", "P19"))
for(d in drops){
  stripchart(threshes[,!colnames(threshes)%in%d],vertical=TRUE,method="jitter",jitter=0.1,pch=16,col=rgb(1,0,0,0.1),ylab="Mutation load threshold (%)",ylim=c(0,100),main="",cex.axis=1.55,cex.lab=1.55)
  boxplot(threshes[,!colnames(threshes)%in%d],add=TRUE,outline=FALSE,pars=list(boxwex=0.3),cex.axis=1.55,cex.lab=1.55)
  print(summary(threshes))
}

 plot.new()

par(op)

dev.off()

getSig2 = function(threshes,A,B,showplot=FALSE){
  teststat = abs((threshes[[A]]-mean(threshes[[A]])) - (threshes[[B]] - mean(threshes[[B]])))
  est = abs(median(threshes[[A]]-threshes[[B]]))
  pval = abs(sum(teststat>est)/length(teststat))
  if(showplot&(sum(is.na(teststat))==0)){
    plot(density(teststat),main=pval)
    abline(v=est,col="red")
  }
  return(pval)
}

combs = combn(colnames(threshes),2)
ncombs = dim(combs)[2]
sigs = replicate(ncombs,NA)
probhigh = replicate(ncombs,NA)
problow = replicate(ncombs,NA)
for(i in 1:ncombs){
  print(paste(i,"of",ncombs))
  sigs[i] = getSig2(threshes,combs[1,i],combs[2,i])
  probhigh[i] = mean(threshes[[combs[2,i]]]-threshes[[combs[1,i]]],na.rm=TRUE)
}

sigs2 = sigs
sigs2[is.na(probhigh)] = NaN

sigsadj = p.adjust(sigs2,method="fdr")

pats = sort(looknums)
prep = matrix(data="",nrow = length(pats), ncol = length(pats))#
stren = matrix(data="",nrow = length(pats), ncol = length(pats))#
rownames(prep)=pats
colnames(prep)=pats
rownames(stren)=pats
colnames(stren)=pats
for(i in seq_along(combs[1,])){
  cmb = combs[,i]
  prep[cmb[1],cmb[2]] = ifelse(is.na(probhigh[i])|is.na(sigs2[i]),"",sprintf("%0.3f",probhigh[i]))
  prep[cmb[2],cmb[1]] = ifelse(is.na(probhigh[i])|is.na(sigs2[i]),"",sprintf("%0.3f",sigs2[i]))
}
print(prep)

write.table(prep,file="probabilities_raw_p.txt",sep="\t",row.names=TRUE,col.names=NA,quote=TRUE)

pats = sort(looknums)
prep = matrix(data="",nrow = length(pats), ncol = length(pats))#
stren = matrix(data="",nrow = length(pats), ncol = length(pats))#
rownames(prep)=pats
colnames(prep)=pats
rownames(stren)=pats
colnames(stren)=pats
for(i in seq_along(combs[1,])){
  cmb = combs[,i]
  prep[cmb[1],cmb[2]] = ifelse(is.na(probhigh[i])|is.na(sigsadj[i]),"",sprintf("%0.3f",probhigh[i]))
  prep[cmb[2],cmb[1]] = ifelse(is.na(probhigh[i])|is.na(sigsadj[i]),"",sprintf("%0.3f",sigsadj[i]))
}
print(prep)

write.table(prep,file="probabilities_adjusted_p.txt",sep="\t",row.names=TRUE,col.names=NA,quote=TRUE)



