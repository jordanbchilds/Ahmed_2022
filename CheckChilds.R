raw = read.delim("DataChildsFormat_raw.txt",sep="\t",stringsAsFactors=FALSE)
adj = read.delim("DataChildsFormat_adj.txt",sep="\t",stringsAsFactors=FALSE)

dat = adj

pats = sort(unique(dat$sampleID[grepl("P",dat$sampleID)]))

transform = identity
xlim = quantile(transform(dat$value[dat$channel=="Porin"]),c(0.0,1.0),na.rm=TRUE)

cairo_pdf("ChildsPlots_adj.pdf",width=14,height=7,onefile=TRUE)

for(pat in pats){

 cx = transform(dat$value[(dat$channel=="Porin")&(dat$sbj_type=="control")])
 px = transform(dat$value[(dat$channel=="Porin")&(dat$sampleID==pat)])

 op = par(mfrow=c(1,2))
 for(prot in c("CI","CIV")){
   if(identical(transform,log)){
     xlab = "log(Porin)"
     ylab = paste0("log(",prot,")")
   }else{
     xlab = "Porin"
     ylab = prot
   }
   ylim = quantile(transform(dat$value[dat$channel==prot]),c(0.0,1.0),na.rm=TRUE)
   cy = transform(dat$value[(dat$channel==prot)&(dat$sbj_type=="control")])
   py = transform(dat$value[(dat$channel==prot)&(dat$sampleID==pat)])
   plot(NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=pat,cex.lab=1.55,cex.main=2.5)
   points(cx,cy,pch=16,cex=0.5,col=rgb(0,0,0,0.1))
   points(px,py,pch=16,cex=0.75,col=rgb(1,0,1,0.25))
 }
 par(op)
}

dev.off()