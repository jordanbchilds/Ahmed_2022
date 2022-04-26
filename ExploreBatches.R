dat = fread("rawdat_a.csv",sep=",",stringsAsFactors=FALSE)
dat = dat[nchar(dat$caseno)==3,]

dat$cellid = paste(dat$caseno,sprintf("%04d",dat$Fibre),sep="_")
dat = dat[match(unique(dat$cellid),dat$cellid),]

print(unique(dat$caseno[dat$controls=="control"]))

looknums = sort(unique(dat$caseno))
print(looknums)
for(snum in looknums){
 dt = dat[(dat$caseno==snum),]
 print(snum)
 print(paste("Batch:",unique(dt$Batch)))
}

batches = sort(unique(dat$Batch))
print(batches)
for(batch in batches){
 db = dat[dat$Batch==batch,]
 print(paste("Batch:",batch))
 print(sort(unique(db$caseno)))
}
