library("devtools")
devtools::install_github("jordanbchilds/analysis2Dmito", force=TRUE)
library("analysis2Dmito")

library("parallel")

dir.create("PLOTS", showWarnings = FALSE)
dir.create("OUTPUT", showWarnings = FALSE)

exampleData = get_exampleData()

exampleData = as.data.frame(exampleData)

colnames(exampleData) = c("sampleID", "fibreID", "subject_type", "channel", "value")

exampleData$sampleID_fibreID = paste(exampleData$sampleID, exampleData$fibreID, sep="_")
problemIDs = exampleData[exampleData$value<=0, "sampleID_fibreID"]
exampleData = exampleData[ !(exampleData$sampleID_fibreID %in% problemIDs), ]

sbj = unique(exampleData$sampleID)
ctrlID = grep("C", sbj, value=TRUE)
pts = sort( sbj[!(sbj %in% ctrlID)] )

# measure of mass (x-axis in 2Dmito)
mitochan = "raw_porin"
# proteins of interest (y-axis in 2Dmito)
channels_all = unique(exampleData$channel)
channels = channels_all[channels_all != mitochan]

# pdf("./PLOTS/raw_data.pdf", width=9, height=6)
# {
#   op = par(mfrow=c(1,1), cex.main=2, cex.lab=2, cex.axis=1.5, mar=c(6,6,6,3))
#   for( chan in channels ){
#     for( pat in pts ){
#       data = getData_mats(exampleData, pts=pat, channels=c(mitochan, chan), ctrlID=ctrlID,
#                           getIndex=TRUE)
#       xlims = range(c(data$ctrl[,1], data$pts[,1]))
#       ylims = range(c(data$ctrl[,2], data$pts[,2]))
#       
#       plot(NA, xlim=xlims, ylim=ylims, 
#            xlab="porin", ylab=gsub("raw_", "", chan), main=pat)
#       points(data$ctrl, col=alphaBlack(0.05), pch=20)
#       points(data$pts, col=alphaGreen(0.1), pch=20)
#     }
#   }
#   par(op)
# }
# dev.off()


data_list = list()
for(chan in channels){
  for( pat in pts ){
    data = getData_mats(exampleData, pts=pat, channels=c(mitochan, chan), 
                        ctrlID=ctrlID, getIndex=TRUE)
    data_list[[paste(chan, pat, sep="__")]] = data

  }
}

ncores = detectCores() - 1
cl  = makeCluster(ncores)
{
  linreg_output = parLapply(cl, data_list, inference, MCMCout=10000)
}
stopCluster(cl)

for(root in names(linreg_output)){
  list_saver(linreg_output[[root]], root=file.path("Output", toupper(root)))
}

