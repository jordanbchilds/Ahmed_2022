library(MASS)
source("helper_functions.R", local = TRUE)

folder = "linReg_classifier"

dir.create("PDF", showWarnings=FALSE)
dir.create(file.path("PDF",folder), showWarnings=FALSE)

cord = c("raw_CI", "raw_CIV")
mitochan = "porin"

dat = read.csv(file.path("..", "rawdat_a.csv"), stringsAsFactors=FALSE)
sbj = unique(dat$caseno)
crl = sbj[grepl("C0.", sbj)]
pts = sort(sbj[grepl("P0.", sbj)])

######################
### the plots
######################

pdf(file.path("PDF", folder, "MCMC.pdf"), width=13, height=8)
{ 
    for(chan in cord){
      outroot_ctrl = paste(chan, "CONTROL", sep="_")
      
      MCMCplot(folder, chan, lag=100, 
               title=paste(chan, "CONTROL"))
      for(pat in pts){
        outroot_pat = paste(froot, chan, pat, sep="_")
        MCMCplot(folder, chan, pat, lag=100,
                 title=paste(chan, pat))
      }
    }
}
dev.off()

pdf(file.path("PDF", folder, "model_post.pdf"), width=13, height=8)
{
    for(chan in cord){
      
    ctrl_data =  getData_mats(chan=chan)$ctrl
    pts_data = getData_mats(chan=chan)$pts
    xlims = range(c(ctrl_data[,1], pts_data[,1]))
    ylims = range(c(ctrl_data[,2], pts_data[,2]))
      
    # priorpost(ctrl_data=ctrl_data, 
    #           priorpred=output_reader0(out_type="PRIORPRED"),
    #           postpred=output_reader0(out_type="POSTPRED"),
    #           chan=chan,title=paste(chan, "CONTROL"),
    #           xlims=xlims, ylims=ylims)
    #   
    for(pat in pts){
      pat_data = getData_mats(chan=chan, pts=pat)$pts
        
      priorpost(ctrl_data=ctrl_data, pat_data=pat_data,
                classif=output_reader(folder, chan, pat=pat, out_type="CLASSIF")[[1]],
                priorpred=output_reader(folder, chan, pat=pat, "PRIORPRED"), 
                postpred=output_reader(folder, chan, pat=pat, "POSTPRED"),
                chan=chan, mitochan="VDAC1", title=paste(chan, pat),
                xlims=xlims, ylims=ylims)
      }
    }
}
dev.off()

pdf(file.path("PDF", folder, "pi_post.pdf"), width=13, height=8)
{
  par(mfrow=c(1,1))
  for(fulldat in fulldats_raw){
    for(chan in c("MTCO1", "NDUFB8")){
      pipost_plotter(fulldat, chan, folder=folder ) 
    }
  }
}
dev.off()

pdf(file.path("./PDF", folder, "classif_plot.pdf"), width=13, height=8)
{
  for(fulldat in fulldats_raw){
    chanpts = getData_chanpats(fulldat)
    
    for(chan in chanpts$channels){
      output_reader0 = function(pat="CONTROL", out_type){
        output_reader(folder=folder, fulldat=fulldat, chan=chan, pat, out_type)
      }
      ctrl_data =  getData_mats(fulldat, chan=chan)$ctrl
      pts_data = getData_mats(fulldat, chan=chan)$pts
      
      for(pat in chanpts$patients){
        pat_data = getData_mats(fulldat, chan=chan, pts=pat)$pts
        
        classif_plot(ctrl_data=ctrl_data, pat_data=pat_data,
                  classifs_pat=output_reader0(pat=pat, out_type="CLASSIF")[[1]],
                  chan=chan, mitochan="VDAC1", 
                  title=paste(chan, pat) )
      }
    }
  }
}
dev.off()

pdf(file.path("./PDF", folder, "pipost_freqComparison.pdf"), width=13, height=8)
{
  op = par(mfrow=c(1,1), mar=c(6,6,6,3), cex.main=2, cex.lab=2, cex.axis=1.5 )
  for(fulldat in fulldats_raw){
    froot = gsub("", ".RAW.txt", fulldat)
    for(fulldat in fulldats_raw){
      for(chan in c("MTCO1", "NDUFB8")){
        pipost_plotter(fulldat, chan, folder=folder, freq_compare=TRUE) 
      }
    }
  }
  par(op)
}
dev.off()







