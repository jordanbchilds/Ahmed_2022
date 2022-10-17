library("MASS")
source("./helper_functions.R", local = TRUE)

folder = "allData_naiveBayes_GMM_hier"

dir.create(file.path("./PDF"), showWarnings = FALSE)
dir.create(file.path("./PDF", folder), showWarnings = FALSE)

fulldats_raw = c(
  "IMV.E02.P01.PMT.M3243AG.QIF.7TR.RAW.txt",
  "IMV.E01.P02.PMT.M3243AG.QIF.TGo.RAW.txt",
  "IMV.E02.P01.PMT.M3243AG.QIF.TGo.RAW.txt",
  "IMV.E02.P01.PMT.M3243AG.QIF.TCF.RAW.txt"
)

prior = read.table(file.path("./Output", folder, "PRIOR.txt"), 
                   header=TRUE, stringsAsFactors=FALSE)

pdf(file.path("./PDF", folder, "MCMC.pdf"), width=13, height=8)
{ 
  for(fulldat in fulldats_raw){
    chanpts = getData_chanpats(fulldat)
    pts = gsub("_", ".", chanpts$patients)
      for(chan in chanpts$channels){
        MCMCplot(folder, fulldat, chan,  lag=100,
                  title=paste(fulldat, "\n", chan))
      }
  }
}
dev.off()

pdf(file.path("./PDF", folder, "model_post_RAW.pdf"), width=13, height=8)
{
  for(fulldat in fulldats_raw){
    chanpts = getData_chanpats(fulldat)
    for(chan in chanpts$channels){
      ctrl_data =  getData_mats(fulldat, chan=chan, ctrl_only=TRUE, 
                                data_transform=myData_transform)

      pat_data = getData_mats(fulldat, chan=chan,
                              data_transform=myData_transform)$pts
      
      classif_allpats = output_reader(folder, fulldat, chan, out_type="CLASS")[[1]]
      post = output_reader(folder, fulldat, chan, out_type="POST")
        
      priorpost(ctrl_data=ctrl_data, pat_data=pat_data,
                  classif=classif_allpats,
                  prior=prior, post=post,
                  chan=chan, mitochan="VDAC1", 
                  title=paste(gsub(".RAW.txt", "", fulldat), "\n", chan))
    }
  }
}
dev.off()

pdf(file.path("./PDF", folder, "model_post.pdf"), width=13, height=8)
{
  for(fulldat in fulldats_raw){
    chanpts = getData_chanpats(fulldat)
    for(chan in chanpts$channels){
      ctrl_data =  getData_mats(fulldat, chan=chan, ctrl_only=TRUE)
      pat_data = getData_mats(fulldat, chan=chan)$pts
        
      priorpost(ctrl_data=ctrl_data, pat_data=pat_data,
                classif=output_reader(folder, fulldat, chan, out_type="CLASS")[[1]],
                prior=prior, post=output_reader(folder, fulldat, chan, out_type="POST"),
                chan=chan, mitochan="VDAC1", title=paste(gsub(".RAW.txt", "", fulldat), "\n", chan),
                reverse_transform=back_transform, fulldat=fulldat, folder=folder)
    }
  }
}
dev.off()

pdf(file.path("./PDF", folder, "pi_post.pdf"), width=13, height=8)
{   
  par(mfrow=c(1,1))
  for(fulldat in fulldats_raw){
    chans = getData_chanpats(fulldat)$channels
      for(chan in chans){
        pipost_plotter(folder, fulldat, chan, allData=TRUE) 
      }
    }
}
dev.off()

pdf(file.path("./PDF", folder, "classifs.pdf"), width=13, height=8)
{
  for(fulldat in fulldats_raw){
    chanpts = getData_chanpats(fulldat, get_Npats=TRUE)
    for(chan in chanpts$channels){
      ctrl_data =  getData_mats(fulldat, chan=chan, ctrl_only=TRUE)
      
      for(pat in chanpts$patients){
        pat_data = getData_mats(fulldat, chan=chan, pts=pat)$pts
        classifs = output_reader(folder, fulldat, chan, out_type="CLASS")[[1]]
        par(cex.lab=2, cex.main=2, cex.axis=1.5, mar=c(6,6,6,3))
        classif_plot(ctrl_data=ctrl_data, pat_data=pat_data,
                     classifs_pat=classifs[index_creator(chanpts$Npats, pat)],
                     chan=chan, mitochan="VDAC1", 
                     title=paste(gsub(".RAW.txt", "", fulldat), "\n", chan, pat))
        
      }
    }
  }
}
dev.off()












