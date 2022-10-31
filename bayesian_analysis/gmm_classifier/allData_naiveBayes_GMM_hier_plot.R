library("MASS")
source("./helper_functions.R", local = TRUE)

folder = "allData_naiveBayes_GMM_hier"

dir.create(file.path("PDF"), showWarnings = FALSE)
dir.create(file.path("PDF", folder), showWarnings = FALSE)

cord = c("raw_CI", "raw_CIV")
mitochan = "porin"

dat = read.csv(file.path("..", "rawdat_a.csv"), stringsAsFactors=FALSE)
sbj = unique(dat$caseno)
crl = sbj[grepl("C0.", sbj)]
pts = sort(sbj[grepl("P0.", sbj)])

prior = read.table(file.path("Output", folder, "PRIOR.txt"), 
                   header=TRUE, stringsAsFactors=FALSE)

pdf(file.path("./PDF", folder, "MCMC.pdf"), width=13, height=8)
{ 
  for(chan in cord){
    MCMCplot(folder, chan,  lag=100,
             title=chan)
  }
}
dev.off()

pdf(file.path("./PDF", folder, "model_post_RAW.pdf"), width=13, height=8)
{
  for(chan in cord){
    ctrl_data =  getData_mats(chan=chan, ctrl_only=TRUE, 
                              data_transform=myData_transform)

    pat_data = getData_mats(chan=chan,
                            data_transform=myData_transform)$pts
      
    classif_allpats = output_reader(folder, chan, out_type="CLASS")[[1]]
    post = output_reader(folder, chan, out_type="POST")
        
    priorpost(ctrl_data=ctrl_data, pat_data=pat_data,
              classif=classif_allpats,
              prior=prior, post=post,
              chan=chan, mitochan="porin", 
              title=chan)
    }
}
dev.off()

pdf(file.path("PDF", folder, "model_post.pdf"), width=13, height=8)
{
  for(chan in cord){
    ctrl_data = getData_mats(chan=chan, ctrl_only=TRUE)
    pat_data = getData_mats(chan=chan)$pts
    
    classif = output_reader(folder, chan, out_type="CLASS")[[1]]
    post = output_reader(folder, chan, out_type="POST")
        
    priorpost(ctrl_data=ctrl_data, pat_data=pat_data,
              classif=classif,
              prior=prior, post=post,
              chan=chan, mitochan="porin", title=chan,
              reverse_transform=back_transform)
  }
}
dev.off()

pdf(file.path("./PDF", folder, "pi_post.pdf"), width=13, height=8)
{   
  op = par(mfrow=c(1,1))
  for(chan in cord){
    pipost_plotter(folder, chan, allData=TRUE) 
  }
  par(op)
}
dev.off()

pdf(file.path("./PDF", folder, "classifs.pdf"), width=13, height=8)
{
  op = par(cex.lab=2, cex.main=2, cex.axis=1.5, mar=c(6,6,6,3))
   for(chan in cord){
    ctrl_data =  getData_mats(chan=chan, ctrl_only=TRUE)
    for(pat in pts){
      pat_id = getData_mats(chan=chan, get_patindex=TRUE)$pat_id
      pat_data = getData_mats(chan=chan, pts=pat)$pts
      classifs = output_reader(folder, chan, out_type="CLASS")[[1]]
      
      classif_pat = classifs[pat_id==pat]
      classif_plot(ctrl_data=ctrl_data, pat_data=pat_data,
                   classifs_pat=classif_pat,
                   chan=chan, mitochan="porin", 
                   title=paste(chan, pat))
    }
   }
  par(op)
}
dev.off()












