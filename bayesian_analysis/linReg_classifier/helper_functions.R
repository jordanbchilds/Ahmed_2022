library(data.table)
library(MASS)

##
## COLOURS
##

pal = palette()
palette(c(pal, "darkblue", "darkorange", "deeppink"))

myBlack = function(alpha) rgb(0,0,0, alpha)
myDarkGrey = function(alpha) rgb(169/255,169/255,159/255, alpha)
myGrey = function(alpha) rgb(66/255,66/255,66/255, alpha)
myBlue = function(alpha) rgb(0,0,128/255, alpha)
myRed = function(alpha) rgb(1,0,0, alpha)

myGreen = function(alpha) rgb(0,100/255,0, alpha)
myYellow = function(alpha) rgb(225/255,200/255,50/255, alpha)
myPink = function(alpha) rgb(255/255,62/255,150/255, alpha)
myPurple = function(alpha) rgb(160/255, 32/255, 240/255, alpha)

cramp = colorRamp(c(myRed(0.2),myBlue(0.2)), alpha=TRUE)
# rgb(...) specifies a colour using standard RGB, where 1 is the maxColorValue
# 0.25 determines how transparent the colour is, 1 being opaque 
# cramp is a function which generates colours on a scale between two specifies colours

classcols = function(classif){
  # A function using the cramp specified colours to generate rgb colour names
  # input: a number in [0,1]
  # output: rgb colour name
  rgbvals = cramp(classif)/255.0
  return(rgb(rgbvals[,1],rgbvals[,2],rgbvals[,3], alpha=rgbvals[,4]))
}

colQuantiles = function(x, probs=0.5){
  quants = matrix(NA, nrow=ncol(x), ncol=length(probs))
  for(i in 1:ncol(x)){
    quants[i,] = quantile(x[,i], probs)
  }
  colnames(quants) = probs
  return(quants)
}

##
## GEN FUNCTIONS
##

vec_rep = function(x, n, byRow=TRUE){
  out = matrix(x, nrow=n, ncol=length(x), byrow=TRUE)
  if(byRow) return( out )
  else return( t(out) )
}

list2matrix = function(X, rowBind=TRUE){
  # X must list of matrices or a matrix
  if(is.matrix(X)) return(X)
  if(!is.list(X)) stop("X must be a list of matrices or a matrix")
  if(length(X)==1) return(X[[1]])
  
  n = length(X)
  out = X[[1]]
  if(rowBind){
    for(i in 2:n){ out = rbind(out, X[[i]]) }
    return(out)
  } else {
    for(i in 2:n){ out = cbind(out, X[[i]]) }
    return(out)
  }
}

log_transform = function(ctrl_mat, pat_mat=NULL){
  log_ctrl = log(ctrl_mat)
  if(!is.null(pat_mat)){
    return( list(ctrl=log_ctrl, pts=log(pat_mat)) )
  } else {
    return( log_ctrl )
  }
}

##
## DATA COLLEECTION
##

getData_mats = function(chan, mitochan="raw_porin", 
                        pts=NULL, ctrl_only=FALSE, 
                        data_transform=NULL,
                        get_patindex=FALSE, one_matrix=FALSE){
  
  data = read.csv(file.path("..", "rawdat_a.csv"), 
                  stringsAsFactors=FALSE, header=TRUE)

  Xctrl = data[grepl("C", data$caseno), mitochan]
  Yctrl = data[grepl("C", data$caseno), chan]
  ctrl_mat = cbind( Xctrl, Yctrl )
  
  if(!ctrl_only){
    pat_index = vector("numeric")
    pat_id = vector("character")
    ind = 1L
    if(is.null(pts)){
      Ypts = list()
      sbj = sort(unique(data$caseno))
      crl = grep("C", sbj, value = TRUE)
      pts = grep("P0.", sbj, value=TRUE)
      
      for(pat in pts){
        Xpat = data[data$caseno == pat, mitochan]
        Ypat = data[data$caseno == pat, chan]
        XY_pat = cbind(Xpat, Ypat)
        
        Ypts[[pat]] = XY_pat
        pat_index = c(pat_index, rep(ind, nrow(XY_pat)))
        pat_id = c(pat_id, rep(pat, nrow(XY_pat)))
        ind = ind + 1L
      }
    } else {
      Ypts = list()
      for(pat in pts){
        Xpat = data[data$caseno == pat, mitochan]
        Ypat = data[data$caseno == pat, chan]
        XY_pat = cbind(Xpat, Ypat)
        Ypts[[pat]] = XY_pat
        pat_index = c(pat_index, rep(ind, nrow(XY_pat)))
        pat_id = c(pat_id, rep(pat, nrow(XY_pat)))
        ind = ind + 1L
      }
    }
    pat_mat = list2matrix(Ypts)
  } else { pat_mat=NULL }
  
  if(one_matrix){
    mat = rbind(ctrl_mat, pat_mat)
    if( is.null(data_transform) ) return( mat )
    else return( data_transform(mat) )
  }
  if(get_patindex){
    if(!is.null(data_transform)) return( c(data_transform(ctrl_mat, pat_mat), list(pat_index=pat_index, pat_id=pat_id) ) )   
    if(ctrl_only) return( list(ctrl=ctrl_mat) )
    return( list(ctrl=ctrl_mat, pts=pat_mat, pat_index=pat_index, pat_id=pat_id) )
  } else {
    if(!is.null(data_transform)) return( data_transform(ctrl_mat, pat_mat) )
    if(ctrl_only) return( ctrl_mat )
    return( list(ctrl=ctrl_mat, pts=pat_mat) )
  }
}

##
## READ & SAVERS AOUTPUT
##

output_saver = function(outroot, output, folder, pat_only=FALSE){
  split = strsplit(outroot, split="_")[[1]]
  froot = split[1]
  chan = split[2]
  pat = split[3]
  
  if(pat_only){
    for(out_type in names(output)){
      filename = paste(froot, chan, pat, toupper(out_type), sep="_")
      write.table(output[[out_type]], paste0(file.path("Output", folder, filename), ".txt"),                    
                  row.names=FALSE, col.names=TRUE)
    }
    
  } else {
    for(ctrl_pat in c("ctrl", "pat")){
      out_ctrlpat = ifelse(ctrl_pat=="ctrl", "CONTROL", pat)
      for(out_type in names(output[[ctrl_pat]])){
        filename = paste(froot, chan, out_ctrlpat, toupper(out_type), sep="_")
        write.table(output[[ctrl_pat]][[out_type]], paste0(file.path("Output", folder, filename), ".txt"),
                    row.names=FALSE, col.names=TRUE)
      }
    }
  }
  
}

output_reader = function(folder, chan, pat="CONTROL", out_type){
  outroot = paste(chan, pat, sep="_")
  fp = file.path("Output", folder, paste0(outroot, "_", out_type, ".txt"))
  if(file.exists(fp)) return( read.table(fp, header=TRUE, stringsAsFactors=FALSE) )
  else stop(paste(fp, "file does not exist"))
}

##
## PLOTTERS
##

  
classif_plot = function(ctrl_data, pat_data, classifs_pat, chan, mitochan="VDAC1", title){
    op = par(mar=c(6,6,6,3), cex.axis=1.5, cex.lab=2, cex.main=2)
    xrange = range(c(ctrl_data[,1], pat_data[,1]))
    yrange = range(c(ctrl_data[,2], pat_data[,2]))
    
    plot(ctrl_data, pch=20, col=myDarkGrey(0.2), xlim=xrange, ylim=yrange, 
         xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"))
    points(pat_data, pch=20, col=classcols(classifs_pat))
    title(main=title, line=-4, outer=TRUE)
    par(op)
}
  
MCMCplot = function(folder, chan, pat="CONTROL", title="", lag=20){
    post = output_reader(folder, chan, pat, out_type="POST")
    
    if(pat=="CONTROL"){
      prior = read.table(file.path("Output", folder, paste(chan, "CONTROL", "PRIOR.txt", sep="_")), header=TRUE, stringsAsFactors = FALSE)
    } else {
      prior = output_reader(folder, chan, pat, out_type="PRIOR")
    }
    
    col.names = colnames(post)
    n.chains = length(post)
    par(mfrow=c(2,3), cex.main=2, cex.lab=2, cex.axis=1.5, mar=c(6,6,6,3))
    for(param in col.names){
      post_vec = post[[param]]
      plot(ts(post_vec), xlab="Iteration", ylab=paste(param), 
           main="", cex.lab=1.2)
      if(sum(post_vec==post_vec[1])!=length(post_vec)){
        acf(post[[param]], xlab="lag index", ylab="ACF", main="",
            cex.lab=1.2, lag=lag)
      } else {
        plot(NA, type='n', xlim=c(0,lag), ylim=c(0,1), 
             xlab="lag index", ylab="ACF", main="")
      }
      plot(density(post[[param]]), lwd=2, col="blue", xlab=paste(param), ylab="Density",
           main="")
      if(param %in% colnames(prior)) lines(density(prior[[param]]), lwd=2, col="green")
      
      title(main=title, line=-4, outer=TRUE)
    }
}

priorpost = function(ctrl_data, pat_data=NULL, priorpred, postpred,
                     classif=NULL, 
                     chan, mitochan="VDAC1", title="", xlims=NULL, ylims=NULL){
  
    op = par(mfrow=c(1,2), mar=c(6,6,6,3), cex.main=2, cex.lab=2, cex.axis=1.5)
    plot(ctrl_data, pch=20, cex=0.7, col=myGrey(0.1),
         xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), 
         main="Prior Predictive", xlim=xlims, ylim=ylims)
    if(!is.null(pat_data)) points(pat_data, pch=20, cex=1.2, col=myYellow(0.2))
    lines(priorpred[,"mitochan"], priorpred[,"lwr"], lty=2, col=myGreen(0.6), lwd=3)
    lines(priorpred[,"mitochan"], priorpred[,"mid"], lty=1, col=myGreen(0.6), lwd=4)
    lines(priorpred[,"mitochan"], priorpred[,"upr"], lty=2, col=myGreen(0.6), lwd=3)
    
    plot(ctrl_data, pch=20, col=myGrey(0.1),
         xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), 
         main="Posterior Predictive", xlim=xlims, ylim=ylims)
    if(!is.null(pat_data)) points(pat_data, pch=20, cex=1.2, col=classcols(classif))
    lines(priorpred[,"mitochan"], postpred[,"lwr"], lty=2, col=myPink(0.6), lwd=3)
    lines(priorpred[,"mitochan"], postpred[,"mid"], lty=1, col=myPink(0.6), lwd=4)
    lines(priorpred[,"mitochan"], postpred[,"upr"], lty=2, col=myPink(0.6), lwd=3)
    
    title(main=title, line=-2, outer=TRUE)
    
    par(op)
}

colvector_gen = function(pts){
  colind = double(length(pts))
  pts_B = unique(gsub("_S.", "", pts))
  for(i in 1:length(pts_B)){
    colind[ gsub("_S.", "", pts) == pts_B[i] ] = i + 1
  }
  colind
}

pipost_plotter = function(chan, folder, alpha=0.05){
  
  dat = read.csv(file.path("..", "rawdat_a.csv"), stringsAsFactors=FALSE)
  sbj = unique(dat$caseno)
  pts = sort(sbj[grepl("P0.", sbj)])
  
  npat = length(pts)
  pis = list()
  
  freq_list = list()
  for(pat in pts){
    pis[[pat]] = 1 - output_reader(folder, chan, pat, out_type="POST")[,"probdiff"]
  }
  title = paste0( chan )
  stripchart(pis, pch=20, method="jitter", vertical=TRUE, 
             col=rgb(t(col2rgb(palette()[colvector_gen(pts)]))/255, alpha=alpha), 
             group.names=pts,
             at = 1:npat,
             main=title, ylim=c(0.0,1.0), ylab="Proportion Deficiency", 
             xlab="Patient Sample")
}














































































