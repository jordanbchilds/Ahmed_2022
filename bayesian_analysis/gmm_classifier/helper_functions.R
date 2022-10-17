library(data.table)
library(MASS)

##
## GET DAT FUNCTIONS
##

getData_mats = function(chan, mitochan="raw_porin", 
                        pts=NULL, ctrl_only=FALSE, 
                        data_transform=NULL,
                        get_patindex=FALSE, one_matrix=FALSE){

  data = read.csv("../rawdat_a.csv", stringsAsFactors=FALSE)

  sbj = sort(unique(data$caseno))
  crl = grep("C", sbj, value = TRUE)
  pts = grep("P", sbj, value=TRUE)

  Xctrl = data[grepl("C", data$caseno), mitochan]
  Yctrl = data[grepl("C", data$caseno), chan]
  ctrl_mat = cbind( Xctrl, Yctrl )
  
  if(!ctrl_only){
    pat_index = vector("numeric")
    pat_id = vector("character")
    ind = 1L
    if(is.null(pts)){
      Ypts = list()
      for(pat in pts){
        Xpat = data[data$caseno==pat, mitochan]
        Ypat = data[data$caseno==pat, chan]
        XY_pat = cbind(Xpat, Ypat)
        Ypts[[pat]] = XY_pat
        pat_index = c(pat_index, rep(ind, nrow(XY_pat)))
        pat_id = c(pat_id, rep(pat, nrow(XY_pat)))
        ind = ind + 1L
      }
    } else {
      Ypts = list()
      for(pat in pts){
        Xpat = data[data$caseno==pat, mitochan]
        Ypat = data[data$caseno==pat, chan]
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
## COLOUR FUNCTIONS

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

##
## GENERAL FUNCTIONS FOR bayes GMM 
##

colQuantiles = function(x, probs=0.5){
  quants = matrix(NA, nrow=ncol(x), ncol=length(probs))
  for(i in 1:ncol(x)){
    quants[i,] = quantile(x[,i], probs)
  }
  colnames(quants) = probs
  return(quants)
}

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

##
## DATA TRANSFORMATION
##

rotate_mat = function(X, R, reverse=FALSE){
  if(!reverse) return( X%*%R )
  X%*%solve(R)
}

centre_mat = function(X, centre=NULL, reverse=FALSE){
  if(!reverse){
    if(!is.null(centre)){
      return( X - vec_rep(centre, nrow(X)))
    } 
    return( X - vec_rep(colMeans(X), nrow(X)) )
  }
  if(is.null(centre)){ stop("Reverse calculations require centre") }
  X + vec_rep(centre, nrow(X))
}

scale_mat = function(X, scale=NULL, reverse=FALSE){
  if(!reverse){
    if(!is.null(scale)){
      return( X / vec_rep(scale, nrow(X)) )
    }
    return( X / vec_rep(sqrt(diag(var(X))), nrow(X)) )
  }
  if(is.null(scale)){ stop("Reverse calculations require scale") }
  X * vec_rep(scale, nrow(X))
}

myData_transform = function(ctrl_mat, pat_mat=NULL){
  pca = prcomp(ctrl_mat, scale=FALSE, center=FALSE)
  ctrl_mat = pca$x
  # ctrl_mat = scale(ctrl_mat, center=TRUE, scale=TRUE)
  ctrl_mean = colMeans(ctrl_mat)
  ctrl_mat = centre_mat(ctrl_mat, ctrl_mean)
  ctrl_sd = sqrt(diag(var(ctrl_mat)))
  ctrl_mat = scale_mat(ctrl_mat, ctrl_sd)
  
  if(!is.null(pat_mat)){
    pat_mat = rotate_mat(pat_mat, pca$rotation)
    pat_mat = centre_mat(pat_mat, ctrl_mean)
    pat_mat = scale_mat(pat_mat, ctrl_sd)
    return(list(ctrl=ctrl_mat, pts=pat_mat))
  }
  ctrl_mat
}

back_transform = function(X, ctrl_mat=NULL, fulldat=NULL, folder=NULL, chan=NULL, pat=NULL, 
                          parameters=NULL){
  if( !is.null(fulldat) | !is.null(folder) | !is.null(chan) | !is.null(pat)){ 
    ctrl_mat = getData_mats(fulldat, chan=chan, pts=pat, ctrl_only=TRUE, data_transform=log_transform)
  }
  if(!is.null(ctrl_mat)){
    pca = prcomp(ctrl_mat, scale=FALSE, center=FALSE)
    ctrl_mat = pca$x
    ctrl_mean = colMeans(ctrl_mat)
    ctrl_sd = sqrt(diag(var(ctrl_mat)))
    
    if(!is.null(parameters)){ list(seq(1,ncol(X))) }
    Xnew = X
    for(i in 1:length(parameters)){
      params = parameters[[i]]
      Xnew[, params] = scale_mat(X[,params], scale=ctrl_sd, reverse=TRUE)
      Xnew[, params] = centre_mat(Xnew[, params], centre=ctrl_mean, reverse=TRUE)  
      Xnew[, params] = rotate_mat(as.matrix(Xnew[,params]), R=pca$rotation, reverse=TRUE)
    }
  } else {
    stop("Must supply ctrl_raw or fulldat, folder, chan and pat")
  }
  Xnew
}

log_transform = function(ctrl_mat, pat_mat=NULL){
  if(!is.null(pat_mat)) return( list(ctrl=log(ctrl_mat), pts=log(pat_mat)) )
  log(ctrl_mat)
}

##
## SAVE & READ OUTPUT
##

output_saver = function(outroot, output, folder){
  split = strsplit(outroot, split="_")[[1]]
  froot = split[1]
  chan = split[2]
  pat = split[3]
  
  for(ctrl_pat in c("ctrl", "pat")){
    out_ctrlpat = ifelse(ctrl_pat=="ctrl", "CONTROL", pat)
    for(out_type in names(output[[ctrl_pat]])){
      filename = paste(froot, chan, out_ctrlpat, toupper(out_type), sep="_")
      write.table(output[[ctrl_pat]][[out_type]], paste0(file.path("Output", folder, filename), ".txt"),
                  row.names=FALSE, col.names=TRUE)
    }
  }
}

output_reader = function(folder, fulldat, chan, pat=NULL, out_type){
  if(is.null(pat)){
    outroot = paste(gsub(".RAW.txt", "", fulldat), chan, sep="_")
  } else {
    outroot = paste(gsub(".RAW.txt", "", fulldat), chan, gsub("_",".", pat), sep="_")
  }
  fp = file.path("./Output", folder, paste0(outroot, "_", out_type, ".txt"))
  if(file.exists(fp)) return( read.table(fp, header=TRUE, stringsAsFactors=FALSE) )
  else stop("file does not exist")
}

##
## PLOTTING FUNCTIONS
##

colvector_gen = function(pts){
  colind = double(length(pts))
  pts_B = unique(gsub("_S.", "", pts))
  for(i in 1:length(pts_B)){
    colind[ gsub("_S.", "", pts) == pts_B[i] ] = i + 1
  }
  colind
}

pipost_plotter = function(folder, fulldat, chan, alpha=0.05, allData=TRUE){
  pts = getData_chanpats(fulldat)$patients
  pts_blocks = unique(gsub("_S.", "", pts))
  npat = length(pts)
  
  pis = list()
  if(allData){
    post = output_reader(folder, fulldat, chan, out_type="POST")
    for(i in 1:npat){
      pis[[pts[i]]] = post[,paste0("probctrl.",i,".")]
    }
  } else {
    for(pat in pts){
      pis[[pat]] = output_reader(folder, fulldat, chan, out_type="POST")[,"probctrl"]
    }
  }

  pat_labels = as.vector(rbind("", unique(gsub("_S.", "", gsub("P._", "", pts))), ""))
  
  op = par(mar=c(6,6,6,3), cex.main=2, cex.lab=2, cex.axis=1.5)
  
  stripchart(pis, pch=20, method="jitter", vertical=TRUE, 
             col=rgb(t(col2rgb(palette()[colvector_gen(pts)]))/255, alpha=alpha), 
             group.names=rep("", length(pts)),
             main="", ylim=c(0,1), ylab="like control proportion", 
             xlab="Patient Sample")
  if(!all(substr(pts, 4,4)==substr(pts,4,4)[1])){
    sep = sum(substr(pts, 4,4)==substr(pts,4,4)[1])+0.5
    arrows(x0=sep, x1=sep, y0=0, y1=1, lwd=3, lty=2, length=0, 
           col=myDarkGrey(1))
  }
  text(1:length(pts), y=0, labels=pat_labels, pos=3, cex=1)
  title(main=paste0(substr(pts[1],1,2), " ", chan, "\n", gsub(".RAW.txt", "", fulldat) ), 
        line=-4, outer=TRUE)
  par(op)
}

percentiles = function(xdat, ydat, probs=c(0.975, 0.5, 0.025)){
  dens = kde2d(xdat, ydat, n=200); ## estimate the z counts
  dx = diff(dens$x[1:2])
  dy = diff(dens$y[1:2])
  sz = sort(dens$z)
  c1 = cumsum(sz) * dx * dy
  levs = sapply(probs, function(x) {
    approx(c1, sz, xout = 1 - x)$y
  })
  return( list(dens=dens, levels=levs, probs=probs))
}

priorpost = function(ctrl_data, pat_data, prior=NULL, post, classifs_pat, title="",
                     mitochan="VDAC1", chan, pat=NULL, reverse_transform=NULL, ...){
  xlims = range(c(ctrl_data[,1], pat_data[,1]))
  ylims = range(c(ctrl_data[,2], pat_data[,2]))
  
  if(!is.null(reverse_transform)){
    prior = back_transform(prior, fulldat=fulldat, folder=folder, chan=chan, pat=pat, 
                           parameters=list(c("comp.1.1.", "comp.2.1."), c("comp.1.2.", "comp.2.2.")))
    post = back_transform(post, fulldat=fulldat, folder=folder, chan=chan, pat=pat, 
                          parameters=list(c("comp.1.1.", "comp.2.1."), c("comp.1.2.", "comp.2.2.")))
  }
  
  op = par(mfrow=c(1,2), mar=c(6,6,6,3), cex.main=2, cex.lab=2, cex.axis=1.5)
  if(!is.null(prior)){
    op = par(mfrow=c(2,2))
    plot(ctrl_data, pch=20, col=myDarkGrey(0.3), 
         xlab=paste0("log(", mitochan,")"), ylab=paste0("log(",chan,")"),
         xlim=xlims, ylim=ylims)
    points(pat_data, pch=20, col=myYellow(0.3))
    contours = percentiles(prior[,"comp.1.1."], prior[,"comp.2.1."])
    contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2 )
    
    plot(ctrl_data, pch=20, col=myDarkGrey(0.3), 
         xlab=paste0("log(", mitochan,")"), ylab=paste0("log(",chan,")"),
         xlim=xlims, ylim=ylims)
    points(pat_data, pch=20, col=myYellow(0.3))
    contours = percentiles(prior[["comp.1.2."]], prior[["comp.2.2."]])
    contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2 )
  }
  plot(ctrl_data, pch=20, col=myDarkGrey(0.3), 
       xlab=paste0("log(", mitochan,")"), ylab=paste0("log(",chan,")"),
       xlim=xlims, ylim=ylims)
  points(pat_data, pch=20, col=classcols(classifs_pat) )
  contours = percentiles(post[,"comp.1.1."], post[,"comp.2.1."])
  contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2 )
  
  plot(ctrl_data, pch=20, col=myDarkGrey(0.3), 
       xlab=paste0("log(", mitochan,")"), ylab=paste0("log(",chan,")"),
       xlim=xlims, ylim=ylims)
  points(pat_data, pch=20, col=classcols(classifs_pat) )
  contours = percentiles(post[,"comp.1.2."], post[,"comp.2.2."])
  contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2 )
  
  title(main=title, line=-4, outer=TRUE)
  
  par(op)
}

MCMCplot = function(folder, fulldat, chan, pat=NULL, title="", lag=20){
  post = output_reader(folder, fulldat, chan, pat, out_type="POST")
  prior = read.table(file.path("./Output", folder, "PRIOR.txt"), 
                     header=TRUE, stringsAsFactors=FALSE)

  col.names = colnames(post)
  n.chains = length(post)
  op = par(mfrow=c(2,3), cex.main=2, cex.lab=2, cex.axis=1.5,
           mar=c(6,6,6,3))
  for(param in col.names){
    post_vec = post[, param]
    plot(ts(post_vec), xlab="Iteration", ylab=paste(param), 
         main="", cex.lab=1.2)
    if(sum(post_vec==post_vec[1])!=length(post_vec)){
      acf(post_vec, xlab="lag index", ylab="ACF", main="",
          cex.lab=1.2, lag=lag)
    } else {
      plot(NA, type='n', xlim=c(0,lag), ylim=c(0,1), 
           xlab="lag index", ylab="ACF", main="")
    }
    plot(density(post_vec), lwd=2, col="blue", xlab=paste(param), ylab="Density",
         main="")
    if(param %in% colnames(prior)) lines(density(prior[,param]), lwd=2, col="green")
    
    title(main=title, line=-4, outer=TRUE)
  }
  par(op)
}

classif_plot = function(ctrl_data=NULL, pat_data, classifs_pat, 
                        title="", mitochan="VDAC1", chan, pat){
  xlims = range(c(ctrl_data[,1], pat_data[,1]))
  ylims = range(c(ctrl_data[,2], pat_data[,2]))
  
  op = par(mfrow=c(1,1), mar=c(6,6,6,3), cex.main=2, cex.lab=2, cex.axis=1.5)
  
  plot(ctrl_data, pch=20, col=myDarkGrey(0.2), 
       xlab=paste0("log(", mitochan,")"), ylab=paste0("log(",chan,")"),
       xlim=xlims, ylim=ylims)
  points(pat_data, pch=20, col=classcols(classifs_pat) )
  
  title(main=title, line=-4, outer=TRUE)
  par(op)
}

##
## all data specialty functions
##

index_creator = function(Npops, ind){
  index = rep(FALSE, sum(Npops))
  if( is.character(ind) ){
    ind = which(names(Npops)==pat)
    names(Npops) = NULL
  } else if( is.null(names(Npops)) & is.character(ind)){
    stop("if ind is a character Npops must be a named vector")
  }
  if(ind>1 & ind<length(Npops)){
    return((sum(Npops[1:(ind-1)])+1):sum(Npops[1:ind]))
  } else if(ind==1){
    return(1:Npops[1])
  } else {
    return( (sum(Npops[1:(ind-1)])+1):sum(Npops) ) 
  }
}

classif_finder = function(classifs_all, fulldat, pat){
  Npats = getData_chanpats(fulldat, get_Npats=TRUE)$Npats
  
  classifs_all[index_creator(Npats, pat)]
}









































































