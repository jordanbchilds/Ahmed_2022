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
  
  Xctrl = data[grepl("C", data$caseno), mitochan]
  Yctrl = data[grepl("C", data$caseno), chan]
  ctrl_mat = cbind( Xctrl, Yctrl )
  
  if(!ctrl_only){
    pat_index = vector("numeric")
    pat_id = vector("character")
    ind = 1L
    if(is.null(pts)){
      
      sbj = sort(unique(data$caseno))
      pts = grep("P", sbj, value=TRUE)
      
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

back_transform = function(X, ctrl_mat){
  pca = prcomp(ctrl_mat, scale=FALSE, center=FALSE)
  ctrl_mat = pca$x
  ctrl_mean = colMeans(ctrl_mat)
  ctrl_sd = sqrt(diag(var(ctrl_mat)))
  
  Xnew = scale_mat(X, scale=ctrl_sd, reverse=TRUE)
  Xnew = centre_mat(Xnew, centre=ctrl_mean, reverse=TRUE)  
  Xnew[, grepl("comp\\.1", colnames(Xnew))] = rotate_mat(as.matrix(Xnew[,grepl("comp\\.1", colnames(Xnew))]), R=pca$rotation, reverse=TRUE)     
  Xnew[, grepl("comp\\.2", colnames(Xnew))] = rotate_mat(as.matrix(Xnew[,grepl("comp\\.2", colnames(Xnew))]), R=pca$rotation, reverse=TRUE)     
  
  return(Xnew)
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
  chan = split[2]
  pat = split[3]
  
  for(ctrl_pat in c("ctrl", "pat")){
    out_ctrlpat = ifelse(ctrl_pat=="ctrl", "CONTROL", pat)
    for(out_type in names(output[[ctrl_pat]])){
      filename = paste(chan, out_ctrlpat, toupper(out_type), sep="_")
      write.table(output[[ctrl_pat]][[out_type]], paste0(file.path("Output", folder, filename), ".txt"),
                  row.names=FALSE, col.names=TRUE)
    }
  }
}

output_reader = function(folder, chan, out_type){

  fp = file.path("./Output", folder, paste0(chan, "_", out_type, ".txt"))
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

pipost_plotter = function(folder, chan, alpha=0.02, allData=TRUE){
  
  dat = read.csv(file.path("..", "rawdat_a.csv"), stringsAsFactors=FALSE)
  sbj = unique(dat$caseno)
  pts = sort(sbj[grepl("P0.", sbj)])
  npat = length(pts)
  
  pis = list()
  if(allData){
    post = output_reader(folder, chan, out_type="POST")
    for(i in 1:npat){
      pis[[pts[i]]] = post[,paste0("probDef.",i,".")]
    }
  } else {
    for(pat in pts){
      pis[[pat]] = output_reader(folder, chan, out_type="POST")[,"probDef"]
    }
  }

  op = par(mar=c(6,6,6,3), cex.main=2, cex.lab=2, cex.axis=1.5)
  
  stripchart(pis, pch=20, method="jitter", vertical=TRUE, 
             col=rgb(t(col2rgb(palette()[colvector_gen(pts)]))/255, alpha=alpha), 
             group.names=pts,
             main="", ylim=c(0,1), ylab="Deficiency Proportion", 
             xlab="Patient Sample")
  title(main=chan, line=-2, outer=TRUE)
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
                     mitochan="porin", chan, pat=NULL, reverse_transform=NULL ){
  
  ctrl_raw = getData_mats(chan, ctrl_only=TRUE)
  
  xlims = range(c(ctrl_data[,1], pat_data[,1]))
  ylims = range(c(ctrl_data[,2], pat_data[,2]))
  
  if(!is.null(reverse_transform)){
    prior = reverse_transform(prior[,grepl("comp", colnames(prior))], ctrl_raw)
    post = reverse_transform(post[grepl("comp", colnames(post))], ctrl_raw)
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

MCMCplot = function(folder, chan, pat=NULL, title="", lag=20){
  
  post = output_reader(folder, chan, out_type="POST")
  prior = read.table(file.path("Output", folder, "PRIOR.txt"), 
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
                        title="", mitochan="raw_porin", chan, pat){
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










































































