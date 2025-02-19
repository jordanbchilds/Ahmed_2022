#library("devtools")
#devtools::install_github("jordanbchilds/analysis2Dmito", force=TRUE)
library("analysis2Dmito")

library("parallel")
library("data.table")
library("rstan")
source("stan_sampler_function.R", local=TRUE)

folder = "stan_sampler_multiChain"

dir.create("Output", showWarnings = FALSE)
dir.create(file.path("Output", folder), showWarnings = FALSE)

data_cnr = as.data.frame(fread("DataChildsFormat_adj.txt"))

mitochan = "Porin"
channels = c("CI", "CIV")
nChan = length(channels)

dat_raw = data_cnr[ data_cnr$channel %in% c(mitochan, channels), ]

quant = 0.025
deltas = c(quantile(dat_raw[dat_raw$channel=="Porin", "value"],quant,na.rm=TRUE),
           quantile(dat_raw[dat_raw$channel=="CI", "value"],quant,na.rm=TRUE),
           quantile(dat_raw[dat_raw$channel=="CIV", "value"],quant,na.rm=TRUE))
names(deltas) = c("raw_porin","raw_CI","raw_CIV")

dat_adj = dat_raw

lower = -9999999
dat_adj[dat_adj$channel=="Porin", "value"] = pmax(lower,dat_raw[dat_raw$channel=="Porin","value"]-deltas["raw_porin"])
dat_adj[dat_adj$channel=="CI", "value"] = pmax(lower,dat_raw[dat_raw$channel=="CI", "value"]-deltas["raw_CI"])
dat_adj[dat_adj$channel=="CIV", "value"] = pmax(lower,dat_raw[dat_raw$channel=="CIV", "value"]-deltas["raw_CIV"])

#dat_adj[dat$channel=="CI", "value"] = dat$raw_CI-deltas["raw_CI"]
#dat_adj[dat$channel=="CIV", "value"] = dat$raw_CIV-deltas["raw_CIV"]

sbjs = unique(dat_adj$sampleID)
ctrlID = grep("C0", sbjs, value = TRUE)
pts = sort( grep("C0", sbjs, value=TRUE, invert=TRUE) )

grad = matrix(NA, nrow=length(channels), ncol=length(ctrlID))
colnames(grad) = ctrlID
rownames(grad) = channels

inter = grad
prec = grad

op = par(mfrow=c(1,1))
for( chan in channels ){
  for( crl in ctrlID ){
    xCtrl = dat_adj[dat_adj$channel==mitochan & dat_adj$sampleID==crl, "value"]
    yCtrl = dat_adj[dat_adj$channel==chan & dat_adj$sampleID==crl, "value"]
    dd = data.frame(mitochan=xCtrl, chan=yCtrl)
    xSyn = data.frame(mitochan=seq(min(dat_adj$value)-1, max(dat_adj$value)+1, length.out=1e3))
    
    mod = lm(chan ~ mitochan, data=dd)
    pred = predict.lm(mod, newdata=xSyn, interval="prediction")
    for( pat in pts ){
      xPat = dat_adj[dat_adj$channel==mitochan & dat_adj$sampleID==pat, "value"]
      yPat = dat_adj[dat_adj$channel==chan & dat_adj$sampleID==pat, "value"]
      
      plot(xCtrl,yCtrl, pch=20, col=alphaBlack(0.1),
           xlab=mitochan, ylab=chan,
           xlim=range(dat_adj$value), ylim=range(dat_adj$value))
      points(xPat, yPat, pch=20, col=alphaBlue(0.2))
      lines(xSyn$mitochan, pred[,"fit"], lty="solid", col="pink", lwd=3)
      lines(xSyn$mitochan, pred[,"lwr"], lty="dashed", col="pink", lwd=3)
      lines(xSyn$mitochan, pred[,"upr"], lty="dashed", col="pink", lwd=3)
      
      grad[chan, crl] = mod$coefficients[2]
      inter[chan, crl] = mod$coefficients[1]
      prec[chan, crl] = 1 / summary(mod)$sigma^2
    }
  }
}
par(op)

grad_mean = apply(grad, 1, mean)
inter_mean = apply(inter, 1, mean)
prec_mean = apply(prec, 1, mean)

grad_var = apply(grad, 1, var)
inter_var = apply(inter, 1, var)
prec_var = apply(prec, 1, var)

tau_c_vars = rep(50, nChan)^2 # c(2.5, 2, 2.5)
names(tau_c_vars) = channels

tau_m_vars = rep(50, nChan)^2 # c(10, 0.5, 1)
names(tau_m_vars) = channels

tau_vars = rep(50, nChan)^2 # c(3, 1.5, 3)
names(tau_vars) = channels

grad_var =  1 / rep(50, nChan)
inter_var = 1 / rep(50, nChan)
names(grad_var) = channels
names(inter_var) = channels

for (chan in channels) {
  mean_mu_m = grad_mean[chan]
  prec_mu_m = 1/0.25^2
  
  mean_mu_c = inter_mean[chan]
  prec_mu_c = 1/0.25^2
  
  tau_m_mode = 1 / grad_var[chan]
  tau_m_var = tau_m_vars[chan] # 1
  rate_tau_m = 0.5*(tau_m_mode + sqrt(tau_m_mode^2+4*tau_m_var)) / tau_m_var
  shape_tau_m = 1 + tau_m_mode*rate_tau_m
  
  tau_c_mode = 1 / inter_var[chan]
  tau_c_var = tau_c_vars[chan] # 1
  rate_tau_c = 0.5*(tau_c_mode + sqrt(tau_c_mode^2+4*tau_c_var)) / tau_c_var
  shape_tau_c = 1 + tau_c_mode*rate_tau_c
  
  tau_mode = prec_mean[chan]
  tau_var = tau_vars[chan] 
  rate_tau = 0.5*(tau_mode + sqrt(tau_mode^2+4*tau_var)) / tau_var
  shape_tau = 1 + tau_mode*rate_tau
  
  tau_def = 0.01
  
  paramVals = list(
    mean_mu_m=mean_mu_m, 
    prec_mu_m=prec_mu_m, 
    mean_mu_c=mean_mu_c, 
    prec_mu_c=prec_mu_c,
    shape_tau_m=shape_tau_m, 
    rate_tau_m=rate_tau_m, 
    shape_tau_c=shape_tau_c, 
    rate_tau_c=rate_tau_c, 
    shape_tau = shape_tau,
    rate_tau = rate_tau, 
    alpha_pi = 0.0,
    beta_pi = 0.5,
    tau_def=tau_def
  )
  
  data_list = list()
  for (i in seq_along(pts)) {
    data_list[[paste(chan, pts[i], sep = "__")]] = getData_mats(
      dat_adj,
      pts = pts[i],
      channels = c(mitochan, chan),
      ctrlID = ctrlID,
      getIndex = TRUE
    )
  }
  
  ncores = detectCores() - 1
  cl  = makeCluster(ncores)
  {
    clusterExport(cl, c("stan"))
    output = parLapply(
      cl,
      data_list,
      stan_inference,
      MCMCburnin = 1000,
      MCMCout = 1000,
      MCMCthin = 1,
      nChains = 10,
      max_logLik = FALSE,
      parameterVals = paramVals
    )
  }
  stopCluster(cl)
  
  for( rt in names(output) ){
    list_saver(output[[rt]], file.path("Output", folder, rt))
  }
}





