library(rjags)
library(MASS)
library(parallel)
source("helper_functions.R", local=TRUE)

folder = "linReg_classifier"

dir.create(file.path("./Output"), showWarnings = FALSE)
dir.create(file.path("./Output", folder), showWarnings = FALSE)

modelstring = "
 model {
  # Likelihood of data given model parameters
  for (i in 1:N){
   Yobs[i] ~ dnorm(m*Xobs[i] + c, tau_hat[i])
   tau_hat[i] = ifelse(class[i]==1, tau_notctrl, tau_ctrl)
   class[i] ~ dbern(probdiff)
  }
  for(j in 1:Nsyn){
   Ysyn_norm[j] ~ dnorm(m*Xsyn[j]+c, tau_ctrl)
   Ysyn_def[j] ~ dnorm(m*Xsyn[j]+c, tau_notctrl)
  }
  # Specify prior beliefs about parameters
  m ~ dnorm(mu_m, tau_m)
  c ~ dnorm(mu_c, tau_c)
  tau_ctrl ~ dgamma(shape_tau, rate_tau)
  notctrl_var ~ dgamma(shape_gamma, rate_gamma)
  tau_notctrl = 1/notctrl_var
  p ~ dbeta(alpha, beta)
  probdiff = ifelse(ctrl_ind==1, 0, p)
 }
"

cord = c("raw_CI", "raw_CIV")
mitochan = "raw_porin"

dat = read.csv("../rawdat_a.csv", stringsAsFactors=FALSE)
sbj = unique(dat$caseno)
crl = sbj[grepl("C0.", sbj)]
pts = sort(sbj[grepl("P0.", sbj)])

# inferenecce 
inference = function(input){
  with(c(input),{
    data_mats = getData_mats(chan=chan, pts=pat)
    ctrl_mat = data_mats$ctrl
    pat_mat = data_mats$pts
    Nctrl = nrow(ctrl_mat)
    Npat = nrow(pat_mat)
    
    # prior parameters for control data
    c_est = 0
    tau_c = 1/2^2
    m_est = 0
    tau_m = 1/2^2
    shape_tau = 3
    rate_tau = 1
    
    gamma_mean = 20
    gamma_var = 10^2
    
    shape_gamma = gamma_mean^2 / gamma_var
    rate_gamma = gamma_mean / gamma_var

    N_syn = 1e4
    X_syn = seq(0, max(c(ctrl_mat[,1], pat_mat[,1]))*1.5, length.out=N_syn)
    
    ## control inference
    data_ctrl = list( Xobs=ctrl_mat[,1], Yobs=ctrl_mat[,2], N=Nctrl, Nsyn=N_syn,
                      Xsyn=X_syn,
                      mu_m=m_est, tau_m=tau_m, mu_c=c_est, tau_c=tau_c,
                      shape_tau=shape_tau, rate_tau=rate_tau,
                      shape_gamma=shape_gamma, rate_gamma=rate_gamma,
                      alpha=1.0, beta=1.0,
                      ctrl_ind=1)
    
    data_ctrl_priorpred = data_ctrl
    data_ctrl_priorpred$Yobs = NULL
    data_ctrl_priorpred$N = 0
    
    model_ctrl = jags.model(textConnection(modelstring), data=data_ctrl, n.chains=n.chains)
    model_ctrl_priorpred = jags.model(textConnection(modelstring), data=data_ctrl_priorpred)
    
    update(model_ctrl,n.iter=MCMCBurnin)
    output_ctrl = coda.samples(model=model_ctrl, n.iter=MCMCOut*MCMCThin,thin=MCMCThin,
                               variable.names=c("m","c","tau_ctrl","tau_notctrl","Ysyn_norm","Ysyn_def", "class","probdiff"))
    output_ctrl_prior = coda.samples(model=model_ctrl_priorpred, n.iter=MCMCOut,thin=1,
                                     variable.names=c("m","c","tau_ctrl","tau_notctrl","Ysyn_norm","Ysyn_def","probdiff"))
    
    posterior_ctrl = as.data.frame(output_ctrl[[1]])
    prior_ctrl = as.data.frame(output_ctrl_prior[[1]])
    summ_ctrl = summary(output_ctrl)
    classifs_ctrl = summ_ctrl$statistics[grepl("class",rownames(summ_ctrl$statistics)),"Mean"]
    
    ### patient inference
    # pateint priors
    flex = 0.01
    c_est = mean(posterior_ctrl$c)
    tau_c = flex / var(posterior_ctrl$c)
    m_est = mean(posterior_ctrl$m)
    tau_m = flex / var(posterior_ctrl$m)
    tau_mean = mean(posterior_ctrl$tau_ctrl) # Precision tau = (1/sd)^2
    tau_var = var(posterior_ctrl$tau_ctrl) # Deviation from prior tau should require a lot of contradictory data
    tau_shape =  tau_mean^2 / tau_var
    tau_rate = tau_mean / tau_var
    
    gamma_mean = 20 
    gamma_var = 10^2
    gamma_shape = gamma_mean^2 / gamma_var
    gamma_rate = gamma_mean / gamma_var
    
    data_pat = list(Xobs=pat_mat[,1], Yobs=pat_mat[,2], N=Npat, Nsyn=N_syn,
                    Xsyn=X_syn, 
                    mu_m=m_est, tau_m=tau_m, 
                    mu_c=c_est, tau_c=tau_c,
                    shape_tau=tau_shape, rate_tau=tau_rate, 
                    shape_gamma=gamma_shape, rate_gamma=gamma_rate, 
                    alpha=1, beta=1,
                    ctrl_ind=0)
    
    data_pat_priorpred = data_pat
    data_pat_priorpred$Yobs = NULL
    data_pat_priorpred$N = 0
    
    model_pat = jags.model(textConnection(modelstring), data=data_pat, n.chains=n.chains)
    model_pat_priorpred = jags.model(textConnection(modelstring), data=data_pat_priorpred)
    update(model_pat, n.iter=MCMCBurnin)

    output_pat_post = coda.samples(model=model_pat, n.iter=MCMCOut*MCMCThin, thin=MCMCThin,
                            variable.names=c("m","c","tau_ctrl","tau_notctrl","Ysyn_norm", "Ysyn_def","class","probdiff"))
    output_pat_prior = coda.samples(model=model_pat_priorpred,n.iter=MCMCOut, thin=1,
                                      variable.names=c("m","c","tau_ctrl","tau_notctrl","Ysyn_norm", "Ysyn_def","probdiff"))
    
    posterior_pat = as.data.frame(output_pat_post[[1]])
    prior_pat = as.data.frame(output_pat_prior[[1]])
    
    summ_pat = summary(output_pat_post)
    classifs_pat = summ_pat$statistics[grepl("class",rownames(summ_pat$statistics)),"Mean"]
    
    posterior_ctrl_names = colnames(posterior_ctrl)
    post_ctrl = posterior_ctrl[,c("m", "c", "tau_ctrl", "tau_notctrl", "probdiff")]
    postpred_ctrl_norm = colQuantiles( posterior_ctrl[,grepl("Ysyn_norm", posterior_ctrl_names)], probs=c(0.025, 0.5, 0.975) )
    postpred_ctrl_def = colQuantiles( posterior_ctrl[,grepl("Ysyn_def", posterior_ctrl_names)], probs=c(0.025, 0.5, 0.975) )
    postpred_ctrl = cbind(X_syn, postpred_ctrl_norm, postpred_ctrl_def)
    colnames(postpred_ctrl) = c("mitochan", "lwr_norm", "mid_norm", "upr_norm", "lwr_def", "mid_def", "upr_def")
    
    prior_ctrl_names = colnames(prior_ctrl)
    priorpred_ctrl_norm = colQuantiles(prior_ctrl[, grepl("Ysyn_norm", prior_ctrl_names)], probs=c(0.025,0.5,0.975))
    priorpred_ctrl_def = colQuantiles(prior_ctrl[, grepl("Ysyn_def", prior_ctrl_names)], probs=c(0.025,0.5,0.975))
    priorpred_ctrl = cbind(X_syn, priorpred_ctrl_norm, priorpred_ctrl_def)
    colnames(priorpred_ctrl) = c("mitochan", "lwr_norm", "mid_norm", "upr_norm", "lwr_def", "mid_def", "upr_def")
    prior_control = prior_ctrl[,c("m", "c", "tau_ctrl", "tau_notctrl", "probdiff")]
    
    posterior_pat_names = colnames(posterior_pat)
    post_pat = posterior_pat[,c("m", "c", "tau_ctrl", "tau_notctrl", "probdiff")]
    postpred_pat_norm = colQuantiles(posterior_pat[,grepl("Ysyn_norm", posterior_pat_names)], probs=c(0.025,0.5,0.975))
    postpred_pat_def = colQuantiles(posterior_pat[,grepl("Ysyn_def", posterior_pat_names)], probs=c(0.025,0.5,0.975))
    postpred_pat = cbind(X_syn, postpred_pat_norm, postpred_pat_def)
    colnames(postpred_pat) = c("mitochan", "lwr_norm", "mid_norm", "upr_norm", "lwr_def", "mid_def", "upr_def")
    
    prior_pat_names = colnames(prior_pat)
    priorpred_pat_norm = colQuantiles(prior_pat[,grepl("Ysyn_norm", prior_pat_names)], probs=c(0.025,0.5,0.975))
    priorpred_pat_def = colQuantiles(prior_pat[,grepl("Ysyn_def", prior_pat_names)], probs=c(0.025,0.5,0.975))
    priorpred_pat = cbind(X_syn, priorpred_pat_norm, priorpred_pat_def)
    colnames(priorpred_pat) = c("mitochan", "lwr_norm", "mid_norm", "upr_norm", "lwr_def", "mid_def", "upr_def")
    
    prior_patient = prior_pat[,c("m", "c", "tau_ctrl", "tau_notctrl", "probdiff")]
    
    ctrl_list = list(post=post_ctrl, postpred=postpred_ctrl, 
                     prior=prior_control, priorpred=priorpred_ctrl,
                     classif=classifs_ctrl)
    pat_list = list(post=post_pat, postpred=postpred_pat,
                    prior=prior_patient, priorpred=priorpred_pat,
                    classif=classifs_pat)
    
    return( list(ctrl=ctrl_list, pat=pat_list) )
  })
}

inputs = list()
{
  input0 = list()
  input0$MCMCOut = 2000
  input0$MCMCBurnin = 2000
  input0$MCMCThin = 100
  input0$n.chains = 1

    for(chan in cord){
      for(pat in pts){
        outroot = paste(chan, pat, sep="_")
        inputs[[outroot]] = input0
        inputs[[outroot]]$chan = chan
        inputs[[outroot]]$pat = pat
      } # pts
    } # chans
}

ncores = detectCores() - 1
cl  = makeCluster(ncores)
{
  clusterExport(cl, c("modelstring"))
  clusterEvalQ(cl, {
    library("rjags")
    library("data.table")
    source("helper_functions.R", local=TRUE)
  })
  linreg_output = parLapply(cl, inputs, inference)
}
stopCluster(cl)

for(outroot in names(linreg_output)){
  output_saver(outroot, linreg_output[[outroot]], folder)
}








