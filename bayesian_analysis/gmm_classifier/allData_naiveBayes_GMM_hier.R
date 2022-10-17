library("parallel")
library("rstan")
source("helper_functions.R", local=TRUE)

folder = "allData_naiveBayes_GMM_hier"
dir.create("Output", showWarnings = FALSE)
dir.create(file.path("Output", folder), showWarnings = FALSE)

cord = c("raw_CI", "raw_CIV")
mitochan = "raw_porin"

# inferenecce 
inference = function(input){
  with(c(input),{
    data_mats = getData_mats(chan=chan, 
                             data_transform=myData_transform, get_patindex=TRUE)
    Yctrl = data_mats$ctrl
    Ypat = data_mats$pts
    Nctrl = nrow(Yctrl)
    Npat = nrow(Ypat)
    pat_index = data_mats$pat_index
    n_pts = length(unique(pat_index))
    data_list = list(Yctr=Yctrl, Ypat=Ypat, Nctrl=Nctrl, Npat=Npat, 
                pat_index=pat_index, n_pts=n_pts)
    
    output = stan("./allData_GMM_hier.stan", data=c(prior_list, data_list), chains=1, 
                  iter=(MCMCOut+MCMCBurnin)*MCMCThin, warmup=MCMCBurnin, thin=MCMCThin )
    
    outmat = as.matrix(output)
    outcols = colnames(outmat)
    
    classifs_df = outmat[, grepl("classif", outcols)]
    classifs_avg = colMeans(classifs_df)
    post = outmat[, !(grepl("classif", outcols)|grepl("lp__", outcols)|
                        grepl("probvec", outcols)|grepl("dens", outcols)|
                        grepl("z", outcols))]
    
    return( list(post=post, classifs=classifs_avg) )
  })
}

{
  p = 2
  # prior parameters
  mean_mu1 = double(p)
  mean_mu2 = double(p)
  var_mu1=  matrix(c(0.5^2, 0, 0, 0.5^2), ncol=p, nrow=p, byrow=TRUE) 
  var_mu2 = 10^2*diag(p)
  
  n_1 = 100
  S_1 = matrix(c(3^2, 0, 0, 1^2), nrow=p, ncol=p, byrow=TRUE)*(n_1-p-1)
  n_2 = 50
  S_2 = matrix(c(20^2,0,0,20^2), nrow=p, ncol=p, byrow=TRUE)*(n_2-p-1)
  
  g_alpha = 2
  h_alpha = 0.05
  g_beta = 2
  h_beta = 0.05
}
prior_list = list(mean_mu1=mean_mu1, var_mu1=var_mu1,
                  mean_mu2=mean_mu2, var_mu2=var_mu2,
                  n_1=n_1, n_2=n_2, S_1=S_1, S_2=S_2,
                  g_alpha=g_alpha, h_alpha=h_alpha, 
                  g_beta=g_beta, h_beta=h_beta,
                  D=2, K=2)

dat = read.csv("../rawdat_a.csv", stringsAsFactors=FALSE)
sbj = unique(dat$caseno)
crl = sort(sbj[grepl("C0.", sbj)])
pts = sort(sbj[grepl("P0.", sbj)])

# inputs for the inference function - which data file name, chain length, thinning etc.
inputs = list()
{
  input0 = list()
  input0$MCMCOut = 2000
  input0$MCMCBurnin = 2000
  input0$MCMCThin = 1
  input0$n.chains = 1
  for(chan in cord){
    outroot = chan
    inputs[[outroot]] = input0
    inputs[[outroot]]$chan = chan
    inputs[[outroot]]$prior_list = prior_list
  } # chans
}

ncores = 2
cl  = makeCluster(ncores) 
{
  clusterEvalQ(cl, {
    library("rstan")
    source("./helper_functions.R")
  })
  gmm_output = parLapply(cl, inputs, inference)
}
stopCluster(cl)

# save output
for(outroot in names(gmm_output)){
  output = gmm_output[[outroot]]
  
  write.table(output$post, file.path("./Output", folder, paste0(outroot,"_POST.txt")), row.names=FALSE, col.names=TRUE)
  write.table(output$classif, file.path("./Output", folder, paste0(outroot,"_CLASS.txt")), row.names=FALSE, col.names=TRUE)
}

### prior beliefs
output_prior = stan("./allData_GMM_prior.stan", data=prior_list, chains=1, 
                    iter=1e4, warmup=0, thin=1)

prior_all= as.matrix(output_prior)
prior = prior_all[, !grepl("_", colnames(prior_all))]

write.table(prior, file.path("./Output", folder, "PRIOR.txt"), 
            row.names=FALSE, col.names=TRUE)











