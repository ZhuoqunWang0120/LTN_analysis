#!/usr/bin/env Rscript


argv=commandArgs(TRUE)
WORK_DIR=argv[1]
lambda=argv[2]
nsim1=as.numeric(argv[3])
source(paste0(WORK_DIR,"/src/utility/utility.R"))
if (lambda=='0'){lambda='hp'}
risk=matrix(-1,4,4,dimnames = list(c("Frobenius","L1","L_inf","spectral"),c("LN-hub","LN-block","LN-sparse","DTM")))
# DTM
result_dir=paste0(WORK_DIR,"/results/covariance_estimation/dtm/")
clrcov0=readRDS(paste0(WORK_DIR,"/cache/dtm/clrcov_true.RData"))
files=grep(paste0('lambda',lambda,'.RData'),list.files(result_dir,full.names = T),value = T)
loss=do.call(rbind,lapply(files, function(x){matloss(cov2cor(readRDS(x))-cov2cor(clrcov0))}))
risk[,'DTM']=(colSums(loss)/nrow(loss))
# LN
result_dir=paste0(WORK_DIR,"/results/covariance_estimation/ln/")
for (m in c('hub','block','sparse')){
  files=grep(paste0(m,'_lambda',lambda,'.RData'),list.files(result_dir,full.names = T),value = T)
  loss=do.call(rbind,lapply(1:nsim1, function(x){
    clrcov0=readRDS(paste0(WORK_DIR,"/cache/ln/",m,"clrcov_",x,".RData"))
    clrcov1=readRDS(grep(paste0('i',x,'_'),files,value=T))
    matloss(cov2cor(clrcov1)-cov2cor(clrcov0))}))
  risk[,paste0("LN-",m)]=(colSums(loss)/nrow(loss))
}
if (lambda=='hp'){lambda='GammaPrior'}
system(paste0('mkdir -p ',WORK_DIR,"/results/covariance_estimation/risk/"))
saveRDS(risk,paste0(WORK_DIR,"/results/covariance_estimation/risk/lambda_",lambda,".RData"))


