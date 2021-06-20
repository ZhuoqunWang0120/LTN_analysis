#!/usr/bin/env Rscript
pkgs=c('GetoptLong','philr','mvtnorm')
lapply(pkgs, require,character.only = TRUE)
niter=10000
nmc=10^6
GetoptLong(
  "SEED=i","random seed of simulation.",
  "modelCov=i","covariance scenario.",
  "niter=i","number of Gibbs iterations.",
  "nmc=i","number of MC samples in estimating clr covariance.",
  "lambda=f","lambda, 0 represents Gamma prior on it.",
  "WORK_DIR=s","working directory"
)

if (lambda==0){
  lambda='hp'
}
modelcovs=c('hub','block','sparse')
modelcov=modelcovs[modelCov]
source(paste0(WORK_DIR,"/src/utility/utility.R"))
source(paste0(WORK_DIR,"/src/experiments/covariance_estimation/gibbs.R"))
input_data=readRDS(paste0(WORK_DIR,"/cache/ps_sim.RData"))
tree=input_data$tree
datadir=paste0(WORK_DIR,'/cache/ln/')
result_dir=paste0(WORK_DIR,"/results/covariance_estimation/ln/")
filenam=paste0('i',SEED,'_',modelcov,'_lambda',lambda,'.RData')
yyl=readRDS(paste0(datadir,modelcov,'yyl_',SEED,'.RData'))
Y=yyl$Y
YL=yyl$YL
N=nrow(Y)
p=ncol(Y)
K=p+1
SSS=1
if (!file.exists(paste0(result_dir,filenam))){
  st1<-system.time(t<-try(gibbs1<-gibbs_glasso(niter=niter,YL=t(YL),Y=t(Y),r=1,s=0.01,SSS=1,lambda=lambda)))
  while("try-error" %in% class(t)) {
    SSS=SSS+1
    warning('numerical issues')
    t<-try(gibbs1<-gibbs_glasso(niter=niter,YL=t(YL),Y=t(Y),r=1,s=0.01,SSS=SSS,lambda=lambda))
  }
  samp_seq=(niter/2+1):niter
  MU=gibbs1$MU
  OMEGA=gibbs1$OMEGA
  mu=apply(MU[,samp_seq],1,mean)
  omega=Reduce('+',OMEGA[samp_seq])/length(OMEGA[samp_seq])
  sigma=solve(omega)
  rm(gibbs1)
  st2<-system.time(clrcov_ltn<-clrcov_sim_log(mu,sigma,tree,nmc,F,NULL))
  if (sum(is.na(clrcov_ltn))+sum(is.infinite(clrcov_ltn))==0){
    saveRDS(clrcov_ltn,paste0(result_dir,filenam))
  } else {
    warning('NA or Inf in clrcov')
  }
}
