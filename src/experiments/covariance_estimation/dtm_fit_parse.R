#!/usr/bin/env Rscript

library(GetoptLong,quietly = T)
library(mvtnorm,quietly = T)
GetoptLong(
  "i=i","random seed of simulation.",
  "niter=i","number of Gibbs iterations.",
  "nmc=i","number of MC samples in estimating clr covariance.",
  "lambda=f","lambda, 0 represents Gamma prior on it.",
  "WORK_DIR=s","working directory."
)
if (lambda==0){
  lambda='hp'
}
filenam=paste0('i',i,'lambda',lambda,'.RData')
source(paste0(WORK_DIR,"/src/utility/utility.R"))
source(paste0(WORK_DIR,"/src/experiments/covariance_estimation/gibbs.R"))
result_dir=paste0(WORK_DIR,"/results/covariance_estimation/dtm/")
system(paste0('mkdir ',result_dir,'gibbs/'))
system(paste0('mkdir ',result_dir,'clrcov/'))
system(paste0('mkdir ',result_dir,'SSS/'))
datadir=paste0(WORK_DIR,'/cache/dtm/')
if (!file.exists(paste0(result_dir,'gibbs/',filenam))){
  input_data=readRDS(paste0(WORK_DIR,"/cache/ps_sim.RData"))
  tree=input_data$tree
  dat_i=readRDS(paste0(datadir,'sim',i,'.RData'))
  yyl=dat_i[1:2]
  Y=yyl$Y
  YL=yyl$YL
  N=nrow(Y)
  p=ncol(Y)
  K=p+1
  SSS=1
  st1<-system.time(t<-try(gibbs1<-gibbs_glasso(niter=niter,YL=t(YL),Y=t(Y),r=1,s=0.01,SSS=1,lambda=lambda)))
  while("try-error" %in% class(t)) {
    SSS=SSS+1
    warning('error')
    t<-try(gibbs1<-gibbs_glasso(niter=niter,YL=t(YL),Y=t(Y),r=1,s=0.01,SSS=SSS,lambda=lambda))
  }
  if (SSS>1){
    saveRDS(SSS,paste0(result_dir,'SSS/',filenam))
  }
  saveRDS(gibbs1,paste0(result_dir,'gibbs/',filenam))
  samp_seq=ceiling(niter/2):niter
  MU=gibbs1$MU
  OMEGA=gibbs1$OMEGA
  mu=apply(MU[,samp_seq],1,mean)
  omega=Reduce('+',OMEGA[samp_seq])/length(OMEGA[samp_seq])
  sigma=solve(omega)
  rm(gibbs1)
  st2<-system.time(clrcov_ltn<-clrcov_sim_log(mu,sigma,tree,nmc,F,NULL))
  if (sum(is.na(clrcov_ltn))+sum(is.infinite(clrcov_ltn))==0){
    saveRDS(clrcov_ltn,paste0(result_dir,'clrcov/',filenam))
  } else {
    warning('NA or Inf in clrcov')
  }
}


