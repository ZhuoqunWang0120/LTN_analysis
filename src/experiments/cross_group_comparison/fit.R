#!/usr/bin/env Rscript
library(GetoptLong)

c0 = 1
d0 = 0.001
c1 = 1
d1 = 0.001
nu = 1
niter=10000
adjust=T
gm=100
pnull=0.5
r=0
a1=1
a2=1
a1a2='none'

GetoptLong(
  "r=i","number of latent factors (unwarranted experimental feature)",
  "c0=f","hyperparameter c0.",
  "d0=f","hyperparameter d0.",
  "c1=f","hyperparameter c1.",
  "d1=f","hyperparameter d1.",
  "nu=i","hyperparameter nu.",
  "niter=i","number of iterations.",
  "reff!","whether to include random effects in addition to factors.",
  "reffcov=i","random effect covariance, 1: diagonal, 2: sparse.",
  "pi_only!","whether to only save gibbs samples of pi.",
  "gm=f","m in gprior.",
  "pnull=f","prior probability of joint H0.",
  "a1=f","hyperparameter a1.",
  "a2=f","hyperparameter a2.",
  "a1a2=s","whether to put hyperprior on a1,a2.",
  "h=i","null or alternative.",
  "lambda=f","lambda, 0 represents Gamma prior on it.",
  "i=i","seed of simulation.",
  "scenario=s","simulation scenario.",
  "WORK_DIR=s","working directory."
)


source(paste0(WORK_DIR,"/src/utility/mixed_effects.R"))
result_dir=paste0(WORK_DIR,"/results/cross_group_comparison/",scenario,"/")
model_index=reffcov
gprior_m=gm
gibbsdir=paste0(result_dir,'/gibbs/')
pidir=paste0(result_dir,'/pi/')
teststatdir=paste0(result_dir,'/teststat/')
sdir=paste0(result_dir,'/SSS/')
randmodel=c('norand','diagonal','sparse')[model_index+1]
try(system(paste0('mkdir -p ',result_dir)))
try(system(paste0('mkdir -p ',gibbsdir)))
try(system(paste0('mkdir -p ',pidir)))
try(system(paste0('mkdir -p ',teststatdir)))
try(system(paste0('mkdir -p ',sdir)))
# filenam=paste0('i',i,'r',r,'H',h,'gm',gprior_m,randmodel,'nu',nu,'a1',a1,'a2',a2,'lambda',lambda,'iter',niter,'.RData')
filenam=paste0('i',i,'H',h,'_',randmodel,'_lambda',lambda,'.RData')
if (!file.exists(paste0(result_dir,'/teststat/',filenam))){
  input_data1=readRDS(paste0(WORK_DIR,"/cache/ps_sim.RData"))
  input_data2=readRDS(paste0(WORK_DIR,"/cache/cross_group_comparison/",scenario,"/sim",i,"H",h,".RData"))
  Y=input_data2$Y
  YL=input_data2$YL
  N=nrow(Y)
  p=ncol(Y)
  K=p+1
  g=input_data1$g
  Xtest=input_data2$Xtest
  Xadjust=input_data1$Xadjust
  grouplabel=input_data1$grouplabel
  SSS=1
  st<-system.time(t<-try(gibbs<-gibbs_crossgroup(N=N,p=p,g=g,r=r,YL=YL,Y=Y,Xtest=Xtest,Xadjust=Xadjust,grouplabel=grouplabel,c0=c0,d0=d0,c1=c1,d1=d1,nu=nu,niter=niter,adjust=adjust,reff=reff,reffcov=reffcov,SEED=SSS,pi_only=pi_only,gprior_m=gm,pnull=pnull,a1=a2,a2=a2,a1a2=a1a2,lambda_fixed=lambda)))
  while("try-error" %in% class(t)) {
    SSS=SSS+1
    warning('numerical issues')
    t<-try(gibbs<-gibbs_crossgroup(N=N,p=p,g=g,r=r,YL=YL,Y=Y,Xtest=Xtest,Xadjust=Xadjust,grouplabel=grouplabel,c0=c0,d0=d0,c1=c1,d1=d1,nu=nu,niter=niter,adjust=adjust,reff=reff,reffcov=reffcov,SEED=SSS,pi_only=pi_only,gprior_m=gm,pnull=pnull,a1=a2,a2=a2,a1a2=a1a2,lambda_fixed=lambda))
  }
  saveRDS(gibbs,paste0(gibbsdir,'/',filenam))
  BETA1=lapply(gibbs$BETA1,function(x){as.vector(x)})
  BETAMAT1=do.call(rbind,BETA1)
  post_pi=matrix(apply(BETAMAT1[(niter/2+1):niter,],2,function(x){sum(x!=0)})/(niter/2),nrow = 1)
  saveRDS(post_pi,paste0(pidir,'/',filenam))
  all0=rowSums(BETAMAT1[(niter/2+1):niter,]!=0)
  all0freq=sum(all0==0)/(niter/2)
  sum0=1-mean(all0)/p
  saveRDS(list(pjap=1-all0freq,sum0=sum0),paste0(teststatdir,'/',filenam))
  if (SSS>1){
    saveRDS(SSS,paste0(sdir,'/',filenam))}
}
