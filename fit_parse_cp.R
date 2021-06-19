#!/usr/bin/env Rscript
library(GetoptLong)

c0 = 1
d0 = 0.001
c1 = 1
d1 = 0.001
niter=10000
adjust=T
gm=100
pnull=0.5

GetoptLong(
  "r=i","number of factors.",
  "c0=f","hyperparameter c0.",
  "d0=f","hyperparameter d0.",
  "c1=f","hyperparameter c1.",
  "d1=f","hyperparameter d1.",
  "nu=i","hyperparameter nu.",
  "niter=i","number of iterations.",
  "reff!","whether to include random effects in addition to factors.",
  "reffcov=i","random effect covariance.",
  "pi_only!","whether to only save gibbs samples of pi.",
  "gm=f","m in gprior.",
  "pnull=f","prior prob of joint H0.",
  "a1=f","hyperparameter a1.",
  "a2=f","hyperparameter a2.",
  "a1a2=s","whether to put hyperprior on a1,a2.",
  "h=i","null or alternative.",
  "lambda=f","lambda, 0 represents Gamma prior on it.",
  "i=i","seed of simulation.",
  "scenario=s","simulation scenario."
)


source('/work/zw122/pg_tree/cross_group/rscripts/cross_group_comparison_unified.R')
resdir=paste0('/work/zw122/pg_tree/cross_group/results/diab/',scenario,'/')
n_or_a=h
model_index=reffcov
gprior_m=gm
gibbsdir=paste0(resdir,'/gibbs/')
pidir=paste0(resdir,'/pi/')
teststatdir=paste0(resdir,'/teststat/')
sdir=paste0(resdir,'/SSS/')
randmodel=c('norand','diagonal','sparse')[model_index+1]
try(system(paste0('mkdir ',resdir)))
try(system(paste0('mkdir ',gibbsdir)))
try(system(paste0('mkdir ',pidir)))
try(system(paste0('mkdir ',teststatdir)))
try(system(paste0('mkdir ',sdir)))
filenam=paste0('i',i,'r',r,'H',n_or_a,'gm',gprior_m,randmodel,'nu',nu,'a1',a1,'a2',a2,'lambda',lambda,'iter',niter,'.RData')
if (!file.exists(paste0(resdir,'/teststat/',filenam)) & !file.exists(paste0(resdir,'/teststat/lambda_auc/',filenam))){
  input_data1=readRDS('/work/zw122/pg_tree/cross_group/diab/para/diab_input_data_otu.RData')
  input_data2=readRDS(paste0('/work/zw122/pg_tree/cross_group/diab/',scenario,'/sim',i,'H',n_or_a,'.RData'))
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
    print('error!')
    print(SSS)
    t<-try(gibbs<-gibbs_crossgroup(N=N,p=p,g=g,r=r,YL=YL,Y=Y,Xtest=Xtest,Xadjust=Xadjust,grouplabel=grouplabel,c0=c0,d0=d0,c1=c1,d1=d1,nu=nu,niter=niter,adjust=adjust,reff=reff,reffcov=reffcov,SEED=SSS,pi_only=pi_only,gprior_m=gm,pnull=pnull,a1=a2,a2=a2,a1a2=a1a2,lambda_fixed=lambda))
  }
  try(print(st))
  saveRDS(gibbs,paste0(gibbsdir,'/',filenam))
  BETA1=lapply(gibbs$BETA1,function(x){as.vector(x)})
  BETAMAT1=do.call(rbind,BETA1)
  post_pi=matrix(apply(BETAMAT1[(niter/2+1):niter,],2,function(x){sum(x!=0)})/(niter/2),nrow = 1)
  post_pi
  saveRDS(post_pi,paste0(pidir,'/',filenam))
  all0=rowSums(BETAMAT1[(niter/2+1):niter,]!=0)
  all0freq=sum(all0==0)/(niter/2)
  sum0=1-mean(all0)/p
  print(all0freq)
  print(sum0)
  saveRDS(list(all0freq=all0freq,sum0=sum0),paste0(teststatdir,'/',filenam))
  if (SSS>1){
    saveRDS(SSS,paste0(sdir,'/',filenam))}
}
