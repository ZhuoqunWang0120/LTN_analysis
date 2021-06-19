#!/usr/bin/env Rscript

source('/work/zw122/pg_tree/cross_group/rscripts/utilities.R')
source('/work/zw122/pg_tree/cross_group/rscripts/cross_group_comparison_unified.R')

argv=commandArgs(TRUE)
r=0
niter=as.numeric(argv[1])
covariate_index=as.numeric(argv[2])
model_index=as.numeric(argv[3]) #1 diagonal, 2 sparse
lambda_fixed=as.numeric(argv[4])
gprior_m=100
modelspec=c('diagonal','sparse')[model_index]
covariates=c("BF","Solid_Food","Eggs","Fish","Soy_Prod","Rye","Barley","Buckwheat_Millet","casecontrol")
covariate=covariates[covariate_index]
dietcovariates=c("BF","Solid_Food","Eggs","Fish","Soy_Prod","Rye","Barley","Buckwheat_Millet")
datadir='/work/zw122/pg_tree/cross_group/diab/para/application50/data/'
resdir=paste0('/work/zw122/pg_tree/cross_group/results/diab/application/application50/')

  system(paste0('mkdir ',resdir))
  system(paste0('mkdir ',resdir,'/gibbs'))
  system(paste0('mkdir ',resdir,'/pi'))
  system(paste0('mkdir ',resdir,'/teststat'))
  system(paste0('mkdir ',resdir,'/SSS'))

input_data=readRDS(paste0(datadir,'application_',covariate,'.RData'))
cnt=input_data$cnt
tree=input_data$tree
yyl=seqtab2y(cnt,tree)
Y=yyl$Y
YL=yyl$YL
N=nrow(Y)
p=ncol(Y)
K=p+1
g=input_data$g
Xtest=input_data$Xtest
Xadjust=input_data$Xadjust
grouplabel=input_data$grouplabel
SSS=1
t=NULL
st=NULL
filenam=paste0('r',r,'iter',niter,covariate,modelspec,'lambda',lambda_fixed,'.RData')
if (!file.exists(paste0(resdir,'/teststat/',filenam))){st<-system.time(t<-try(gibbs<-gibbs_crossgroup(N=N,p=p,g=g,r=r,YL=YL,Y=Y,Xtest=Xtest,Xadjust=Xadjust,grouplabel=grouplabel,c0=1,d0=0.001,c1=1,d1=0.001,niter=niter,adjust=T,reff=T,SEED=SSS,pi_only=F,gprior_m=gprior_m,reffcov=model_index,lambda_fixed=lambda_fixed)))
while("try-error" %in% class(t)) {
  warning('seed ',SSS,' failed')
  SSS=SSS+1
  print('error!')
  print(SSS)
  t<-try(gibbs<-gibbs_crossgroup(N=N,p=p,g=g,r=r,YL=YL,Y=Y,Xtest=Xtest,Xadjust=Xadjust,grouplabel=grouplabel,c0=1,d0=0.001,c1=1,d1=0.001,niter=niter,adjust=T,reff=T,SEED=SSS,pi_only=F,gprior_m=gprior_m,reffcov=model_index,lambda_fixed=lambda_fixed,SEED=SSS))
}
try(print(st))
saveRDS(gibbs,paste0(resdir,'/gibbs/',filenam))
BETA1=lapply(gibbs$BETA1,function(x){as.vector(x)})
BETAMAT1=do.call(rbind,BETA1)
post_pi=matrix(apply(BETAMAT1[(niter/2):niter,],2,function(x){sum(x!=0)})/(niter/2+1),nrow = 1)
post_pi
saveRDS(post_pi,paste0(resdir,'/pi/',filenam))
all0=rowSums(BETAMAT1[(niter/2):niter,]!=0)
all0freq=sum(all0==0)/(niter/2+1)
sum0=1-mean(all0)/p
print(all0freq)
print(sum0)
saveRDS(list(all0freq=all0freq,sum0=sum0),paste0(resdir,'/teststat/',filenam))
warnings()
if (SSS>1){
saveRDS(SSS,paste0(resdir,'/SSS/',filenam))
}
}
