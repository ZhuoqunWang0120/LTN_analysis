#!/usr/bin/env Rscript


argv=commandArgs(TRUE)
WORK_DIR=argv[1]
source(paste0(WORK_DIR,"/src/utility/utility.R"))
source(paste0(WORK_DIR,"/src/utility/mixed_effects.R"))
r=0
niter=as.numeric(argv[2])
covariate_index=as.numeric(argv[3])
model_index=as.numeric(argv[4]) 
lambda_fixed=as.numeric(argv[5])
gprior_m=100
modelspec=c('diagonal','sparse')[model_index]
covariates=c("BF","Solid_Food","Eggs","Fish","Soy_Prod","Rye","Barley","Buckwheat_Millet","casecontrol")
covariate=covariates[covariate_index]
dietcovariates=c("BF","Solid_Food","Eggs","Fish","Soy_Prod","Rye","Barley","Buckwheat_Millet")
datadir=paste0(WORK_DIR,"/cache/")
resdir=paste0(WORK_DIR,"/results/application/")
system(paste0('mkdir -p ',resdir,'/pmap'))
system(paste0('mkdir -p ',resdir,'/pjap'))
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
filenam=paste0(covariate,'_',modelspec,'_lambda',lambda_fixed,'.RData')
if (!file.exists(paste0(resdir,'/pjap/',filenam))){
  st<-system.time(t<-try(gibbs<-gibbs_crossgroup(N=N,p=p,g=g,r=r,YL=YL,Y=Y,Xtest=Xtest,Xadjust=Xadjust,grouplabel=grouplabel,c0=1,d0=0.001,c1=1,d1=0.001,niter=niter,adjust=T,reff=T,SEED=SSS,pi_only=F,gprior_m=gprior_m,reffcov=model_index,lambda_fixed=lambda_fixed)))
while("try-error" %in% class(t)) {
  SSS=SSS+1
  warning('numerical issues')
  t<-try(gibbs<-gibbs_crossgroup(N=N,p=p,g=g,r=r,YL=YL,Y=Y,Xtest=Xtest,Xadjust=Xadjust,grouplabel=grouplabel,c0=1,d0=0.001,c1=1,d1=0.001,niter=niter,adjust=T,reff=T,SEED=SSS,pi_only=F,gprior_m=gprior_m,reffcov=model_index,lambda_fixed=lambda_fixed,SEED=SSS))
}
BETA1=lapply(gibbs$BETA1,function(x){as.vector(x)})
BETAMAT1=do.call(rbind,BETA1)
post_pi=matrix(apply(BETAMAT1[(niter/2):niter,],2,function(x){sum(x!=0)})/(niter/2+1),nrow = 1)
saveRDS(post_pi,paste0(resdir,'/pmap/',filenam))
all0=rowSums(BETAMAT1[(niter/2):niter,]!=0)
all0freq=sum(all0==0)/(niter/2+1)
sum0=1-mean(all0)/p
saveRDS(list(pjap=1-all0freq,sum0=sum0),paste0(resdir,'/pjap/',filenam))
}
