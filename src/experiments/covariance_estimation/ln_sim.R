#!/usr/bin/env Rscript
library(philr)
library(mvtnorm)
argv=commandArgs(TRUE)
WORK_DIR=argv[1]
SEED=as.numeric(argv[2])
modelCov=as.numeric(argv[3])
#muvar=as.numeric(argv[3])
muvar=16
modelcovs=c('hub','block','sparse')
modelcov=modelcovs[modelCov]
source(paste0(WORK_DIR,"/src/experiments/covariance_estimation/ln_graph.R"))
source(paste0(WORK_DIR,"/src/utility/utility.R"))
input_data=readRDS(paste0(WORK_DIR,"/cache/ps_sim.RData"))
tree=input_data$tree
datadir=paste0(WORK_DIR,"/cache/ln/")
system(paste0('mkdir -p ',datadir))
set.seed(SEED)
parasimu=generateDataCell(n=200,modelCov=modelCov,p1=49,p2=100,p3=200,nRep=1)
mu=rnorm(49,0,sqrt(muvar))
saveRDS(list(mu=mu,sigma=solve(parasimu$sigma1)),paste0(datadir,modelcov,'para_',SEED,'.RData'))
# use sigma1 as omega!
eta=rmvnorm(200,mu,solve(parasimu$sigma1))
colnames(eta)=paste0('n',1:49)
tree$node.label=paste0('n',1:49)
philinv=philrInv(df.ilrp = eta,tree=tree,V=NULL)
cnt=t(apply(philinv,1,function(x){rmultinom(1,10^5,x)}))
colnames(cnt)=tree$tip.label
saveRDS(cnt,paste0(datadir,modelcov,'cnt_',SEED,'.RData'))
yyl=seqtab2y(cnt,tree)
Y=yyl$Y
YL=yyl$YL
saveRDS(list(Y=Y,YL=YL),paste0(datadir,modelcov,'yyl_',SEED,'.RData'))

# clrcov via transformation
sig=solve(parasimu$sigma1)
K=50
p=K-1
ilrE=diag(p)
colnames(ilrE)=tree$node.label
E=philrInv(ilrE,tree)
H=t(apply(E, 1, function(x){clrp(x,rep(1,K))})) # clrp=clr when p=1
G=t(H)%*%sig%*%H
saveRDS(G,paste0(datadir,modelcov,'clrcov_',SEED,'.RData'))
