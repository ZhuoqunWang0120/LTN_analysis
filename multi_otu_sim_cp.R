#!/usr/bin/env Rscript
argv=commandArgs(TRUE)
i=as.numeric(argv[1])
n_or_a=as.numeric(argv[2])
source('/work/zw122/pg_tree/cross_group/rscripts/utilities.R')
if (!file.exists(paste0('/work/zw122/pg_tree/cross_group/diab/multi_otu/sim',i,'H',n_or_a,'.RData'))){
input_data=readRDS('/work/zw122/pg_tree/cross_group/diab/para/diab_input_data_otu.RData')
cnt=input_data$cnt
tree=input_data$tree
K=ncol(cnt)
N=nrow(cnt)
set.seed(i)
group2=sample(N,ceiling(N/2),replace = F)
otu20=input_data$top20_ra
if (n_or_a==1){
  otu_modify=sample(otu20,8)
  cnt[group2,otu_modify]=ceiling(1.5*cnt[group2,otu_modify])
}
Xtest=matrix(0,N,1)
Xtest[group2,]=1
yyl=seqtab2y(cnt,tree)
saveRDS(list(Y=yyl$Y,YL=yyl$YL,cnt=cnt,Xtest=Xtest,group2=group2),paste0('/work/zw122/pg_tree/cross_group/diab/multi_otu/sim',i,'H',n_or_a,'.RData'))
}

