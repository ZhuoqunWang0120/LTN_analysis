#!/usr/bin/env Rscript
library(ape)
library(data.tree)
library(phyloseq)
argv=commandArgs(TRUE)
WORK_DIR=argv[1]
dietcovariates=c("BF","Solid_Food","Eggs","Fish","Soy_Prod","Rye","Barley","Buckwheat_Millet")
ps=readRDS(paste0(WORK_DIR,"/LTN_analysis/cache/ps_otu50.RData"))
cnt=otu_table(ps)
tree=phy_tree(ps)
tree=as.Node(tree)
tree$Set(name=51:99,traversal = 'pre-order',filterFun = isNotLeaf)
tree=ape::as.phylo(tree)
sampdat=sample_data(ps)
grouplabel=as.numeric(sampdat$Subject_ID)
g=max(grouplabel)
for (i in seq_along(dietcovariates)){
  covariate=dietcovariates[i]
  Xtest=matrix(as.numeric(sampdat[,covariate]=='true'),ncol=1)
  Xadjust=model.matrix(as.formula(paste0('~Age_at_Collection+Country+T1D_Status+Gender+',paste(dietcovariates[-i],collapse = '+'))),data.frame(sampdat))
  if (sum(is.na(Xadjust))==0 & det(t(Xadjust)%*%Xadjust)!=0){
    saveRDS(list(cnt=cnt,Xtest=Xtest,Xadjust=Xadjust,grouplabel=grouplabel,g=g,tree=tree),paste0(WORK_DIR,'/LTN_analysis/cache/application_',covariate,'.RData'))
  } 
  else{print(c(covariate,' collinearity'))}
}

Xtest=matrix(as.numeric(sampdat[,'Case_Control']=='case'),ncol=1)
Xadjust=model.matrix(as.formula(paste0('~Age_at_Collection+Country+Gender+',paste(dietcovariates,collapse = '+'))),data.frame(sampdat))
if (sum(is.na(Xadjust))==0 & det(t(Xadjust)%*%Xadjust)!=0){
  saveRDS(list(cnt=cnt,Xtest=Xtest,Xadjust=Xadjust,grouplabel=grouplabel,g=g,tree=tree),paste0(WORK_DIR,'/LTN_analysis/cache/application_','casecontrol','.RData'))
} 

