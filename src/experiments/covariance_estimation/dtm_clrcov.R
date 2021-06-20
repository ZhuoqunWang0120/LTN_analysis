#!/usr/bin/env Rscript

argv=commandArgs(TRUE)
WORK_DIR=argv[1]
nmc=as.numeric(argv[2])
source(paste0(WORK_DIR,"/src/utility/utility.R"))
input_data1=readRDS(paste0(WORK_DIR,"/cache/dtm_diab_MoM.RData"))
input_data2=readRDS(paste0(WORK_DIR,"/cache/ps_sim.RData"))
theta=as.vector(input_data1[[1]])
tau=as.vector(input_data1[[2]])
tree=input_data2$tree
clrcov_true=clrcov_dtm_sim_log(nmc,tree,theta,tau,SSS=1,savesamp=F,dir=NULL)
system(paste0('mkdir -p ',WORK_DIR,"/cache/dtm/"))
if (sum(is.na(clrcov_true))+sum(is.infinite(clrcov_true))==0){
  saveRDS(clrcov_true,paste0(WORK_DIR,"/cache/dtm/clrcov_true.RData"))
  print('finished calculating the true clr covariance in DTM example')
} else{
  warning('NA in clr covariance')
}