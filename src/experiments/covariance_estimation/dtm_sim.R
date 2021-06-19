#!/usr/bin/env Rscript

argv=commandArgs(TRUE)
WORK_DIR=argv[1]
source(paste0(WORK_DIR,"/LTN_analysis/src/utility/utility.R"))
i=as.numeric(argv[2])
simdir=paste0(WORK_DIR,'/LTN_analysis/cache/dtm/')
system(paste0('mkdir ',simdir))
if (!file.exists(paste0(simdir,'sim',i,'.RData'))){
  nsim=200
  total=10^5
  input_data1=readRDS(paste0(WORK_DIR,"/LTN_analysis/cache/dtm_diab_MoM.RData"))
  input_data2=readRDS(paste0(WORK_DIR,"/LTN_analysis/cache/ps_sim.RData"))
  theta=as.vector(input_data1[[1]])
  tau=as.vector(input_data1[[2]])
  tree=input_data2$tree
  set.seed(i)
  sim_i=dtm_sim(nsim,tree,theta,tau,total)
  saveRDS(sim_i,paste0(simdir,'sim',i,'.RData'))
}
