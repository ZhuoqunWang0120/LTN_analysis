#!/usr/bin/env Rscript

argv=commandArgs(TRUE)
WORK_DIR=argv[1]
source(paste0(WORK_DIR,"/LTN_analysis/src/utility/utility.R"))
if(!file.exists(paste0(WORK_DIR,"/LTN_analysis/cache/dtm_diab_MoM.RData"))){
  inputdata=readRDS(paste0(WORK_DIR,"/LTN_analysis/cache/ps_sim.RData"))
  cnt=inputdata$cnt
  tree=inputdata$tree
  est=dtm_mom(cnt,tree)
  saveRDS(est,paste0(WORK_DIR,"/LTN_analysis/cache/dtm_diab_MoM.RData"))
}
