#!/usr/bin/env Rscript

library(ROCR)
argv=commandArgs(TRUE)
WORK_DIR=argv[1]
lambda=argv[2]
for (s in c('single_otu','multi_otu')){
  system(paste0('mkdir -p ',WORK_DIR,'/results/cross_group_comparison/',s,'/summary/'))
  for (m in c('sparse','diagonal')){
    result_dir=paste0(WORK_DIR,'/results/cross_group_comparison/',s,'/pjap/')
    files0=grep(paste0('H0_',m,'_lambda',lambda,'.RData'),list.files(result_dir,full.names = T),value = T)
    files1=grep(paste0('H1_',m,'_lambda',lambda,'.RData'),list.files(result_dir,full.names = T),value = T)
    pjaps=c(sapply(files0, function(x){readRDS(x)[[1]]}),sapply(files1, function(x){readRDS(x)[[1]]}))
    labels=c(rep(0,length(files0)),rep(1,length(files1)))
    if (length(unique(labels))>1){
      if (lambda=='0' & m=='sparse'){lambda='GammaPrior'}
      if (m=='diagonal'){lambda='NA'}
      saveRDS(data.frame(pjap=pjaps,label=labels,row.names = NULL),paste0(WORK_DIR,'/results/cross_group_comparison/',s,'/summary/',m,'_lambda',lambda,'.RData'))
      png(paste0(WORK_DIR,'/results/cross_group_comparison/',s,'/summary/',m,'_lambda',lambda,'.png'))
      plot(performance(prediction(pjaps,labels),'tpr','fpr'),main='ROC curve')
      dev.off()
    }
  }
}