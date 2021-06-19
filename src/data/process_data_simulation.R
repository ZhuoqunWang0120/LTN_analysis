#!/usr/bin/env Rscript

argv=commandArgs(TRUE)
WORK_DIR=argv[1]
source(paste0(WORK_DIR,"/src/utility/utility.R"))
library(phyloseq)
library(ggplot2)
library(reshape2)
library(gridExtra)
psotu50=readRDS(paste0(WORK_DIR,"/cache/ps_otu50.RData"))
md_16S=sample_data(psotu50)
# for simulation, keep finland samples in control group, only adjust for age
samples_finland_control=rownames(md_16S)[md_16S$Case_Control=='control' & md_16S$Country=='Finland']
ps_finland_control=prune_samples(samples_finland_control,psotu50)
tree=phy_tree(ps_finland_control)
cnt=matrix(otu_table(ps_finland_control),ncol=ncol(otu_table(ps_finland_control)),dimnames = dimnames(otu_table(ps_finland_control)))
sampdat=sample_data(ps_finland_control)
age=sampdat$Age_at_Collection
age_s=(age-mean(age))/sd(age)
Xadjust=matrix(age_s,ncol=1)
Xadjust=cbind(1,Xadjust)
grouplabel=as.numeric(sampdat$Subject_ID)
g=max(grouplabel)
gprior_m=1
tree=as.Node(phy_tree(ps_finland_control))
tree$Set(name=51:99,traversal = 'pre-order',filterFun = isNotLeaf)
tree=as.phylo(tree)
ra50=apply(cnt, 1, function(x){x/sum(x)})
top20_ra=names(sort(rowSums(ra50),decreasing = T)[1:20])
saveRDS(list(cnt=cnt,Xadjust=Xadjust,grouplabel=grouplabel,g=g,gprior_m=gprior_m,tree=tree,top20_ra=top20_ra),paste0(WORK_DIR,"/cache/ps_sim.RData"))
