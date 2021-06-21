#!/usr/bin/env Rscript

library(data.tree)
library(phyloseq)
argv=commandArgs(TRUE)
WORK_DIR=argv[1]
lambda=argv[2]
covariates=c("BF","Solid_Food","Eggs","Fish","Soy_Prod","Rye","Barley","Buckwheat_Millet","casecontrol")
result_dir=paste0(WORK_DIR,"/results/application/")
ps50=readRDS(paste0(WORK_DIR,"/cache/ps_otu50.RData"))
tree=phy_tree(ps50)
pjapmat=matrix(-1,9,1,dimnames = list(c('Breastfeeding',covariates[2:8],'T1D'),'PJAP'))
system(paste0('mkdir -p ',result_dir,'/PMAP_plot/'))
system(paste0('mkdir -p ',result_dir,'/node_table/'))
system(paste0('mkdir -p ',result_dir,'/PJAP_table/'))
threshold_node=function(pmap,threshold){
  return(which(pmap>threshold)+50)
}
tree_copy=tree
tree_copy$tip.label=1:50
tree_copy$node.label=51:99
tn=as.Node(tree_copy)
tn$Do(function(x) {
  leftchild=names(x$children[1])
  rightchild=names(x$children[2])
  x$left <-leftchild
  x$right <-rightchild
}, traversal = "post-order",filterFun = isNotLeaf)
left=tn$Get('left', traversal = "pre-order", filterFun=isNotLeaf)
right=tn$Get('right', traversal = "pre-order", filterFun=isNotLeaf)
nodedf=data.frame(node=51:99,left=left,right=right,row.names = NULL)
nodepaths=nodepath(tree_copy)
taxtab50=tax_table(ps50)
descendant_otu=function(node){
  leftchild=nodedf[nodedf$node==node,'left']
  rightchild=nodedf[nodedf$node==node,'right']
  leftotu=unlist(lapply(nodepaths[which(unlist(lapply(nodepaths, function(x){leftchild%in%x})))],function(x){x[length(x)]}))
  rightotu=unlist(lapply(nodepaths[which(unlist(lapply(nodepaths, function(x){rightchild%in%x})))],function(x){x[length(x)]}))
  return(list(leftotu=leftotu,rightotu=rightotu))
}
common_diff=function(taxtableft,taxtabright,taxlevel=T,showpath=F){
  taxtabmerge=rbind(taxtableft,taxtabright)
  level_common=max(which(apply(rbind(taxtableft,taxtabright), 2, function(x){length(unique(x))})==1))
  if (showpath){
    taxa_common=paste0(taxtableft[1,1:level_common],collapse = '->')
  }else{
    names_split=sapply(taxtableft[1,1:level_common], function(x){strsplit(x,'__')[[1]][2]})
    taxa_common=names(names_split[max(which(!is.na(names_split)))])
  }
  if (!taxlevel){
    taxa_common=sapply(taxa_common,function(x){strsplit(x,'__')[[1]][2]})
  }
  if (level_common<ncol(taxtab50)){
    leftvec=as.character(unique(taxtableft[,level_common+1]))
    rightvec=as.character(unique(taxtabright[,level_common+1]))
    if (!taxlevel){
      leftvec=sapply(leftvec, function(x){strsplit(x,'__')[[1]][2]})
      rightvec=sapply(rightvec, function(x){strsplit(x,'__')[[1]][2]})
    }
    lefttaxa=paste0(leftvec,collapse = ', ')
    righttaxa=paste0(rightvec,collapse=', ')
  }else{
    lefttaxa=paste0(rownames(taxtableft),collapse = ', ')
    righttaxa=paste0(rownames(taxtabright),collapse=', ')
  }
  return(list(common=taxa_common,left=lefttaxa,right=righttaxa))
}
for (i in 1:9){
  covariate=covariates[i]
  filenam=paste0(covariate,'_sparse_lambda',lambda,'.RData')
  # PJAP
  pjapmat[i,1]=readRDS(paste0(result_dir,'/pjap/',filenam))[[1]]
  # PMAP plot
  png(paste0(result_dir,'/PMAP_plot/PMAP_',rownames(pjapmat)[i],'_lambda',lambda,'.png'))
  par(mar=c(1,1,1,1))
  layout(mat = matrix(1:2,ncol = 1),heights = c(85,15))
  col_Pal = colorRampPalette(c('white', 'red'))
  pmap=readRDS(paste0(result_dir,'/pmap/',filenam))
  node_col = col_Pal(500)[as.numeric(cut(c((pmap),0,1), breaks = 500)) ]
  plot(tree,main=paste0('PMAPs of ',rownames(pjapmat)[i]),show.node.label=FALSE,direction="downwards",show.tip.label = FALSE)
  nodelabels(bg=node_col,frame='none',cex=2,pch=21,col ='black')
  nodelabels(text=1:49,frame='none',cex=1)
  legend_image <- as.raster(matrix(col_Pal(500), nrow=1))
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
  text(x=seq(0,1,l=5), y =0.4 , labels = seq(0,1,l=5),cex=1)
  rasterImage(legend_image, 0, 0.5, 1,0.7)
  dev.off()
  # node-level table
  nodes=threshold_node(pmap,0.5)
  otulist=lapply(nodes,function(x){descendant_otu(x)})
  taxtablist=lapply(otulist,function(x){list(taxtab50[x$leftotu,],taxtab50[x$rightotu,])})
  nodetaxa=lapply(taxtablist,function(x){common_diff(x[[1]],x[[2]],taxlevel=F)})
  common=unlist(lapply(nodetaxa, function(x){x[[1]]}))
  left=unlist(lapply(nodetaxa, function(x){x[[2]]}))
  right=unlist(lapply(nodetaxa, function(x){x[[3]]}))
  markdf=data.frame(taxon_in_common=common,taxa_left=left,taxa_right=right)
  postmean=readRDS(paste0(result_dir,'/alpha/',filenam))
  markdf=cbind(node=nodes-50,markdf,alpha_hat=round(postmean[nodes-50],2))
  if (nrow(markdf)>0){
    write.csv(markdf,paste0(result_dir,'/node_table/',rownames(pjapmat)[i],'_lambda',lambda,'.csv'))
  }
}
saveRDS(pjapmat,paste0(result_dir,'/PJAP_table/lambda',lambda,'.RData'))








