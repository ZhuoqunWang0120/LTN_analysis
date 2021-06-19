#!/usr/bin/env Rscript
library(phyloseq)
library(ape)
library(data.tree)
argv=commandArgs(TRUE)
WORK_DIR=argv[1]
if (!file.exists(paste0(WORK_DIR,"/cache/ps_otu50.RData"))){
  raw_data = read.csv(paste0(WORK_DIR,"/res/diabimmune_t1d_16s_otu_table.txt"), sep = "\t", row.names = 1 )
  otutab0=apply(as.matrix(raw_data[,-ncol(raw_data)]),1:2,as.numeric)
  rownames(otutab0)=rownames(raw_data)
  colnames(otutab0)=colnames(raw_data)[-ncol(raw_data)]
  taxtab0=do.call(rbind,sapply(as.character(raw_data[,ncol(raw_data)]),function(x){strsplit(x,';  ')}))
  rownames(taxtab0)=rownames(otutab0)
  colnames(taxtab0)=c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
  ps0=phyloseq(otu_table(otutab0,taxa_are_rows = T),tax_table(taxtab0))
  ra=apply(otutab0,2,function(x){x/sum(x)})
  otu100=names(sort(rowSums(ra),decreasing = T)[1:100])
  ps100=prune_taxa(otu100,ps0)
  tax100=tax_table(ps100)
  # build tree
  taxpath=matrix(rep('ROOT',nrow(tax100)),ncol=1)
  tax100_otu=cbind(tax100,rownames(tax100))
  taxranks=c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
  taxranks=c(taxranks,'OTUlevel')
  colnames(tax100_otu)=taxranks
  set.seed(1)
  for (i in 1:8){ # i=1: ROOT, i=2: kingdom
    taxmat=cbind(taxpath,tax100_otu[,i])
    nodes=unique(taxpath)
    taxpath=apply(taxmat,1,function(x){paste(x,collapse = '/')})
    nchild=sapply(nodes, function(x){length(unique(taxmat[taxmat[,1]==x,2]))})
    children=lapply(nodes, function(x){unique(taxmat[taxmat[,1]==x,2])})
    for (j in 1:length(nchild)){
      if (nchild[j]==1){
        taxpath[taxmat[,1]==nodes[j]]=nodes[j]
      } else if (nchild[j]>2){
        randphylo=rtree(n=nchild[j],tip.label = children[[j]],rooted = T) 
        randNode=as.Node(randphylo)
        randpath=nodepath(randphylo) # pre-order, [ NOT the order of children[[j]] !!! ]
        leaf_preorder=randNode$Get('name',traversal = 'pre-order',filterFun = isLeaf)
        leaf_preorder=gsub(' ','_',leaf_preorder)
        randpath_named=lapply(1:nchild[j], function(x){paste0(paste(paste0(taxranks[i],randpath[[x]][-c(1,(length(randpath[[x]])))]),collapse = '/'),'/',leaf_preorder[x])})
        randpath_named=lapply(randpath_named,function(x){gsub(paste0(taxranks[i],'/'),'',x)})
        for (k in 1:nchild[j]){
          taxpath[taxpath==paste0(nodes[j],'/',leaf_preorder[k])]=paste0(taxmat[taxpath==paste0(nodes[j],'/',leaf_preorder[k]),1],'/',randpath_named[[k]])
        }
      }
    }
  }
  pathdf=data.frame(pathString=taxpath)
  pathlist=apply(pathdf,1,function(x){strsplit(x,split='/')})
  pathdf_otu=data.frame(pathString=do.call(c,lapply(1:100,function(x){paste(c(pathlist[[x]][[1]][-length(pathlist[[x]][[1]])],rownames(pathdf)[x]),collapse = '/')})))
  tree=as.Node(pathdf_otu)
  phylotree=as.phylo(tree)
  # pre-order traversal
  otu_preorder=tree$Get('name',traversal = 'pre-order',filterFun = isLeaf)
  otutab100=otu_table(ps100)[otu_preorder,]
  taxtab100=tax100[otu_preorder,]
  ps100=phyloseq(otu_table(t(otutab100),taxa_are_rows = F),tax_table(taxtab100),phylotree)
  load(paste0(WORK_DIR,"/res/diabimmune_t1d_16s_metadata.rdata")) # md_16S
  rownames(md_16S)=as.character(md_16S$G_id)
  md_16S=md_16S[rownames(otu_table(ps100)),]
  ps100=merge_phyloseq(ps100,sample_data(md_16S))
  saveRDS(ps100,paste0(WORK_DIR,"/cache/ps_otu100.RData"))
  # ra=apply(otutab100,2,function(x){x/sum(x)})
  otu50=names(sort(rowSums(ra),decreasing = T)[1:50])
  ps50=prune_taxa(otu50,ps100)
  saveRDS(ps50,paste0(WORK_DIR,"/cache/ps_otu50.RData"))
}
