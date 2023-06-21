# correlations - mainly on TM data
# as less noise 
library(Seurat)
library(Hmisc)
library("BiocNeighbors")
TM <- readRDS('TM_norm.rds')
# correlations 

info <- as.data.frame(table(TM$tissue, TM$cell_ontology_class))
info.sub <- info[info$Freq>100,]
tissue_list <- list()
for (tissue in unique(info.sub$Var1)){
  cell.type <- droplevels(info.sub[info.sub$Var1==tissue,2])
  tissue_list[[tissue]] <- TM[,TM$tissue==tissue & TM$cell_ontology_class %in% cell.type]
  
}

celltype_list_all <- list()
for (tissue in names(tissue_list)){
  for (type in unique(tissue_list[[tissue]]$cell_ontology_class)){
    name <- paste(tissue,type,sep='.' )
    dat1 <- subset(tissue_list[[tissue]],subset=cell_ontology_class==type)
    celltype_list_all[[name]] <- dat1
  }
}

mmu_es <- read.csv("scEssential_mmu.csv", row.names = 1)
celltype_list_es <- lapply(celltype_list_all, function(mat)subset(mat, features=mmu_es$x))
celltype_list_es_cor <- lapply(celltype_list_es,function(mat){
  mat= as.matrix(mat@assays$SCT@data)
  mat = mat[rowSums(mat>0)>10,]
  res2 <- rcorr(t(mat), type = "spearman")
  cor_met <- flattenCorrMatrix(res2$r, res2$P)})

# correlatePairs
library(Seurat)
library(scran)
celltype_list_es_cor1 <- lapply(celltype_list_es,function(mat){
  mats = as.SingleCellExperiment(mat)
  mats = mats[rowSums(counts(mats)>0)>10,]
  correlatePairs(mats)})

# knn method
x = celltype_list_es[[1]]@assays$SCT@count
x = cor(t(as.matrix(x)))
pseudobulk.df <- data.frame()
celltype_list_es_cor2 <- lapply(celltype_list_es,function(mat){
  dat = as.matrix(mat@assays$SCT@data)
  dat = dat[rowSums(dat>0)>10,]
  matrix.knn <- findKNN(t(dat), k=10, BNPARAM=KmknnParam(), get.distance=FALSE)
  raw = as.matrix(mat@assays$SCT@counts)
  for(i in 1:ncol(raw)){
    pseudobulk.exp <- rowSums(raw[,matrix.knn$index[i,]])
    pseudobulk.df <- rbind(pseudobulk.df,pseudobulk.exp)
  }
  pseudobulk.df <- t(pseudobulk.df)
  rownames(pseudobulk.df) <- rownames(raw)
  colnames(pseudobulk.df) <- colnames(raw)
  pseudobulk.seurat <- CreateSeuratObject(pseudobulk.df)
  pseudobulk.seurat <- SCTransform(pseudobulk.seurat,return.only.var.genes = F)
  res  <- rcorr(as.matrix(t(as.matrix(pseudobulk.seurat@assays$SCT@data))))
  cor_met <- flattenCorrMatrix(res$r, res$P)
  
   })
 
boxplot(x[1,],each.exp)

hist(celltype_list_es_cor2[[1]])

# how about random?
set.seed(123441)  
random_sample <- sample(nrow(TM),nrow(mmu_es))
celltype_list_random <- lapply(celltype_list_all, function(mat)subset(mat, features=rownames(TM)[random_sample]))
celltype_list_random_cor <- lapply(celltype_list_random,function(mat){
  mat= as.matrix(mat@assays$SCT@data)
  mat = mat[rowSums(mat>0)>10,]
  res2 <- rcorr(t(mat), type = "spearman")
  cor_met <- flattenCorrMatrix(res2$r, res2$P)})


celltype_list_random_cor1 <- lapply(celltype_list_random,function(mat){
  mats = as.SingleCellExperiment(mat)
  mats = mats[rowSums(counts(mats)>0)>10,]
  correlatePairs(mats)})

pseudobulk2.df <- data.frame()
celltype_list_random_cor2 <- lapply(celltype_list_random,function(mat){
  dat = as.matrix(mat@assays$SCT@data)
  dat = dat[rowSums(dat>0)>10,]
  matrix.knn <- findKNN(t(dat), k=10, BNPARAM=KmknnParam(), get.distance=FALSE)
  raw = as.matrix(mat@assays$SCT@counts)
  for(i in 1:ncol(raw)){
    pseudobulk.exp <- rowSums(raw[,matrix.knn$index[i,]])
    pseudobulk2.df <- rbind(pseudobulk2.df,pseudobulk.exp)
  }
  pseudobulk2.df <- t(pseudobulk2.df)
  rownames(pseudobulk2.df) <- rownames(raw)
  colnames(pseudobulk2.df) <- colnames(raw)
  pseudobulk.seurat <- CreateSeuratObject(pseudobulk2.df)
  pseudobulk.seurat <- SCTransform(pseudobulk.seurat,return.only.var.genes = F)
  res  <- rcorr(as.matrix(t(as.matrix(pseudobulk.seurat@assays$SCT@data))))
  cor_met <- flattenCorrMatrix(res$r, res$P)
  
})


# More significant correlation? #######
 
Num_cor_sig <- sapply(celltype_list_es_cor1, function(mm){
  (sum(mm$FDR<0.05,na.rm=T))
  }) 
Num_cor_sig_random <- sapply(celltype_list_random_cor1, function(mm){
  (sum(mm$FDR<0.05,na.rm=T))
  })
wilcox.test(Num_cor_sig, Num_cor_sig_random)
boxplot(Num_cor_sig, Num_cor_sig_random, names=c("Essential genes",'Random genes'), ylab =  "Number of correlations", main='Significantly correlated gene pairs (p-adj<0.05')

barplot(Num_cor_sig,Num_cor_sig_random)
par(mfrow=c(8,9))
for (i in 1:68){
  #barplot(c(Num_cor_sig[i],Num_cor_sig_random[i]), names.arg = c("Essential genes",'Random genes'))
  print(chisq.test(c(Num_cor_sig[i],Num_cor_sig_random[i])))
}

# more strong correlation?
Num_cor_sig_strong <- sapply(celltype_list_es_cor1, function(mm){
  (sum(abs(mm$rho)>0.3,na.rm=T))
}) 
Num_cor_sig_random_strong <- sapply(celltype_list_random_cor1, function(mm){
  (sum(abs(mm$rho)>0.3,na.rm=T))
})
max(Num_cor_sig_random_strong)/550000
max(Num_cor_sig_strong)/550000
wilcox.test(Num_cor_sig_strong, Num_cor_sig_random_strong)
boxplot(log1p(Num_cor_sig_strong), log1p(Num_cor_sig_random_strong), names=c("Essential genes",'Random genes'), ylab =  "Log(Number of strong  correlations)", main='|Correlation coefficient|>0.3')


# igraph
library(igraph)
cor_sig_strong <- lapply(celltype_list_es_cor1,function(x)x[abs(x$rho)>0.3,])
cor_sig_random_strong <- lapply(celltype_list_random_cor1,function(x)x[abs(x$rho)>0.3,])
cor_sig_sig <- lapply(celltype_list_es_cor1,function(x)x[x$FDR<0.05,])
cor_sig_random_sig <- lapply(celltype_list_random_cor1,function(x)x[x$FDR<0.05,])

gf <- lapply(cor_sig_sig, function(ss)graph_from_edgelist(as.matrix(ss[,1:2])))
gf_random <- lapply(cor_sig_random_sig, function(ss)graph_from_edgelist(as.matrix(ss[,1:2])))
gf_betweeness <-   sapply(gf,betweenness)  
gf_random_betweeness <-   sapply(gf_random,betweenness)  
boxplot(gf_betweeness,gf_random_betweeness)
boxplot( log1p(gf_betweeness) , log1p(gf_random_betweeness),names=c('Essential genes','Random genes'), ylab = 'Log(betweenness)', main= 'Betweenness of cell type coexpression network')
par(mfrow=c(8,9))
for (i in 1:68){
  #barplot(c(Num_cor_sig[i],Num_cor_sig_random[i]), names.arg = c("Essential genes",'Random genes'))
  boxplot(gf_betweeness[[i]],gf_random_betweeness[[i]])
  print(wilcox.test(gf_betweeness[[i]],gf_random_betweeness[[i]])$p.value)
}

 
gf_degree <-   unlist(sapply(gf,function(net)degree(net,mode = "all", normalized = T)))
gf_random_degree <-   unlist(sapply(gf_random,function(net)degree(net,mode = "all", normalized = T)))
boxplot(gf_degree,gf_random_degree,names=c('Essential genes','Random genes'), ylab = 'Normalised node degree', main= 'Normalised node degree across cell types coexpression network')
wilcox.test(gf_degree,gf_random_degree)


gf_eigen_centrality <- unlist(sapply(gf,eigen_centrality))
gf_random_eigen_centrality <- unlist(sapply(gf_random, eigen_centrality))
boxplot(as.numeric(gf_eigen_centrality),as.numeric(gf_random_eigen_centrality))

# Any robust correlation gene pairs??
library(stringr)
cor_sig_gps <- (sapply(cor_sig_sig,function(x){apply(x,1,function(dd)paste(dd[1],dd[2],sep='_'))}))
cor_sig_gps.1 <- unlist(cor_sig_gps)
xx <- table(cor_sig_gps.1)[order(table(cor_sig_gps.1))] 
xx_top30 <- xx[xx>30]
# edger table list for network 
core_edges <- names(xx_top30)
core_edges1 <- sapply(core_edges,function(s)str_split(s,'_')[[1]][1])
core_edges2 <- sapply(core_edges,function(s)str_split(s,'_')[[1]][2])
core_edges_df <- cbind(core_edges1, core_edges2)
colnames(core_edges_df) <- c('gene1','gene2')
# update the core_edges dataframe Feb 2023
index <- apply(core_edges_df,1,function(s)
  (s[1]%in%rownames(mmu_es))&(s[2]%in%rownames(mmu_es)))
core_edges_df1 <- core_edges_df[index,]
library(igraph)
core_gf <- graph_from_edgelist(as.matrix(core_edges_df1))
xx_top30 <- xx_top30[index]
E(core_gf)$width <- as.numeric(xx_top30)/58

net.sp <- delete_vertices(core_gf, c('Jund','Mcl1'))


plot(net.sp, vertex.label.color="black", vertex.label.dist=1, vertex.size=7, edge.arrow.size	=0 )
#cluster_fast_greedy
clp <- cluster_fast_greedy(as.undirected(net.sp))
l <- layout_with_fr(net.sp)
plot(clp, net.sp,vertex.label.color="black", vertex.label.dist=1.5, vertex.size=8, edge.arrow.size	=0,  edge.color="gray70")
module1 <- clp$names[clp$membership==1] #MYC targets
module2 <- clp$names[clp$membership==2]#OXIDATIVE_PHOSPHORYLATION
module3 <- clp$names[clp$membership==3] #CELLULAR_RESPONSE_TO_STRESS
module4 <- clp$names[clp$membership==4]#APOPTOTIC_PROCESS
write.csv(module1,'module1_coexp.csv')
write.csv(module2,'module2_coexp.csv')
write.csv(module3,'module3_coexp.csv')
write.csv(module4,'module4_coexp.csv')
# KEGG pathways 
library(KEGGgraph)
library(KEGGREST)
library(Rgraphviz)
library(stringr)
library(RCy3)
library(msigdbr)

KEGG_names <-  names(table(msigdbr(species = "Mus musculus", category = "C2",subcategory = "CP:KEGG")%>% dplyr::select(gs_name)))

Wnt_DE <- ekk1@result[7,"geneID"]
Wnt_DE_name <- CYCNB_DE_genes[CYCNB_DE_genes$ENTREZID%in% str_split(Wnt_DE,"/")[[1]],1] #65
Wnt_DE_name.fc <- CYCNB_DE.sig[rownames(CYCNB_DE.sig)%in%Wnt_DE_name,1]
names(Wnt_DE_name.fc) <- Wnt_DE_name

CYC_wnt <- pathway_graph(pathway = 'mmu04310.xml',genes = Wnt_DE_name.fc ,org.db = org.Mm.eg.db)
V(LIN_wnt)$name
write.graph(CYC_wnt,"CYC_wnt.gml",format = "gml" )
library(KEGG.db)
pName <- "p53 signaling pathway"
pId <- mget(pName, KEGGPATHNAME2ID)[[1]]
tmp <- tempfile()
retrieveKGML(pId, organism="mmu", destfile=tmp, method="wget", quiet=TRUE)

# functions ####
 
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  = (cormat)[ut],
    p = pmat[ut]
  )
}

makeAttr <- function(graph, default, valNodeList) {
  tmp <- nodes(graph)
  x <- rep(default, length(tmp)); names(x) <- tmp
  
  if(!missing(valNodeList)) {
    stopifnot(is.list(valNodeList))
    allnodes <- unlist(valNodeList)
    stopifnot(all(allnodes %in% tmp))
    for(i in seq(valNodeList)) {
      x[valNodeList[[i]]] <- names(valNodeList)[i]
    }
  }
  return(x)
}


pathway_graph <- function(pathway, genes=NULL, org.db){
  # pathway- pathway.xml file
  # genes- DE genes with FC info
  # org.db- annotation reference org.db
  `%notin%` <- Negate(`%in%`)
  library(KEGGgraph)
  pwayG <- parseKGML2Graph(pathway,expandGenes=TRUE,genesOnly=T)
  outs.all <- sapply(edges(pwayG), length) > 0
  ins.all <- sapply(inEdges(pwayG), length) > 0
  ios.all <- outs.all | ins.all
  
  if(require(org.Mm.eg.db)) {
    ioGeneID.all <- translateKEGGID2GeneID(names(ios.all))
    nodesNames.all <- sapply(mget(ioGeneID.all, org.Mm.egSYMBOL, ifnotfound=NA), "[[",1)
  } else {
    nodesNames.all <- names(ios.all)
  }
  names(nodesNames.all) <- names(ios.all)
  nAttrs.all <- list()
  nAttrs.all$fillcolor <- makeAttr(pwayG, "lightgrey", list(orange=names(ios.all)[ios.all]))
  nAttrs.all$label <- nodesNames.all
  
  library(igraph)
  g <- graph_from_graphnel(pwayG)
  if (is.null(genes)){
    V(g)$name <- nodesNames.all
    return(g)
  }
  else{
    V(g)$name <- nodesNames.all
    labels <-  V(g)$name[which(V(g)$name%notin%names(genes))]
    #print(labels)
    gsub <- delete.vertices(g , labels)
    gsub.fc <- genes
    #node.color <-ifelse(gsub.fc>0,"#a62117","#076e0f") # red- upregulated
    names(node.color) <- names(gsub.fc)
    
    Lout = layout_with_kk(gsub)
    Isolated = which(degree(gsub)==0)
    g.sub = delete.vertices(gsub, Isolated)
    Lout = Lout[-Isolated,]
    node.color.sub <- node.color[-Isolated]
    E(g.sub)$color <- "lightgrey"
    V(g.sub)$label.cex =1 
    plot(g.sub, edge.arrow.size=0.5,vertex.size=12, vertex.label.color= "white", vertex.color=node.color.sub,layout=Lout )
    return(g.sub)
  }
  
  
}

