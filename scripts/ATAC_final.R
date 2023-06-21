# Find accessibility values for essential genes
library(monocle3)
library(cicero)
library(data.table)
library(Seurat)
`%notin%` <- Negate(`%in%`)

metadata = fread('cell_metadata.txt')
metadata = metadata[match(colnames(atac_matrix.binary),metadata$cell),]
row.names(metadata)= metadata$cell
rowdata = fread('peak_promoter_intersections.txt')
input_cds <- new_cell_data_set(atac_matrix.binary ,  cell_metadata = metadata)
input_cds1 <- input_cds[match(rowdata$peak_id,rownames(input_cds)),] 
rowData(input_cds1) <- rowdata
Marrow_cds <- input_cds1[,pData(input_cds1)$tissue=="BoneMarrow"]
Marrow_list <- list()
for(i in unique(pData(Marrow_cds)$cell_label)){
  if (sum(pData(Marrow_cds)$cell_label==i)>100){
    Marrow_list[i]= Marrow_cds[,pData(Marrow_cds)$cell_label==i]
  }
}
scEssential_mmu  <- readRDS("D:/PhD/Essential genes/Revisions for thesis/revised figures/scEssential_mmu.rds")
Marrow_list <- Marrow_list[-2] # remove unknown
Marrow_es <- sapply(Marrow_list, function(x)exprs(x)[(rowData(x)$gene_short_name%in%rownames(scEssential_mmu)),])
Marrow_es_gene <- rowdata[rowdata$peak_id%in%rownames(Marrow_es[[1]]),]

Marrow_es_no <- sapply(Marrow_list, function(x){
  set.seed(123441)
  x = exprs(x)[(rowData(x)$gene_short_name%notin%rownames(scEssential_mmu)),]
  x[sample(nrow(x),nrow(Marrow_es[[1]])),]
  })

boxplot(sapply(Marrow_es_no,rowMeans))
boxplot(sapply(Marrow_es_no,function(y)rowMeans(y>0)))

# compared with random genes
ES_atac = sapply(Marrow_es,rowMeans)
ES_no_atac = sapply(Marrow_es_no,rowMeans)
par(mfrow=c(3,3))
for (plot in 1:ncol(ES_atac)){
  boxplot(ES_atac[,plot],ES_no_atac[,plot],notch = T,names=c('Essential genes',"random genes"), ylab = "Accessibility from scATACseq data",main=paste(colnames(ES_no_atac)[plot]))
  print(wilcox.test(ES_atac[,plot],ES_no_atac[,plot])$p.value )
}

# gene level, average the peaks for one gene.
Marrow_gene_list <- sapply(Marrow_list, function(ct)aggregate(exprs(ct), list(rowData(ct)$gene_short_name), mean))
Marrow_es_gene <- sapply(Marrow_gene_list, function(x) x[x$Group.1%in% rownames(scEssential_mmu),])
Marrow_es_gene_no <- sapply(Marrow_gene_list, function(x){
  set.seed(123441)
  x = x[x$Group.1%notin%rownames(scEssential_mmu),]
  x[sample(nrow(x),nrow(Marrow_es_gene[[1]])),]
})
# compared with random genes
ES_gene_atac = sapply(Marrow_es_gene,function(y)rowMeans(y[,-1]))
ES_gene_no_atac = sapply(Marrow_es_gene_no,function(y)rowMeans(y[,-1]))
par(mfrow=c(3,3))
for (plot in 1:ncol(ES_gene_atac)){
  boxplot(ES_gene_atac[,plot],ES_gene_no_atac[,plot],notch = T,names=c('Essential genes',"random genes"), ylab = "Accessibility from scATACseq data",main=paste(colnames(ES_gene_no_atac)[plot]), col=c("#c37ba6", "#f4bd8a"))
  print(wilcox.test(ES_gene_atac[,plot],ES_gene_no_atac[,plot])$p.value )
}
 
# Does accessibility correlate with essentiality score?
rownames(ES_gene_atac) <- Marrow_es_gene[["Hematopoietic progenitors"]][["Group.1"]]
scEssential_mmu$genename <- rownames(scEssential_mmu)
scEssential_mmu <- scEssential_mmu[order(scEssential_mmu$ES_score),]
ES_score <- scEssential_mmu[scEssential_mmu $genename%in%rownames(ES_gene_atac),]
par(mfrow=c(2,3))
for (gx in c(1,3,4,7,8)){
  plot(ES_score$ES_score, ES_gene_atac[,gx],xlab ="Essentiality score", ylab = "Accessibility from scATACseq data",main=paste(colnames(ES_gene_no_atac)[gx]))
  #print(cor.test(ES_score$ES_score, ES_gene_atac[,gx])$p.value )
  print(cor(ES_score$ES_score, ES_gene_atac[,gx], method = "spearman") )
}
plot(rowMeans(ES_gene_atac) ,ES_score$ES_score,xlab ="Essentiality score", ylab = "Accessibility from scATACseq data", main="Average across cell types")
cor(rowMeans(ES_gene_atac) ,ES_score$ES_score, method = 'spearman') 

# Select five cell types that overlap with TMS
# accessibility for cell types overlap with TMS
ES_atac.5 <- ES_atac[,c(1,3,4,7,8)]
ES_atac.5 <- reshape2::melt(ES_atac.5)
ES_atac_no.5 <- ES_no_atac[,c(1,3,4,7,8)]
ES_atac_no.5 <- reshape2::melt(ES_atac_no.5)
ES_gene_atac.5 <- ES_gene_atac[,c(1,3,4,7,8)]
ES_gene_atac.5 <- reshape2::melt(ES_gene_atac.5)
ES_gene_no_atac.5 <- ES_gene_no_atac[,c(1,3,4,7,8)]
ES_gene_no_atac.5 <- reshape2::melt(ES_gene_no_atac.5)

ES_all <- rbind(ES_atac.5[,-1],ES_atac_no.5[,-1],ES_gene_atac.5[,-1],ES_gene_no_atac.5[,-1])
ES_all$type <- c(rep("scEssentials",nrow(ES_atac.5[,-1])),rep("Random genes",nrow(ES_atac_no.5[,-1])),rep("scEssentials",nrow(ES_gene_atac.5[,-1])),rep("Random genes",nrow(ES_gene_no_atac.5[,-1])))

ES_all$category <- c(rep("Peak", 2*nrow(ES_atac.5[,-1])), rep('Gene', 2*nrow(ES_gene_atac.5[,-1]))) 
ES_all$type <- factor(ES_all$type ,     
                         c("scEssentials", "Random genes"))
library(ggplot2)
library(ggpubr)
ggplot(ES_all, aes(x=Var2, y=value )) + 
  geom_boxplot(aes(fill=type),outlier.size = 0.3)+ facet_wrap(~ category)+xlab('')+ylab('')+stat_compare_means(aes(group = type),label = "p.signif")+theme_bw(base_size = 12)+
  scale_fill_manual(values=c("#c37ba6", "#f4bd8a")) +
  theme(legend.title=element_blank())

# library(dplyr)
# xx <- as_tibble(ES_gene_atac.5) %>%
#   group_by(Gene) %>% dplyr::summarise(mean = mean(Accessibility))
ES_gene_atac.x <- ES_gene_atac.5[ES_score$genename%in%ES_gene_atac.5$Var1, ]
ESs <- ES_score[match(ES_gene_atac.x$Var1,ES_score$genename),]
ES_all_2 <- cbind(ES_gene_atac.x,ESs$ES_score)
colnames(ES_all_2) <- c('Gene','Celltype','Accessibility','ES')
ggplot(ES_all_2, aes(x=log(ES), y=Accessibility))+ geom_point(alpha=0.6)+ stat_smooth(method="lm" )+ theme_bw() +theme(legend.position="none")+ylab("Accessibility level from scATAC-seq") + xlab('Essentiality score')+ facet_wrap(~ Celltype)+
  stat_cor(method='spearman',cor.coef.name="rho",label.y = 0.5)

# Are there more TFs?
Mouse_TF <- fread('Mus_musculus_TF.txt')
Mouse_TF_ES <- Mouse_TF[Mouse_TF$Symbol%in%rownames(scEssential_mmu),]
phyper(nrow(Mouse_TF_ES)-1, length(unique(Mouse_TF$Symbol)), 24351-length(unique(Mouse_TF$Symbol)), 733,lower.tail= FALSE) 
TF_family <- as.data.frame(table(Mouse_TF$Family))
TF_ES_family <- as.data.frame(table(Mouse_TF_ES$Family))
 
TF_family_all <- merge(TF_family,TF_ES_family, by="Var1")
barplot(TF_family_all[,3]/TF_family_all[,2])

Marrow_assay <- CreateChromatinAssay(counts = counts(Marrow_cds), min.cells = 3, min.features = 200,sep = c("_", "_") )
Marrow_seu <-  CreateSeuratObject(
  counts = Marrow_assay,
  assay = "peaks",
  meta.data = as.data.frame(pData(Marrow_cds))
)
library(EnsDb.Mmusculus.v79)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
Marrow_seu <- RunTFIDF(Marrow_seu)
Marrow_seu <- FindTopFeatures(Marrow_seu, min.cutoff = 'q0')
Marrow_seu <- RunSVD(object = Marrow_seu)
DepthCor(Marrow_seu) # remove first two component >0.5
Marrow_seu <- RunUMAP(
  object = Marrow_seu,
  reduction = 'lsi',
  dims = 3:30
)
DimPlot(object = Marrow_seu, group.by = "cell_label")
VlnPlot(Marrow_seu,features = "chr7-29704333-29705255",group.by = 'cell_label')
DotPlot(Marrow_seu,features = c("chr11-93855923-93858996","chr7-29704333-29705255"),group.by = 'cell_label')
FeaturePlot(
  object = Marrow_seu,
  features = "chr11-93855923-93858996",
  pt.size = 0.1 
)
boxplot( Marrow_seu@assays$peaks@data[3,])
Marrow_cds_norm <- new_cell_data_set(Marrow_seu@assays$peaks@data,cell_metadata =  Marrow_seu@meta.data)
 