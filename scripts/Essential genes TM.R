# Essential genes in Tabula Muris 
library("TabulaMurisData")
library(ExperimentHub)
library(SingleCellExperiment)
library(Seurat)
library(UpSetR)
library(clusterProfiler)
library(org.Mm.eg.db)
library(data.table)
library(umap)
library(corrplot)
library(ggfortify)
library(corrplot)
library(RColorBrewer)
library(ggplot2)
library(ggridges)
library(stringi)
library(dplyr)
library(msigdbr)
library(biomaRt)
memory.limit(230*1024*1024)
eh <- ExperimentHub()
query(eh, "TabulaMurisData")
facs <- eh[["EH1618"]]
 
info <- as.data.frame(table(facs$tissue, facs$cell_ontology_class))

# Does essential genes express in all the tissues and have a relatively higher expression? #######
# subset the data that only cell type has more than 100 cells are included - 68 tissue-cell type
info.sub <- info[info$Freq>100,]
tissue_list <- list()
for (tissue in unique(info.sub$Var1)){
  cell.type <- droplevels(info.sub[info.sub$Var1==tissue,2])
  tissue_list[[tissue]] <- facs[,facs$tissue==tissue & facs$cell_ontology_class %in% cell.type]
  
}

# tissue level #
# tissue_norm_list <- lapply(tissue_list,function(x){
#   scran::computeSumFactors(x)
#   scater::logNormCounts(x)
# })
# keep consistent with previous - use scTransform 

tissue_normed_list <- lapply(tissue_list, function(met){
  met <- met[!grepl("^ERCC-", rownames(met)),]
  met <- as.Seurat(met, data = NULL)
  RenameAssays(met,originalexp = 'RNA' )
  met <- SCTransform(met, assay = 'originalexp',vars.to.regress = "mouse_id",verbose=TRUE)
})
par(mfrow=c(3,6))
lapply(tissue_normed_list,function(x)boxplot(rowMeans(as.matrix(x@assays$SCT@data)), main=unique(x$tissue)))

celltype_normed_list <- list()
for (tissue in names(tissue_normed_list)){
  for (type in unique(tissue_normed_list[[tissue]]$cell_ontology_class)){
    name <- paste(tissue,type,sep='.' )
    dat <- subset(tissue_normed_list[[tissue]],subset=cell_ontology_class==type)
    celltype_normed_list[[name]] <- dat}
}
 
par(mfrow=c(8,9))
lapply(celltype_normed_list,function(x)boxplot(rowMeans(as.matrix(x@assays$SCT@data)), main=paste(unique(x$tissue),unique(x$cell_ontology_class),sep='.')))

# extract for ES genes 
ES_genes_mmu <- read.csv( "essential_genes_mmu_final.csv")
colnames(ES_genes_mmu) <- c("Name","Rank")
ES_genes_tbl <- bitr(ES_genes_mmu$Name,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Mm.eg.db)

tissue_ES_list <- lapply(tissue_normed_list,function(data){
  data <- as.matrix(data@assays$SCT@data)
  data <- data[rownames(data)%in%ES_genes_tbl$SYMBOL,]
})

barplot(sapply(tissue_ES_list,function(x) sum(rowSums(x)>0)))
# how many genes shared in the end?
ES_genes_mmu_final_list <- lapply(tissue_ES_list, rownames)
upset(fromList(ES_genes_mmu_final_list), order.by = "freq",nsets = 18)
ES_genes_mmu_final <- Reduce(intersect,ES_genes_mmu_final_list)


celltype_ES_list <- lapply(celltype_normed_list,function(data){
  data <- as.matrix(data@assays$SCT@data)
  index <-  apply(data,1,function(y)sum(y>0))
  data <- data[index>3,]
  data <- data[rownames(data)%in%ES_genes_mmu_final,]
})
# ES genes are detectable
barplot(sapply(celltype_ES_list,function(x) sum(rowSums(x)>0))[order(sapply(celltype_ES_list,function(x) sum(rowSums(x)>0)),decreasing = T)], names.arg  = '',main = 'Number of essential genes across 68 cell types (TM)')
median(sapply(celltype_ES_list,nrow))
ES_genes_mmu_final_update <- Reduce(intersect,lapply(celltype_ES_list,rownames)) #775

# update the tissue ES list, remove some variable genes
tissue_ES1_list <- lapply(tissue_ES_list,function(x)x[rownames(x)%in%ES_genes_mmu_final_update,])
# update so have the same dimension 
celltype_ES1_list <- lapply(celltype_ES_list,function(x)x[rownames(x)%in%ES_genes_mmu_final_update,])
celltype_ES1_tbl <-  do.call(cbind,lapply(celltype_ES1_list,rowMeans))
hist(rowSds(celltype_ES1_tbl),xlab = '', main='Standard deviation of the distribution of the ranks across cell types (TM)') 
ES_genes_mmu_final_update1 <- ES_genes_mmu_final_update[rowSds(celltype_ES1_tbl)< 4*sd(rowSds(celltype_ES1_tbl))]
write.csv(ES_genes_mmu_final_update1,"scEssential_mmu.csv")
saveRDS(ES_genes_mmu_final_update1,'scEssential_mmu.rds')

# update the list for tissues and cell types 
# update the tissue ES list, remove some variable genes
tissue_ES1_list <- lapply(tissue_ES_list,function(x)x[rownames(x)%in%ES_genes_mmu_final_update1,])
# update so have the same dimension 
celltype_ES1_list <- lapply(celltype_ES_list,function(x)x[rownames(x)%in%ES_genes_mmu_final_update1,])

# overall ES genes have higher exp
par(mfrow=c(3,6))
sapply(tissue_normed_list, plot_mean_exp)

par(mfrow=c(5,6))
sapply(celltype_normed_list[1:30], plot_mean_exp)#1600*1000
par(mfrow=c(5,6))
sapply(celltype_normed_list[31:60], plot_mean_exp)
sapply(celltype_normed_list[61:68], plot_mean_exp)
# par(mfrow=c(8,9))
# sapply(celltype_normed_list, plot_mean_exp)

tissue_ES_rank_list <- lapply(tissue_normed_list, get_ES_rank)
par(mfrow=c(3,6))
sapply(tissue_ES_rank_list,boxplot)

celltype_ES_rank_list <- lapply(celltype_normed_list, get_ES_rank)
par(mfrow=c(3,6))
lapply(tissue_ES1_list,function(x)boxplot(rowMeans(x)))

# heat map and correlation ##########
##just with the expression level ####
tissue_ES1_tbl <-  do.call(cbind,lapply(tissue_ES1_list,rowMeans))

tissue_mat <- cor(tissue_ES1_tbl)
# distribution of the correlations between tissue
hist(tissue_mat[upper.tri(tissue_mat)])
corrplot(tissue_mat, method = 'circle',tl.pos = 'l',addCoef.col = 'black', tl.col = 'black',order = 'hclust',cl.pos = 'b',cl.ratio = 0.1, tl.cex = 0.6, col = colorRampPalette(c("blue","white","yellow","red"))(100))
 
heatmap(t(tissue_ES1_tbl),scale ='row')
# distribution of the correlations between cell types
celltype_ES1_tbl <-  do.call(cbind,lapply(celltype_ES1_list,rowMeans))

celltype_mat <- cor(celltype_ES1_tbl)
# distribution of the correlations between cell types
hist(celltype_mat[upper.tri(celltype_mat)])
colnames(celltype_mat) <- rep(" ", NROW(celltype_mat))
corrplot(celltype_mat, method = 'circle',tl.pos = 'l',tl.col = 'black',order = 'hclust',cl.pos = 'b',cl.ratio = 0.1, tl.cex = 0.6, col = colorRampPalette(c("blue","white","yellow","red"))(100)) #1600
heatmap(t(celltype_ES1_tbl),scale ='row')
median(celltype_mat)
 
## With rank level #########
tissue_ES_rank_tbl <-  do.call(cbind,tissue_ES_rank_list)
tissue_rank_mat <- cor(tissue_ES_rank_tbl)

# distribution of the correlations between tissue
hist(tissue_rank_mat[upper.tri(tissue_rank_mat)])
corrplot(tissue_rank_mat, method = 'circle',tl.pos = 'l',addCoef.col = 'black', tl.col = 'black',order = 'hclust',cl.pos = 'b',cl.ratio = 0.1, tl.cex = 0.8, col = colorRampPalette(c("blue","white","yellow","red"))(100))
heatmap(t(tissue_ES_rank_tbl),scale ='row')
# distribution of the correlations between cell types
celltype_ES_rank_tbl <-  do.call(cbind,celltype_ES_rank_list)
 
celltype_rank_mat <- cor(celltype_ES_rank_tbl)
# distribution of the correlations between cell types
hist(celltype_rank_mat[upper.tri(celltype_rank_mat)])
corrplot(celltype_rank_mat, method = 'circle',tl.pos = 'l',tl.col = 'black',order = 'hclust',cl.pos = 'b',cl.ratio = 0.1, tl.cex = 0.6, col = colorRampPalette(c("blue","white","yellow","red"))(100))
 

## PCA plots #####
library(ggfortify)
autoplot(prcomp(t(tissue_ES1_tbl)),label = TRUE, label.size = 3)
autoplot(prcomp(t(celltype_ES1_tbl)),label = TRUE, label.size = 3)

# tissue_ES_exp <- do.call(cbind, tissue_ES1_list)
# tissue_umap <- umap(t(tissue_ES_exp)) 
# plot(tissue_umap$layout   )

## ridgeplots ####
# To show the expression distribution variations for essential genes
# mean based
tissue_ES1_tbl1 <- melt(tissue_ES1_tbl)
ggplot(tissue_ES1_tbl1, aes(x = value, y = Var2,fill=as.factor(Var2))) + geom_density_ridges()+theme_bw()
ggplot(tissue_ES1_tbl1, aes(x = value, y = Var2,fill=as.factor(Var2))) + geom_boxplot()+theme_bw()

celltype_ES1_tbl1 <- melt(celltype_ES1_tbl) 
celltype_ES1_tbl1$tissue <- sapply(celltype_ES1_tbl1$Var2,function(cx)stringr::str_split(cx,"\\.")[[1]][1])
celltype_ES1_tbl1$tissue <- factor(celltype_ES1_tbl1$tissue , levels=levels(as.factor(tissue_ES1_tbl1$Var2)))
 
ggplot(celltype_ES1_tbl1, aes(x = value, y = Var2,fill=tissue )) + geom_density_ridges( )+theme_bw()

ggplot(celltype_ES1_tbl1, aes(x = value, y = Var2,fill=as.factor(tissue))) + geom_boxplot()+theme_bw()
 
barplot(sapply(tissue_normed_list, function(x)max(rowMeans(x))), horiz = T)
barplot(sapply(celltype_normed_list, function(x)  max(rowMeans(x))),horiz=T)

# sd based 
tissue_ES_sd_tbl <-  do.call(cbind,lapply(tissue_ES1_list,rowSds))
tissue_ES_sd_tbl1 <- melt(tissue_ES_sd_tbl)
ggplot(tissue_ES_sd_tbl1, aes(x = value, y = Var2,fill=as.factor(Var2))) + geom_density_ridges()+theme_bw()

celltype_ES_sd_tbl <-  do.call(cbind,lapply(celltype_ES_list,rowSds))
celltype_ES_sd_tbl1 <- melt(celltype_ES_sd_tbl) 
celltype_ES_sd_tbl1$tissue <- sapply(celltype_ES_sd_tbl1$Var2,function(cx)str_split(cx,"\\.")[[1]][1])
celltype_ES_sd_tbl1$tissue <- factor(celltype_ES_sd_tbl1$tissue , levels=levels(as.factor(tissue_ES1_tbl1$Var2)))

ggplot(celltype_ES_sd_tbl1, aes(x = value, y = Var2,fill=tissue )) + geom_density_ridges( )+theme_bw()

barplot(sapply(tissue_normed_list, function(x)max(rowSds(as.matrix(x@assays$SCT@data)))), horiz = T)
barplot(sapply(celltype_normed_list, function(x)max(rowSds(as.matrix(x@assays$SCT@data)))),horiz=T)

# Percentage of expression  #########
tissue_ES_pct <- sapply(tissue_ES1_list,get_ES_pct)
tissue_ES_pct.1 <- melt(tissue_ES_pct)
gene_pct_order <- rowMeans(tissue_ES_pct)[order(rowMeans(tissue_ES_pct),decreasing = T)]
ggplot(tissue_ES_pct.1, aes(x = value, y = Var2,fill=as.factor(Var2))) + geom_density_ridges()+theme_bw()
ggplot(tissue_ES_pct.1[tissue_ES_pct.1$Var1%in%names(gene_pct_order)[1:50],], aes(x = value, y = Var1)) + geom_boxplot()+theme_bw()


celltype_ES_pct <- sapply(celltype_ES_list,get_ES_pct)
celltype_ES_pct.1 <- melt(celltype_ES_pct)
gene_pct_order1 <- rowMeans(celltype_ES_pct)[order(rowMeans(celltype_ES_pct),decreasing = T)]

ggplot(celltype_ES_pct.1, aes(x = value, y = Var2,fill=as.factor(Var2))) + geom_density_ridges()+theme_bw()
ggplot(celltype_ES_pct.1[celltype_ES_pct.1$Var1%in%names(gene_pct_order1)[1:50],], aes(x = value, y = Var1)) + geom_boxplot()+theme_bw()


# Packages
library(hexbin)
library(RColorBrewer)
# Make the plot
bin<-hexbin(tissue_ES1_tbl1$value, tissue_ES_pct.1$value, xbins=100)
my_colors=colorRampPalette(rev(brewer.pal(11,'Spectral')))
plot(bin, main="" , colramp=my_colors , legend=F, xlab='ES expression (tissue)', ylab='Percentage of ES expression (tissue)') 

bin <- hexbin(celltype_ES1_tbl1$value, celltype_ES_pct.1$value, xbins=100)
my_colors=colorRampPalette(rev(brewer.pal(11,'Spectral')))
plot(bin, main="" , colramp=my_colors , legend=F, xlab='ES expression (cell type)', ylab='Percentage of ES expression (cell type)') 

# more pct than random?
par(mfrow=c(3,6))
sapply(tissue_normed_list, plot_pct)

par(mfrow=c(5,6))
sapply(celltype_normed_list[1:30], plot_pct)#1600*1000
par(mfrow=c(5,6))
sapply(celltype_normed_list[31:60], plot_pct)
sapply(celltype_normed_list[61:68], plot_pct)
# Marker genes ########
tissue_markers <- lapply(tissue_normed_list,function(x){
  Idents(x) <- x$cell_ontology_class
  markers <- FindAllMarkers(x, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use="MAST")
})
tissue_markers.733 <- lapply(tissue_markers[1:15],function(x){x%>%
  group_by(cluster) %>%
  slice_max(n = 733, order_by = avg_log2FC)%>% 
    dplyr::select(gene) })

tissue_markers.733[["Pancreas"]] <- tissue_markers[[18]]%>%
  group_by(cluster) %>%
  slice_max(n = 733, order_by = avg_log2FC)%>% 
  dplyr::select(gene)

celltype_markers_list <- list()
for (tissue in names(tissue_markers.733)){
  for (type in unique(tissue_markers.733[[tissue]]$cluster)){
    name <- paste(tissue,type,sep='.' )
    dat <-  tissue_markers.733[[tissue]][tissue_markers.733[[tissue]]$cluster==type,'gene']
    celltype_markers_list[[name]] <- dat}
}
upset(fromList(ES_markers_overlap), nsets = 10)
ES_markers_overlap <- sapply(celltype_markers_list,function(x)intersect(x[["gene"]],ES_genes_mmu_final_update1))
barplot(sapply(ES_markers_overlap,length), horiz = T)
boxplot(sapply(ES_markers_overlap,length))
library(stringr)
ES_markers_overlap.1 <- as.data.frame(sapply(ES_markers_overlap,length))
colnames(ES_markers_overlap.1) <- 'value'
ES_markers_overlap.1$celltype <- rownames(ES_markers_overlap.1)
ES_markers_overlap.1$type <- sapply(ES_markers_overlap.1$celltype,function(cx)str_split(cx,"\\.")[[1]][1])

ggplot(ES_markers_overlap.1, aes(x = celltype, y = value ,fill=type )) + geom_bar(stat="identity")+theme_bw()+ coord_flip()+ggtitle("Number of scEssentials overlapped with cell-type markers")+labs(fill='Source tissue')+xlab('')+ylab('')

median(sapply(ES_markers_overlap,length))
# when compared to large marker database
library(data.table)
Panglao <- as.data.frame(fread("PanglaoDB_markers_27_Mar_2020.tsv"))
mouse_Panglao <- Panglao[Panglao$species=='Mm'|Panglao$species=='Mm Hs',]
mouse_Panglao1 <- mouse_Panglao 
mouse_Panglao1$`official gene symbol` <- stringr::str_to_title(mouse_Panglao1$`official gene symbol`)

mouse_Panglao.overlap <- mouse_Panglao1[mouse_Panglao1$`official gene symbol`%in%ES_genes_mmu_final_update1,]
unique(mouse_Panglao.overlap$`cell type`)

mouse_markers_Panglao <- unique(mouse_Panglao$`official gene symbol`)
library(stringr)
mouse_markers_Panglao <- str_to_title(mouse_markers_Panglao)
library(ggvenn)
x = intersect(mouse_markers_Panglao,ES_genes_mmu_final_update1)
mouse_sensitivity <- mouse_Panglao[str_to_title(mouse_Panglao$`official gene symbol`)%in%x,]
mouse_sensitivity_score <- aggregate(mouse_sensitivity$sensitivity_mouse,by=list(mouse_sensitivity$`official gene symbol`) , max)
mouse_sensitivity_score$Group.1 <- str_to_title(mouse_sensitivity_score$Group.1)
colnames(mouse_sensitivity_score) = c('gene','sensitivity')
Number_of_panglao_ES <- sapply(celltype_list_final,function(o)length(intersect(rownames(o),human_markers_Panglao)))
median(Number_of_panglao_ES)/48

# compute the essentially score
xx = as.data.frame(ES_genes_mmu_final_update1)
colnames(xx) = 'gene'
S = dplyr::full_join(xx,mouse_sensitivity_score,by='gene')
S[is.na(S)] <- 0 

num_marker =  unlist(celltype_markers_list)
num_marker <- table(num_marker)
Gx <- as.data.frame(num_marker)
colnames(Gx) <- c("gene","freq")
Gx <- Gx[Gx$gene%in%xx$gene,]
G = dplyr::full_join(xx,Gx, by='gene')
G[is.na(G)] <- 0 
G$freq <- G$freq/63

w = as.data.frame((S$sensitivity + G$freq)/2)
w$gene = S$gene

pct = rowMeans(unlist(celltype_ES_pct))*100
# update pct
pct <- pct[names(pct)%in%ES_genes_mmu_final_update1]
ES_score = (1-w$`(S$sensitivity + G$freq)/2`)*pct
ES_score <- as.data.frame(ES_score)
#write.csv(ES_score,'ES_score_mmu_final.csv')
write.csv(ES_score,'scEssential_mmu.csv')
saveRDS(ES_score,'scEssential_mmu.rds')
# HVG genes ########
celltype_HVGs <- list()
for (tissue in names(tissue_normed_list)){
  for (type in unique(tissue_normed_list[[tissue]]$cell_ontology_class)){
    name <- paste(tissue,type,sep='.' )
    dat <- subset(tissue_normed_list[[tissue]],subset=cell_ontology_class==type)
    dat <- FindVariableFeatures(dat,assay = "originalexp", selection.method = "vst", nfeatures = 733)
    hvgs <-  VariableFeatures(dat,assay = "originalexp",selection.method = "vst")
    celltype_HVGs[[name]] <- hvgs}
}
 
ES_HVGs <- sapply(celltype_HVGs,function(s) intersect(s,ES_genes_mmu_final_update1 ))
barplot(unique(sapply(ES_HVGs,length)), horiz = T)

ES_HVGs.1 <- as.data.frame(sapply(ES_HVGs, length))
colnames(ES_HVGs.1) <- 'value'
ES_HVGs.1$celltype <- rownames(ES_HVGs.1)
ES_HVGs.1$type <- sapply(ES_HVGs.1$celltype,function(cx)str_split(cx,"\\.")[[1]][1])

ggplot(ES_HVGs.1, aes(x = celltype, y = value ,fill=type )) + geom_bar(stat="identity")+theme_bw()+ coord_flip()+ggtitle("Number of scEssentials with HVGs")+labs(fill='Source tissue')
median(sapply(ES_HVGs,length))

###markers+HVGs######
ES_all_overlap <- merge(ES_markers_overlap.1,ES_HVGs.1, by="celltype")
#saveRDS(ES_all_overlap,"TM_overlap_markerHVG.rds")

library(ggalt)
p1 <- ggplot(ES_all_overlap, aes(y=celltype, x=value.x/733, xend=value.y/733)) +
  geom_dumbbell(size=3, color="#e3e2e1",
                colour_x = "#5b8124", colour_xend = "#bad744",
                dot_guide=TRUE, dot_guide_size=0.25) +
  labs(x=NULL, y=NULL, title="Overlap of scEssentials and markers/HVGs") +theme_minimal()

library(scales)
ttt = as.data.frame(table(ES_all_overlap$type.y))
x_cols <- rep(hue_pal()(length(ttt$Var1)),ttt$Freq) p1+
 theme(axis.text.y = element_text(colour=x_cols),
        axis.ticks.y=element_blank())
saveRDS(p1,'overlap_mmu_ggplot.rds')

# correlations #####
tissue_ES_cor <- lapply(tissue_normed_list, get_ES_cor)
tissue_random_cor <-  lapply(tissue_normed_list, get_cor_random)
# Select only significant correlations 
tissue_ES_cor_sig <- lapply(tissue_ES_cor, function(x){x[x$p<0.05,]})
tissue_random_cor_sig<- lapply(tissue_random_cor, function(x){x[x$p<0.05,]})
# are there more signfincant terms than random?
par(mfrow=c(3,6))

sapply(tissue_ES_cor_sig,function(x)boxplot(x$zcor))

apply(cbind(sapply(tissue_ES_cor_sig,function(x) dim(x)[1]),sapply(tissue_random_cor_sig, function(y)dim(y)[1])),1,function(value)barplot(value, names.arg=c('Essential genes',"random genes"), main=names(value)))
apply(cbind(sapply(tissue_ES_cor_sig,function(x) dim(x)[1]),sapply(tissue_random_cor_sig, function(y)dim(y)[1])),1, chisq.test)

for (i in 1:length(tissue_ES_cor_sig)){
  # boxplot(tissue_ES_cor_sig[[i]]$zcor,tissue_random_cor_sig[[i]]$zcor)
  print(t.test(tissue_ES_cor_sig[[i]]$zcor,tissue_random_cor_sig[[i]]$zcor))
}

# subset based on 'strong' correlation > |0.5|
tissue_ES_cor_selected <- lapply(tissue_ES_cor_sig,function(x)x[abs(x$zcor)>0.5,])
tissue_random_cor_selected <- lapply(tissue_random_cor_sig,function(x)x[abs(x$zcor)>0.5,])


apply(cbind(sapply(tissue_ES_cor_selected,function(x) dim(x)[1]),sapply(tissue_random_cor_selected, function(y)dim(y)[1])),1,function(value)barplot(value, names.arg=c('Essential genes',"random genes"), main=names(value)))

# are there more nodes than random? OR node degree
par(mfrow=c(3,6))
for(dat in 1:length(tissue_ES_cor_sig)){
  plot_node_degree(tissue_ES_cor_sig[dat],tissue_random_cor_sig[dat])
}

# cell type level
celltype_ES_cor <- lapply(celltype_normed_list, get_ES_cor)
celltype_random_cor <-  lapply(celltype_normed_list, get_cor_random)
# Select only significant correlations 
celltype_ES_cor_sig <- lapply(celltype_ES_cor, function(x){
  x <- x[x$p<0.05,]
  x <- x[!is.na(x$zcor),]
  x <- x[!is.infinite(x$zcor),]})
celltype_random_cor_sig <- lapply(celltype_random_cor, function(x){
  x <- x[x$p<0.05,]
  x <- x[!is.na(x$zcor),]
  x <- x[!is.infinite(x$zcor),]})
# are there more significant terms than random?
dev.off()
par(mfrow=c(5,6))

apply(cbind(sapply(celltype_ES_cor_sig,function(x) dim(x)[1]),sapply(celltype_random_cor_sig, function(y)dim(y)[1])),1,function(value)barplot(value, names.arg=c('Essential genes',"random genes"), main=names(value)))
apply(cbind(sapply(celltype_ES_cor_sig,function(x) dim(x)[1]),sapply(celltype_random_cor_sig, function(y)dim(y)[1])),1, chisq.test)

dev.off()
par(mfrow=c(5,6))
for (i in 1:length(celltype_ES_cor_sig)){
  boxplot(as.numeric(celltype_ES_cor_sig[[i]]$cor)
          ,as.numeric(celltype_random_cor_sig[[i]]$cor), names=c("ES genes","Random genes"))
  print(t.test(as.numeric(celltype_ES_cor_sig[[i]]$cor)
               ,as.numeric(celltype_random_cor_sig[[i]]$cor))$p.value)
}

# # subset based on 'strong' correlation > |0.5|
# tissue_ES_cor_selected <- lapply(tissue_ES_cor_sig,function(x)x[abs(x$zcor)>0.5,])
# tissue_random_cor_selected <- lapply(tissue_random_cor_sig,function(x)x[abs(x$zcor)>0.5,])
# 
# 
# apply(cbind(sapply(tissue_ES_cor_selected,function(x) dim(x)[1]),sapply(tissue_random_cor_selected, function(y)dim(y)[1])),1,function(value)barplot(value, names.arg=c('Essential genes',"random genes"), main=names(value)))

# are there more nodes than random? OR node degree
dev.off()
par(mfrow=c(5,6))
for(dat in 1:length(celltype_ES_cor_sig)){
  plot_node_degree(celltype_ES_cor_sig[dat],
                   celltype_random_cor_sig[dat])
}

# correlations with Bayesian corrrelation######
 # how to determine the significance?
# ES genes within pathways ######
# For Hallmark 
ES_per_pway_hallmark <- ES_per_pathway(ES_genes_mmu_final_update,database = "Hallmark")
hallmark_size <- c(table(msigdbr(species = "Mus musculus", category = "H")%>% dplyr::select(gs_name)))
ES_per_pway_hallmark_pct <- sapply(ES_per_pway_hallmark,function(x)length(x))/hallmark_size*100
allgenes <- Reduce(intersect,lapply(tissue_normed_list, function(x)rownames(x)) )
allgenes_random_hallmark <- pathway_random(allgenes,database = "Hallmark")
allgenes_per_pway_hallmark_pct <- sapply(allgenes_random_hallmark,function(x)length(x))/hallmark_size*100

boxplot(ES_per_pway_hallmark_pct,allgenes_per_pway_hallmark_pct, notch = T,names=c('Essential genes',"random genes" ), main="Hallmark pathways")
plot(ES_per_pway_hallmark_pct,allgenes_per_pway_hallmark_pct, pch =19,xlim=c(0,35),ylim=c(0,35), xlab="Percentage of ES genes per pathway (Hallmark)", ylab="Percentage of random genes per pathway (Hallmark)")

# For KEGG
ES_per_pway_kegg <- ES_per_pathway(ES_genes_mmu_final_update,database = "KEGG")
kegg_size <- c(table(msigdbr(species = "Mus musculus", category = "C2",subcategory = "CP:KEGG")%>% dplyr::select(gs_name)))
ES_per_pway_kegg_pct <- sapply(ES_per_pway_kegg,function(x)length(x))/kegg_size*100
allgenes_random_kegg <- pathway_random(allgenes,database = "KEGG")
allgenes_per_pway_kegg_pct <- sapply(allgenes_random_kegg,function(x)length(x))/kegg_size*100
boxplot(ES_per_pway_kegg_pct,allgenes_per_pway_kegg_pct, notch = T,names=c('Essential genes',"random genes" ), main="KEGG pathways")
plot(ES_per_pway_kegg_pct,allgenes_per_pway_kegg_pct, pch =19,xlim=c(0,70),ylim=c(0,70), xlab="Percentage of ES genes per pathway (KEGG)", ylab="Percentage of random genes per pathway (KEGG)")


# Chromosomal location #####
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
getBM(attributes = c('affy_hg_u133_plus_2', 'hgnc_symbol', 'chromosome_name',
                     'start_position', 'end_position', 'band'),
      filters = 'affy_hg_u133_plus_2', 
      values = affyids, 
      mart = ensembl)


# functions ########
get_ES_rank <- function(full_metric){
  metric_rank <- rank(rowMeans(as.matrix(full_metric@assays$SCT@data)))
  ES_rank <- metric_rank[names(metric_rank)%in%ES_genes_mmu_final_update1]
  ES_rank <- ES_rank/dim(full_metric)[1]*100
}
get_ES_pct <- function(metric){
  pct <- apply(metric,1, function(x) sum(x>0))
  pct <- pct/dim(metric)[2]
}
get_ES_cor <- function(full_metric){
  met <-  as.matrix(full_metric@assays$SCT@data)
  ES_met <- met[rownames(met)%in%ES_genes_mmu_final_update1,]
  library(Hmisc)
  res2 <- rcorr(as.matrix(t(ES_met)))
  cor_met <- flattenCorrMatrix(res2$r, res2$P)
  cor_met$zcor <- atanh(cor_met$cor)
  cor_met
}
get_ES_cor1 <- function(full_metric){
  met <-  as.matrix(full_metric@assays$SCT@data)
  ES_met <- met[rownames(met)%in%ES_genes_mmu_final_update1,]
  #compute the Bayesian correlation matrix
  B <- BaCo(ES_met)   
  ut <- upper.tri(B)
  data.frame(
    row = rownames(B)[row(B)[ut]],
    column = rownames(B)[col(B)[ut]],
    cor  =(B)[ut] 
  )
  
}
get_cor_random <- function(full_metric){
  set.seed(123441)  
  random_sample <- sample(dim(full_metric)[1],length(ES_genes_mmu_final_update1))
  random_metric <- subset(full_metric,features=rownames(full_metric)[random_sample])
  random_metric <- as.matrix(random_metric@assays$SCT@data)
  library(Hmisc)
  res  <- rcorr(as.matrix(t(random_metric)))
  cor_met <- flattenCorrMatrix(res$r, res$P)
  cor_met$zcor <- atanh(cor_met$cor)
  cor_met
  
}
plot_mean_exp <- function(full_metric){
  set.seed(123441)  
  ES_metric <- as.matrix(full_metric@assays$SCT@data)
  ES_metric <- ES_metric[rownames(ES_metric)%in%ES_genes_mmu_final_update1,]
  random_sample <- sample(dim(full_metric)[1],dim(ES_metric)[1])
  random_metric <- subset(full_metric,features=rownames(full_metric)[random_sample])
  boxplot(rowMeans(ES_metric),rowMeans(as.matrix(random_metric@assays$SCT@data)),notch = T,names=c('scEssentials',"Random"),ylab = 'Log(mean)', main = paste0(unique(full_metric$tissue),'.',unique(full_metric$cell_ontology_class)), cex.main=0.8)
  t.test(rowMeans(ES_metric),rowMeans(as.matrix(random_metric@assays$SCT@data)))$p.value
  
}
plot_pct <- function(full_metric){
  set.seed(123441) 
  ES_metric <- as.matrix(full_metric@assays$SCT@data)
  ES_metric <- ES_metric[rownames(ES_metric)%in%ES_genes_mmu_final_update1,]
  pct <- get_ES_pct(ES_metric)
  
  random_sample <- sample(dim(full_metric)[1],dim(ES_metric)[1])
  random_metric <- subset(full_metric,features=rownames(full_metric)[random_sample])
  random_metric <- as.matrix(random_metric@assays$SCT@data)
  
  pct.random <- get_ES_pct(random_metric)
  boxplot(pct,pct.random,notch = T,names=c('scEssential',"Random"),ylab = 'Percentage of expression', main = paste0(unique(full_metric$tissue),'.',unique(full_metric$cell_ontology_class)), cex.main=0.8)
  wilcox.test(pct,pct.random)$p.value
  
}
plot_node_degree <- function(ES, random){
  library(igraph)
  df.g <- graph.data.frame(d = ES, directed = FALSE)
  df.g1 <- graph.data.frame(d = random, directed = FALSE)
  boxplot(degree(df.g), degree(df.g1), names=c('scEssentials',"Random"))
  print(t.test(degree(df.g), degree(df.g1))$p.value)
}
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  = (cormat)[ut],
    p = pmat[ut]
  )
}
BaCo <- function(X){
  
  alpha0 <- rep(1/nrow(X),ncol(X))
  beta0=1-alpha0
  nrowsX <- nrow(X)
  k <- ncol(X)
  cs <- colSums(X)
  alphas <- alpha0 + X
  betas  <- matrix(rep(beta0,nrowsX), nrow=nrowsX, byrow=TRUE) + matrix(rep(cs,nrowsX), nrow=nrowsX, byrow=TRUE) - X
  alphasPLUSbetas <- alphas + betas
  Psi <- alphas/alphasPLUSbetas - matrix(rep(rowSums(alphas/alphasPLUSbetas)/k, k), ncol=k, byrow=FALSE) 
  var_vec <- as.matrix( ( rowSums( (alphas*betas)/( (alphasPLUSbetas^2)*(alphasPLUSbetas+1) ) ) + rowSums(Psi^2) )/k )
  cov_mtrx <- (Psi %*% t(Psi))/k
  Bcorrvals <- cov_mtrx / sqrt( var_vec %*% t(var_vec) )
  diag(Bcorrvals) <- 1
  Bcorrvals
}
ES_per_pathway <- function(ES_list, database){
  if (database=="Hallmark"){
    hallmarks <-  msigdbr(species = "Mus musculus", category = "H")%>% dplyr::select(gs_name,gene_symbol)
    pways_list <- list()
    for (pways in unique(hallmarks$gs_name)){
      pways_list[[pways]] <- hallmarks[hallmarks$gs_name== pways,2]
    }
    overlap <- lapply(pways_list, function(x){x[["gene_symbol"]][x[["gene_symbol"]]%in%ES_list]})
    overlap
  }
  if (database=='KEGG'){
    keggs <-  msigdbr(species = "Mus musculus", category = "C2",subcategory = "CP:KEGG")%>% dplyr::select(gs_name,gene_symbol)
    pways_list <- list()
    for (pways in unique(keggs$gs_name)){
      pways_list[[pways]] <- keggs[keggs$gs_name== pways,2]
    }
    overlap <- lapply(pways_list, function(x){x[["gene_symbol"]][x[["gene_symbol"]]%in%ES_list]})
    overlap
    
  }
  overlap
}
pathway_random <- function(full_gene_name,database){
  set.seed(123441)
  random_sample <- sample(full_gene_name,length(ES_genes_mmu_final_update1))
  random_metric_per_pways <- ES_per_pathway(random_sample,database)
  random_metric_per_pways
}
