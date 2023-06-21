library(sceasy)
library(reticulate)
library(plyr)
library(stringr)
library(biomaRt)
library(UpSetR)
library(genefilter)
library(corrplot)
library(ggfortify)
library(reshape2)
library(ggridges)
library(ggplot2)
library(dplyr)
library(msigdbr)
library(org.Hs.eg.db)
library(caret)
library(car)
memory.limit(23000*1024*1024)
# data preparation #####
# subset based on immune - for five donors and FACs
x1@assays$RNA@key <- "rna_"
x1=subset(tabula_sapiens_immune,subset=donor%in%c("TSP4","TSP5","TSP9","TSP10","TSP13"))
x1 = subset(x1,subset=assay=="Smart-seq2")

# cell types from endothelia, epithelia and stromal cell types
# extract more than 100 cell tyoes
TS_all <- merge(immune_sub,y=c(endo_selected,epith_selected,stromal_selected))
# convert from ENSEMBLE to gene symbol
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# scTransform normalisation and convert ensembl id to gene name
TS_all <- convert_norm(TS_all)
#saveRDS(TS_all,"TS_merge_all.rds")
info <- as.data.frame(table(TS_all$tissue, TS_all$cell_ontology_class))
info.sub <- info[info$Freq>100,]

tissue_list <- list()
for (tissue in unique(info.sub$Var1)){
  cell.type <- droplevels(info.sub[info.sub$Var1==tissue,2])
  tissue_list[[tissue]] <- TS_all[,TS_all$tissue==tissue & TS_all$cell_ontology_class %in% cell.type]
  
}

par(mfrow=c(4,4))
lapply(tissue_list,function(x)boxplot(rowMeans(as.matrix(x@assays$SCT@data)), main=paste(x$tissue[1],x$cell_ontology_class[1],sep='.')))

celltype_list <- list()
for (tissue in names(tissue_list)){
  for (type in unique(tissue_list[[tissue]]$cell_ontology_class)){
    name <- paste(tissue,type,sep='.' )
    dat <- subset(tissue_list[[tissue]],subset=cell_ontology_class==type)
    celltype_list[[name]] <- dat
    }
}

# extract for ES genes 
ES_genes_hsa <- read.csv( "essential_genes_hsa_final.csv")
colnames(ES_genes_hsa) <- c("Name","Rank")

tissue_ES_list <- lapply(tissue_list,function(data){
  data <- as.matrix(data@assays$SCT@data)
  data <- data[rownames(data)%in%ES_genes_hsa$Name,]
})
barplot(sapply(tissue_ES_list,function(x) sum(rowSums(x)>0)))

# how many genes shared in the end?
# 1. able to detect at celltype  level
celltype_list.1 <- lapply(celltype_list,function(data){
  data <- as.matrix(data@assays$SCT@data)
  index <-  apply(data,1,function(y)sum(y>0))
  data <- data[index>3,]
  data <- data[rownames(data)%in%ES_genes_hsa$Name,]
})
upset(fromList(celltype_list.1), order.by = "freq",nsets = 8)
 
# ES genes are detectable
barplot(sapply(celltype_list.1,function(x) sum(rowSums(x)>0))[order(sapply(celltype_list.1,function(x) sum(rowSums(x)>0)),decreasing = T)], names.arg  = '',main = 'Number of essential genes across 53 cell types (TS)')
median(sapply(celltype_list.1,nrow))
ES_genes_hsa_final <- Reduce(intersect,lapply(celltype_list.1,rownames)) #2048

# # 2. less variation at cell type level 
celltype_ES_list <- lapply(celltype_list,function(data){
  data <- as.matrix(data@assays$SCT@data)
  data <- data[rownames(data)%in%ES_genes_hsa_final,]
})
celltype_ES_tbl <-  do.call(cbind,lapply(celltype_ES_list,rowMeans))
hist(rowSds(celltype_ES_tbl),xlab = '', main='Standard deviation of the distribution of the ranks across cell types (TS)')
ES_genes_hsa_final_update <- ES_genes_hsa_final[rowSds(celltype_ES_tbl)< 4*sd(rowSds(celltype_ES_tbl))] #1969

write.csv(ES_genes_hsa_final_update,'essential_genes_hsa_final.csv')
saveRDS(ES_genes_hsa_final_update,'scEssential_hsa.rds')

# overall ES genes have higher exp #######
#update the data with only 1969 genes
tissue_list_final <- lapply(tissue_ES_list,function(data){
  data <- data[rownames(data)%in%ES_genes_hsa_final_update,]
  })
celltype_list_final <- lapply(celltype_ES_list,function(data){
  data <- data[rownames(data)%in%ES_genes_hsa_final_update,]
}) 
par(mfrow=c(4,4))
sapply(tissue_list, plot_mean_exp)
par(mfrow=c(6,9))
sapply(celltype_list, plot_mean_exp)
par(mfrow=c(5,6))
sapply(celltype_list[1:30], plot_mean_exp)#1600*1000
par(mfrow=c(5,6))
sapply(celltype_list[31:53], plot_mean_exp)#1600*1000
# correlation and heatmaps ######
tissue_ES_tbl <-  do.call(cbind,lapply(tissue_list_final,rowMeans))

tissue_mat <- cor(tissue_ES_tbl)
# distribution of the correlations between tissue
corrplot(tissue_mat, method = 'circle',tl.pos = 'l',addCoef.col = 'black', tl.col = 'black',order = 'hclust',cl.pos = 'b',cl.ratio = 0.1, tl.cex = 0.8, col = colorRampPalette(c("blue","white","yellow","red"))(100))

heatmap(t(tissue_ES_tbl),scale ='row')

# distribution of the correlations between cell types
celltype_ES_tbl <-  do.call(cbind,lapply(celltype_list_final,rowMeans))
celltype_mat <- cor(celltype_ES_tbl)
# distribution of the correlations between cell types
hist(celltype_mat[upper.tri(celltype_mat)])
colnames(celltype_mat) <- rep(" ", NROW(celltype_mat))
corrplot(celltype_mat, method = 'circle',tl.pos = 'l',tl.col = 'black',order = 'hclust',cl.pos = 'b',cl.ratio = 0.1, tl.cex = 0.6, col = colorRampPalette(c("blue","white","yellow","red"))(100))
median(celltype_mat)
autoplot(prcomp(t(tissue_ES_tbl)),label = TRUE, label.size = 3)
autoplot(prcomp(t(celltype_ES_tbl)),label = TRUE, label.size = 3)

# Percentage of expression  #########
tissue_ES_pct <- sapply(tissue_list_final,get_ES_pct)
tissue_ES_pct.1 <- melt(tissue_ES_pct)
gene_pct_order <- rowMeans(tissue_ES_pct)[order(rowMeans(tissue_ES_pct),decreasing = T)]
ggplot(tissue_ES_pct.1, aes(x = value, y = Var2,fill=as.factor(Var2))) + geom_density_ridges()+theme_bw()
ggplot(tissue_ES_pct.1[tissue_ES_pct.1$Var1%in%names(gene_pct_order)[1:50],], aes(x = value, y = Var1)) + geom_boxplot()+theme_bw()


celltype_ES_pct <- sapply(celltype_list_final,get_ES_pct)
celltype_ES_pct.1 <- melt(celltype_ES_pct)
gene_pct_order1 <- rowMeans(celltype_ES_pct)[order(rowMeans(celltype_ES_pct),decreasing = T)]

ggplot(celltype_ES_pct.1, aes(x = value, y = Var2,fill=as.factor(Var2))) + geom_boxplot()+theme_bw()
ggplot(celltype_ES_pct.1[celltype_ES_pct.1$Var1%in%names(gene_pct_order1)[1:50],], aes(x = value, y = Var1)) + geom_boxplot()+theme_bw()
 
# more pct than random?
par(mfrow=c(4,4))
sapply(tissue_list, plot_pct)

par(mfrow=c(5,6))
sapply(celltype_list[1:30], plot_pct)#1600*1000
celltype_GDI
sapply(celltype_list[31:53], plot_pct)#1600*1000


# Marker genes ########
# some of the tissues does not contain more than two cell types 
tissue_markers <- lapply(tissue_list,function(x){
  Idents(x) <- x$cell_ontology_class
  markers <- FindAllMarkers(x, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
})
# for now focus on cell types that have markers 
tissue_markers.1969 <- lapply(tissue_markers[c(1:7,9,11:13)],function(x){x%>%
    group_by(cluster) %>%
    slice_max(n = 1969, order_by = avg_log2FC)%>% 
    dplyr::select(gene) })

celltype_markers_list <- list()
for (tissue in names(tissue_markers.1969)){
  for (type in unique(tissue_markers.1969[[tissue]]$cluster)){
    name <- paste(tissue,type,sep='.' )
    dat <-  tissue_markers.1969[[tissue]][tissue_markers.1969[[tissue]]$cluster==type,'gene']
    celltype_markers_list[[name]] <- dat}
}

ES_markers_overlap <- sapply(celltype_markers_list,function(x)intersect(x[["gene"]],ES_genes_hsa_final_update))
barplot(sapply(ES_markers_overlap,length), horiz = T)
boxplot(sapply(ES_markers_overlap,length))
upset(fromList(ES_markers_overlap), nsets = 10)

ES_markers_overlap.1 <- as.data.frame(sapply(ES_markers_overlap,length))
colnames(ES_markers_overlap.1) <- 'value'
ES_markers_overlap.1$celltype <- rownames(ES_markers_overlap.1)
ES_markers_overlap.1$type <- sapply(ES_markers_overlap.1$celltype,function(cx)str_split(cx,"\\.")[[1]][1])

ggplot(ES_markers_overlap.1, aes(x = celltype, y = value ,fill=type )) + geom_bar(stat="identity")+theme_bw()+ coord_flip()+ggtitle("Number of scEssentials overlapped with cell-type markers")+labs(fill='Source tissue')

ggplot(ES_markers_overlap.1, aes(x=celltype, y=value/1969,group = 1)) + geom_line() 

median(sapply(ES_markers_overlap,length)/1969)

# when compared to large marker database
library(data.table)
Panglao <- as.data.frame(fread("PanglaoDB_markers_27_Mar_2020.tsv"))
human_Panglao <- Panglao[Panglao$species=='Hs'|Panglao$species=='Mm Hs',]
human_markers_Panglao <- unique(human_Panglao$`official gene symbol`)
library(ggvenn)
x = intersect(human_markers_Panglao,ES_genes_hsa_final_update)
human_sensitivity <- human_Panglao[human_Panglao$`official gene symbol`%in%x,]
human_sensitivity_score <- aggregate(human_sensitivity$sensitivity_human,by=list(human_sensitivity$`official gene symbol`) , max)
colnames(human_sensitivity_score) = c('gene','sensitivity')
Number_of_panglao_ES <- sapply(celltype_list_final,function(o)length(intersect(rownames(o),human_markers_Panglao)))
median(Number_of_panglao_ES)/48
# compute the essentially score
xx = as.data.frame(ES_genes_hsa_final_update)
colnames(xx) = 'gene'
S = dplyr::full_join(xx,human_sensitivity_score,by='gene')
S[is.na(S)] <- 0 

num_marker =  unlist(celltype_markers_list)
num_marker <- table(num_marker)
Gx <- as.data.frame(num_marker)
colnames(Gx) <- c("gene","freq")
Gx <- Gx[Gx$gene%in%xx$gene,]
G = dplyr::full_join(xx,Gx, by='gene')
G[is.na(G)] <- 0 
G$freq <- G$freq/48

w = as.data.frame((S$sensitivity + G$freq)/2)
w$gene = S$gene

pct = rowMeans(unlist(celltype_ES_pct))*100

ES_score = (1-w$`(S$sensitivity + G$freq)/2`)*pct
ES_score <- as.data.frame(ES_score)
 
#HVG genes ########
library(Seurat)
celltype_HVGs <- list()
for (tissue in names(tissue_list)){
  for (type in unique(tissue_list[[tissue]]$cell_ontology_class)){
    name <- paste(tissue,type,sep='.' )
    dat <- subset(tissue_list[[tissue]],subset=cell_ontology_class==type)
    dat <- FindVariableFeatures(dat, assay = "RNA",selection.method = "vst", nfeatures = 1969)
    hvgs <-  head(VariableFeatures(dat,assay = "RNA",selection.method = "vst"), 1969)
    celltype_HVGs[[name]] <- hvgs}
}

ES_HVGs <- sapply(celltype_HVGs,function(s) intersect(s,ES_genes_hsa_final_update ))
barplot(unique(sapply(ES_HVGs,length)), horiz = T)
max(unique(sapply(ES_HVGs,length)))

ES_HVGs.1 <- as.data.frame(sapply(ES_HVGs,length))
colnames(ES_HVGs.1) <- 'value'
ES_HVGs.1$celltype <- rownames(ES_HVGs.1)
ES_HVGs.1$type <- sapply(ES_HVGs.1$celltype,function(cx)str_split(cx,"\\.")[[1]][1])

ggplot(ES_HVGs.1, aes(x = celltype, y = value ,fill=type )) + geom_bar(stat="identity")+theme_bw()+ coord_flip()+ggtitle("Number of scEssentials with HVGs")+labs(fill='Source tissue')

## combined DEG and HVG ######
ES_all_overlap <- merge(ES_markers_overlap.1,ES_HVGs.1, by="celltype")
#ES_all_overlap$gene_category <- c(rep('DEGs',63), rep("HVGs",68))
saveRDS(ES_all_overlap,"TS_overlap_markerHVG.rds")

library(ggalt)
library(scales)
ttt = as.data.frame(table(ES_all_overlap$type.y))
x_cols <- rep(hue_pal()(length(ttt$Var1)),ttt$Freq) 
p2 <- ggplot(ES_all_overlap, aes(y=celltype, x=value.x/1969, xend=value.y/1969)) +
  geom_dumbbell(size=3, color="#e3e2e1",
                colour_x = "#5b8124", colour_xend = "#bad744",
                dot_guide=TRUE, dot_guide_size=0.25) +
  labs(x=NULL, y=NULL, title="Overlap of scEssentials and markers/HVGs (human)") +
  theme_minimal()+
  theme(axis.text.y = element_text(),
        axis.ticks.y=element_blank())
p2
p1 <- readRDS('overlap_mmu_ggplot.rds')
library(ggpubr)
p1|p2
ggarrange(p1+
            labs(x=NULL, y=NULL, title="Overlap of scEssentials and markers/HVGs (Mouse)")+theme(axis.text.y = element_text(colour="black")) , p2+theme(axis.text.y = element_text(colour="black")),  
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

# SEG genes #######
# many source 
SEG_list <- list("SEGs"= scHK_human$`Gene Symbol`, "Housekeeping genes"=Housekeeping_TranscriptsHuman$Gene_symbol,'Cytosolic ribosomal genes'= Gene_Summaries[Gene_Summaries$CytosolicRibosome_indicator,'gene']$gene,'scEssentials'= ES_genes_hsa_final_update)
unique_essential <- intersect(Housekeeping_TranscriptsHuman$Gene_symbol,ES_genes_hsa_final_update)
library(clusterProfiler)
unique_essential1 <- bitr(unique_essential,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)
universe <- rownames(celltype_list[[1]])
universe1 <- bitr(universe,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene          = unique_essential1$ENTREZID,
                universe      = universe1$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                minGSSize = 5,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
R.utils::setOption("clusterProfiler.download.method","auto")
kk <- enrichKEGG(gene         = unique_essential1$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05 )

library(UpSetR)
upset(fromList(SEG_list), order.by = 'freq',nsets =4 ) 
library(ggvenn)
ggvenn(SEG_list,fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
       stroke_size = 0.5, set_name_size = 6)
# correlations #####
tissue_ES_cor <- lapply(tissue_list, get_ES_cor)
tissue_random_cor <-  lapply(tissue_list, get_cor_random)
# Select only significant correlations 
tissue_ES_cor_sig <- lapply(tissue_ES_cor, function(x){x[x$p<0.05,]})
tissue_random_cor_sig<- lapply(tissue_random_cor, function(x){x[x$p<0.05,]})
# are there more signfincant terms than random?
par(mfrow=c(4,4))

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
dev.off()
par(mfrow=c(4,4))
for(dat in 1:length(tissue_ES_cor_sig)){
  plot_node_degree(tissue_ES_cor_sig[dat],tissue_random_cor_sig[dat])
}

# cell type level
celltype_ES_cor <- lapply(celltype_list, get_ES_cor)
celltype_random_cor <-  lapply(celltype_list, get_cor_random)
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




# ES genes within pathways ######
# For Hallmark 
ES_per_pway_hallmark <- ES_per_pathway(ES_genes_hsa_final_update,database = "Hallmark")
hallmark_size <- c(table(msigdbr(species = "Homo sapiens", category = "H")%>% dplyr::select(gs_name)))
ES_per_pway_hallmark_pct <- sapply(ES_per_pway_hallmark,function(x)length(x))/hallmark_size*100
allgenes <- Reduce(intersect,lapply(tissue_list, function(x)rownames(x)) )
allgenes_random_hallmark <- pathway_random(allgenes,database = "Hallmark")
allgenes_per_pway_hallmark_pct <- sapply(allgenes_random_hallmark,function(x)length(x))/hallmark_size*100
dev.off()
boxplot(ES_per_pway_hallmark_pct,allgenes_per_pway_hallmark_pct, notch = T,names=c('Essential genes',"random genes" ), main="Hallmark pathways",ylab='Percentage of genes found in each hallmark pathway')
wilcox.test(ES_per_pway_hallmark_pct,allgenes_per_pway_hallmark_pct)
plot(ES_per_pway_hallmark_pct,allgenes_per_pway_hallmark_pct, pch =19,xlim=c(0,35),ylim=c(0,35), xlab="Percentage of ES genes per pathway (Hallmark)", ylab="Percentage of random genes per pathway (Hallmark)")
par(mfrow=c(1,2))
# For KEGG
ES_per_pway_kegg <- ES_per_pathway(ES_genes_hsa_final_update,database = "KEGG")
kegg_size <- c(table(msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")%>% dplyr::select(gs_name)))
ES_per_pway_kegg_pct <- sapply(ES_per_pway_kegg,function(x)length(x))/kegg_size*100
allgenes_random_kegg <- pathway_random(allgenes,database = "KEGG")
allgenes_per_pway_kegg_pct <- sapply(allgenes_random_kegg,function(x)length(x))/kegg_size*100
boxplot(ES_per_pway_kegg_pct,allgenes_per_pway_kegg_pct, notch = T,names=c('Essential genes',"random genes" ), main="KEGG pathways", ylab='Percentage of genes in each KEGG pathway')

plot(ES_per_pway_kegg_pct,allgenes_per_pway_kegg_pct, pch =19,xlim=c(0,70),ylim=c(0,70), xlab="Percentage of ES genes per pathway (KEGG)", ylab="Percentage of random genes per pathway (KEGG)")


# gene damage index with gene expression variation ####
gdi <- fread("Gene damage index/Gene damage index.txt")
GDI_es_low <- GDI_es[GDI_es$`Gene damage prediction (all disease-causing genes)`=="Low",]
GDI_es_Medium <- GDI_es[GDI_es$`Gene damage prediction (all disease-causing genes)`=="Medium",]
GDI_es_High <- GDI_es[GDI_es$`Gene damage prediction (all disease-causing genes)`=="High",]

# chi-square test 
Numbers <- as.data.frame(rbind(c(91,336),c(1728,16596), c(39,768)))
chisq.test(Numbers)
colnames(Numbers) <- c("scEssential","Non-scEssential")
rownames(Numbers) <- c('Low',"Medium","High")
Numbers.1 <- melt(Numbers)
Numbers.1$level <- rep(rownames(Numbers),2)
ggplot(data=Numbers.1, aes(x=variable, y=value, fill=level)) + geom_bar(stat="identity",position="fill")


# measuring the gene expression variability 
library(scran)
celltype_list_EV <- lapply(celltype_list,function(w)get_var(w@assays$SCT@counts))
celltype_GDI_low <- sapply(celltype_list_EV,function(mat)x2= mat[names(mat)%in%GDI_es_low$Gene])
celltype_GDI_medium <- sapply(celltype_list_EV,function(mat)x2= mat[names(mat)%in%GDI_es_Medium$Gene])
celltype_GDI_high <- sapply(celltype_list_EV,function(mat)x2= mat[names(mat)%in%GDI_es_High$Gene])

celltype_GDI <- rbind(celltype_GDI_low,celltype_GDI_medium, celltype_GDI_high)

library(caret)
par(mfrow=c(5,6))
celltype_GDI_list <- list()
for (i in 1:10){
  celltype_GDI_p <- apply(celltype_GDI, 2,function(ctype){
    ctype_df <- melt(ctype)
    ctype_df$group <- c(rep('low',111),rep('medium',4203), rep('high',112))
    ctype_df$group <- as.factor(ctype_df$group)
    down_df <- downSample(ctype_df$value, ctype_df$group)
    #boxplot(down_df[down_df$Class=='low',1],down_df[down_df$Class=='medium',1],down_df[down_df$Class=='high',1], names=c("Low","Medium","High"),ylab="Gene expression variability")
    kruskal.test(x~Class, data = down_df)$p.value })
  celltype_GDI_p.adj <- p.adjust(celltype_GDI_p)
  celltype_GDI_list[[i]] <- celltype_GDI_p.adj
}

# celltype_GDI_p <- apply(celltype_GDI, 2,function(ctype){
#   ctype_df <- melt(ctype)
#   ctype_df$group <- c(rep('low',111),rep('medium',4315), rep('high',112))
#   ctype_df$group <- as.factor(ctype_df$group)
#   down_df <- downSample(ctype_df$value, ctype_df$group)
#   boxplot(down_df[down_df$Class=='low',1],down_df[down_df$Class=='medium',1],down_df[down_df$Class=='high',1], names=c("Low","Medium","High"),ylab="Gene expression variability")
#   leveneTest(x~Class, data = down_df)[[3]][1]$p.value })
# celltype_GDI_p.adj <- p.adjust(celltype_GDI_p)
#celltype_GDI_p1 <- apply(celltype_GDI, 2,function(ctype)levene.test(ctype[1:283], ctype[284:3699])$estimate)
aa <- Reduce(rbind,celltype_GDI_list)
set.seed(123411)
par(mfrow=c(5,6),
    cex.axis=0.8)
GDI.ps <- c()
for (n in 1:ncol(celltype_GDI)){
  ctype = celltype_GDI[,n]
  ctype_df <- reshape2::melt(ctype)
  ctype_df$group <- c(rep('low',111),rep('medium',4203), rep('high',112))
  ctype_df$group <- as.factor(ctype_df$group)
  down_df <- downSample(ctype_df$value, ctype_df$group)
  boxplot(down_df[down_df$Class=='low',1],down_df[down_df$Class=='medium',1],down_df[down_df$Class=='high',1], names=c("Low","Medium","High"),ylab="Gene expression variability", main=colnames(celltype_GDI)[n],cex.main=0.8)
  GDI.ps <- c(GDI.ps,kruskal.test(x~Class, data = down_df)$p.value)
}
 
GDI.ps.adj <- p.adjust(GDI.ps)
sum(GDI.ps.adj<0.05)

#GDI var all
GDI_var_all <- cbind(colMeans(celltype_GDI_low),
                     colMeans(celltype_GDI_medium),
                     colMeans(celltype_GDI_high))
colnames(GDI_var_all) <- c("Low","Medium","High")
GDI_var_all_1 <- melt(GDI_var_all)
library(ggpubr)
p_GDI <- ggplot(GDI_var_all_1, aes(x=Var2, y=value, fill=Var2)) + 
  geom_violin()+xlab('')+ylab("Gene expression variability")+scale_fill_manual(values=c("#D35E40", "#D3A840", "#B5D340"))+theme_classic() +  theme(legend.position="none")+ stat_compare_means(label.y = 0.8)+ theme(text = element_text(size = 16))  +
  stat_compare_means(label = "p.signif", method = "wilcox.test",ref.group = "Low")+xlab('Gene damage prediction (all disease-causing genes)')+ geom_boxplot(width=0.1)
# what are the top variable genes??
#top_genes <- apply(celltype_GDI_medium, 2,rank)
top_genes_avg <- rowMeans(celltype_GDI) 
boxplot(top_genes_avg[1:111],top_genes_avg[112:4426],top_genes_avg[4427:4538])
kruskal.test(top_genes_avg[1:111],top_genes_avg[112:4426],top_genes_avg[4427:4538])
top_genes_avg <- top_genes_avg[order(top_genes_avg, decreasing = T)]
xx = GDI_es[match(names(top_genes_avg),GDI_es$Gene),]
plot(xx$`GDI-Phred`,top_genes_avg, pch = 19,cex=0.5, col=as.factor(xx$`Gene damage prediction (all disease-causing genes)`))

top_genes_avg[1:10]
top_genes_avg[4528:4538]

# variability of other high damage genes
library(data.table)
gdi <- fread("Gene damage index/Gene damage index.txt")
gdi_high <- gdi[gdi$`Gene damage prediction (all disease-causing genes)`=='High','Gene']
GDI_other_high <- sapply(celltype_list_EV,function(mat)x2= mat[names(mat)%in%setdiff(gdi_high$Gene,GDI_es_High$Gene)])
par(mfrow=c(5,6),
    cex.axis=0.8)
p.others <- c()
for(i in 1:53){
  aaa=wilcox.test(celltype_GDI_high[,i],GDI_other_high[,i],alternative = 'less')$p.value
  p.others=c(p.others,aaa)
  boxplot(celltype_GDI_high[,i],GDI_other_high[,i], names=c("scEssentials","Non-scEssentials"), ylab = "Gene expression variability", main=names(celltype_list_EV)[i],cex.main=0.8)
}

sum(colMeans(celltype_GDI_high)<colMeans(GDI_other_high))
mat = cbind(colMeans(celltype_GDI_high),colMeans(GDI_other_high),p.adjust(p.others))

p_GDI1 <- ggplot(GDI_ES_score.df,aes(x=(`Gene damage prediction (all disease-causing genes)`), y=log2(ES_score),fill=(`Gene damage prediction (all disease-causing genes)`))) + geom_violin()+scale_fill_manual(values=c("#D35E40", "#D3A840", "#B5D340")) + xlab("Gene damage prediction (all disease-causing genes)")+ylab("Essentiality score")+theme_classic()+theme(legend.position="none")+stat_compare_means(label.y = 7)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",ref.group = "Low")+ theme(text = element_text(size = 16))+ geom_boxplot(width=0.1)
ggarrange(p_GDI,p_GDI1,labels = 'AUTO')

# does it related to the essentially score?###
ES_score$Gene <- rownames(ES_score)
GDI_ES_score.df <- dplyr::left_join(ES_score,GDI_es, by='Gene')
GDI_ES_score.df <- na.omit(GDI_ES_score.df)
write.csv(ES_score,'scEssential_hsa.csv')
saveRDS(ES_score,'scEssential_hsa.rds')
# correlation

cor.test(log2(GDI_ES_score.df$ES_score), GDI_ES_score.df$`GDI-Phred`, method = 'spearman')

plot(log2(GDI_ES_score.df$ES_score),log1p(GDI_ES_score.df$`GDI-Phred`))

# Default plot
library(ggplot2)
ggplot(GDI_ES_score.df, aes(x=log2(ES_score), y=log1p(`GDI-Phred`)))+ geom_density_2d_filled(alpha = 0.8, colour = "black")+ theme_classic() +theme(legend.position="none")+ylab("Normalised GDI score") + xlab('Essentiality score')
GDI_ES_score.df$`Gene damage prediction (all disease-causing genes)` <- as.factor(GDI_ES_score.df$`Gene damage prediction (all disease-causing genes)`)
GDI_ES_score.df$`Gene damage prediction (all disease-causing genes)` <- factor(GDI_ES_score.df$`Gene damage prediction (all disease-causing genes)`,levels = c("Low", "Medium", "High"))
 
kruskal.test(log2(ES_score)~`Gene damage prediction (all disease-causing genes)`, data = GDI_ES_score.df)$p.value
wilcox.test(GDI_ES_score.df[GDI_ES_score.df$`Gene damage prediction (all disease-causing genes)`=='High','ES_score'],GDI_ES_score.df[GDI_ES_score.df$`Gene damage prediction (all disease-causing genes)`=='Medium','ES_score'])

# ES score vs variability 
ES_GDI_EV <- sapply(celltype_list_EV,function(s)s[names(s)%in%GDI_ES_score.df$Gene])
par(mfrow=c(5,6),
    cex.axis=0.8)
GDI_cor <- c()
for(i in 1:53){
  aaa=cor(ES_GDI_EV[,i],GDI_ES_score.df$ES_score) 
  GDI_cor=c(GDI_cor,aaa)
  plot(ES_GDI_EV[,i],GDI_ES_score.df$ES_score, xlab ="Gene expression variability",ylab='Essentiality score' ,main=names(celltype_list_EV)[i],cex.main=0.8)
}
boxplot(GDI_cor, main= ' ')
median(GDI_cor)
#functions #####
convert_norm <- function(mat){
  genes <- rownames(mat)
  gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                    values = genes, mart= mart)
  gene_IDs <- na.omit(gene_IDs)
  gene_IDs <- gene_IDs[gene_IDs$hgnc_symbol!='',]
  gene_IDs <- gene_IDs[!duplicated(gene_IDs$ensembl_gene_id),]
  gene_IDs <- gene_IDs[!duplicated(gene_IDs$hgnc_symbol),]

  mat <- subset(x = mat, features = gene_IDs$ensembl_gene_id )
  mat1 <- GetAssayData(mat, slot = 'count')
  rownames(mat1) <-   gene_IDs[match(rownames(mat),gene_IDs$ensembl_gene_id),2]
  sobj <- CreateSeuratObject(counts=mat1, meta.data = mat@meta.data)
  sobj <- SCTransform(sobj, vars.to.regress = "donor",verbose=TRUE)

}
plot_mean_exp <- function(full_metric){
  set.seed(123441)  
  ES_metric <- as.matrix(full_metric@assays$SCT@data)
  ES_metric <- ES_metric[rownames(ES_metric)%in%ES_genes_hsa_final_update,]
  random_sample <- sample(dim(full_metric)[1],dim(ES_metric)[1])
  random_metric <- subset(full_metric,features=rownames(full_metric)[random_sample])
  boxplot(rowMeans(ES_metric),rowMeans(as.matrix(random_metric@assays$SCT@data)),names=c('scEssential',"Random"),ylab = 'Log(mean)', main = paste0(unique(full_metric$tissue),'.',unique(full_metric$cell_ontology_class)), cex.main=0.8)
  wilcox.test(rowMeans(ES_metric),rowMeans(as.matrix(random_metric@assays$SCT@data)))$p.value
  
}
get_ES_pct <- function(metric){
  pct <- apply(metric,1, function(x) sum(x>0))
  pct <- pct/dim(metric)[2]
}
plot_pct <- function(full_metric){
  set.seed(123441) 
  ES_metric <- as.matrix(full_metric@assays$SCT@data)
  ES_metric <- ES_metric[rownames(ES_metric)%in%ES_genes_hsa_final_update,]
  pct <- get_ES_pct(ES_metric)
  
  random_sample <- sample(dim(full_metric)[1],dim(ES_metric)[1])
  random_metric <- subset(full_metric,features=rownames(full_metric)[random_sample])
  random_metric <- as.matrix(random_metric@assays$SCT@data)
  
  pct.random <- get_ES_pct(random_metric)
  boxplot(pct,pct.random,names=c('scEssential',"Random"),ylab = 'Percentage of expression', main = paste0(unique(full_metric$tissue),'.',unique(full_metric$cell_ontology_class)), cex.main=0.8)
  wilcox.test(pct,pct.random)$p.value
  
}
get_ES_cor <- function(full_metric){
  met <-  as.matrix(full_metric@assays$SCT@data)
  ES_met <- met[rownames(met)%in%ES_genes_hsa_final_update,]
  library(Hmisc)
  res2 <- rcorr(as.matrix(t(ES_met)))
  cor_met <- flattenCorrMatrix(res2$r, res2$P)
  #cor_met$zcor <- atanh(cor_met$cor)
  cor_met
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
get_cor_random <- function(full_metric){
  set.seed(123441)  
  random_sample <- sample(dim(full_metric)[1],length(ES_genes_hsa_final_update))
  random_metric <- subset(full_metric,features=rownames(full_metric)[random_sample])
  random_metric <- as.matrix(random_metric@assays$SCT@data)
  library(Hmisc)
  res  <- rcorr(as.matrix(t(random_metric)))
  cor_met <- flattenCorrMatrix(res$r, res$P)
  cor_met$zcor <- atanh(cor_met$cor)
  cor_met
  
}
plot_node_degree <- function(ES, random){
  library(igraph)
  df.g <- graph.data.frame(d = ES, directed = FALSE)
  df.g1 <- graph.data.frame(d = random, directed = FALSE)
  boxplot(log(degree(df.g)), log(degree(df.g1)), names=c('ES genes',"Random genes"), ylab="Log(Network node degree)")
  print(wilcox.test(degree(df.g), degree(df.g1))$p.value)
}
ES_per_pathway <- function(ES_list, database){
  if (database=="Hallmark"){
    hallmarks <-  msigdbr(species = "Homo sapiens", category = "H")%>% dplyr::select(gs_name,gene_symbol)
    pways_list <- list()
    for (pways in unique(hallmarks$gs_name)){
      pways_list[[pways]] <- hallmarks[hallmarks$gs_name== pways,2]
    }
    overlap <- lapply(pways_list, function(x){x[["gene_symbol"]][x[["gene_symbol"]]%in%ES_list]})
    overlap
  }
  if (database=='KEGG'){
    keggs <-  msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")%>% dplyr::select(gs_name,gene_symbol)
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
  random_sample <- sample(full_gene_name,length(ES_genes_hsa_final_update))
  random_metric_per_pways <- ES_per_pathway(random_sample,database)
  random_metric_per_pways
}
get_var <- function(metric){
  library(SingleCellExperiment)
  library(scran)
  sce <- SingleCellExperiment(list(counts=metric))
  sce <- computeSumFactors(sce)
  sce <- logNormCounts(sce)
  dec <- modelGeneVar(sce)
  # plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
  # curve(metadata(dec)$trend(x), col="blue", add=TRUE)
  x <- dec$bio
  names(x) <- dec@rownames
  x
}
