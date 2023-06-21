# normal B lymphnode cell lines sequenced by six platforms 
# load("gene_counts.rds")
# samples <- gene_counts[c("10X_LLU_B_cellranger3.1",
#                          "10X_NCI_B_cellranger3.1",
#                          "10X_NCI_M_B_cellranger3.1",
#                          "C1_FDA_HT_B_featureCounts",
#                          "C1_LLU_B_featureCounts",
#                          "ICELL8_PE_B_featureCounts",
#                          "ICELL8_SE_B_featureCounts")]
library(biomaRt)
library(Seurat)
library(readxl)
library(stringr)
library(org.Hs.eg.db)
library(corrplot)
library(RColorBrewer)
# ES for human 
esg <- readxl::read_excel("deg_annotation_clean.xlsx")
esg.hsa <- esg[esg$organisms=="Homo sapiens",]
esg.hsa <- esg.hsa[!duplicated(esg.hsa$symbol),c("symbol","...9")]
esg.hsa <- esg.hsa[order(esg.hsa$symbol),]
esg.hsa <- esg.hsa[!duplicated(esg.hsa$...9),]
# Collapse with EXCEL gene names, manually changing it
esg.hsa[1:7,1] <- c("MARCHF5","SEPTIN5","MARCHF6",'SEPTIN6','MARCHF7','SEPTIN7','SEPTIN8')
# change to capital letters
esg.hsa$symbol <- toupper(esg.hsa$symbol)  
esg.hsa <- esg.hsa[!duplicated(esg.hsa$symbol),]
# convert from ENSEMBLE to gene symbol
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# data normalised 
sample_normed_list <- lapply(samples , sample_norm)
# genes that expressed across all platforms
# for plotting 
# fetch essesntial genes for human
ES_normed_list1 <- lapply(sample_normed_list,function(x){
  x <- x[rownames(x)%in%esg.hsa$symbol,]
  index <-  apply(GetAssayData(object = x, 
                               assay = "RNA", slot = "data"),1,function(y)sum(y>0))
  rownames(x)[index>3]})
ES_normed_share_number <- sapply(ES_normed_list1,length)
ES_normed_share_number <- ES_normed_share_number[order(ES_normed_share_number, decreasing = T)]
barplot( ES_normed_share_number, names.arg = c("ICELL8_SE","ICELL8_PE","10X_NCI_M","10X_NCI","10X_LLU",'C1_LLU',"C1_FDA_HT"), main = 'Number of essential genes across sequencing centers (human)')

essential_genes_hsa <-  Reduce(intersect,ES_normed_list1) 
ES_normed_shared_list <- lapply(sample_normed_list,function(x){x[rownames(x)%in%essential_genes_hsa,]}) # 5928 genes

# get the ranks
ES_rank_list <- lapply(sample_normed_list,function(s)get_ES_rank(essential_genes_hsa,s))

ES_rank_df <-  t(data.frame(matrix(unlist(ES_rank_list), nrow=length(ES_rank_list), byrow=TRUE)))
rownames(ES_rank_df) <- names(ES_rank_list[[1]])
colnames(ES_rank_df) <- stringi::stri_join(sapply(names(ES_rank_list),function(x)str_split(x,"_")[[1]][1]),sapply(names(ES_rank_list),function(x)str_split(x,"_")[[1]][2]),sep="_")
# visualise the ranks 
heatmap(t(ES_rank_df),scale ='row')
corrplot(cor(ES_rank_df),addCoef.col = 'black', tl.pos = 'd',tl.col = 'black',cl.pos = 'n',col = colorRampPalette(c("blue","orange"))(20))
 

# remove the "fluctuated genes"
ES_ranks_diff <- apply(ES_rank_df,1,function(x)sd(x))
hist(ES_ranks_diff,xlab = '', main='Standard deviation of the distribution of the ranks across sequencing methods (human)')

ES_ranks_names <- ES_ranks_diff[ES_ranks_diff<sd(ES_ranks_diff)*4]
# update the correlation list 
ES_ranks_final <- ES_rank_df[rownames(ES_rank_df)%in%names(ES_ranks_names),]
corrplot(cor(ES_ranks_final),addCoef.col = 'black', tl.pos = 'd',tl.col = 'black',cl.pos = 'n',col = colorRampPalette(c("blue","orange"))(20))

# whether it is significantly more expressed 
par(mfrow=c(2,4))
base::mapply(plot_mean_exp,ES_normed_list,sample_normed_sharegene_list)
# 
# plot_mean_exp(ES_normed_list[[7]],sample_normed_sharegene_list[[7]])+ title(colnames(ES_rank_df)[7])
write.csv(ES_ranks_names,"essential_genes_hsa_final.csv")
saveRDS(ES_ranks_names,"essential_genes_hsa_final.rds")
# functions ####
sample_norm <- function(counts){
  # consistency, use Seurat 
  # remove ERCC if applicable
  counts <- counts[!grepl("^gERCC", rownames(counts)),]
  genes <-  rownames(counts)
  gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                    values = genes, mart= mart)
  gene_IDs <- na.omit(gene_IDs)
  gene_IDs <- gene_IDs[gene_IDs$hgnc_symbol!='',]
  gene_IDs <- gene_IDs[!duplicated(gene_IDs$ensembl_gene_id),]
  gene_IDs <- gene_IDs[!duplicated(gene_IDs$hgnc_symbol),]
  counts <- counts[rownames(counts)%in%gene_IDs$ensembl_gene_id,]
  rownames(counts) <- gene_IDs[match(rownames(counts),gene_IDs$ensembl_gene_id ),2]
  
  counts <- CreateSeuratObject(counts = counts)
  counts[["percent.mt"]] <- PercentageFeatureSet(counts, pattern = "^MT-")
  counts <- SCTransform(counts, vars.to.regress = "percent.mt",verbose=TRUE)
  
}
get_ES_rank <- function(ES_list, full_metric){
  metric_rank <- rank(rowMeans(as.matrix(full_metric@assays$SCT@data)))
  ES_rank <- metric_rank[names(metric_rank)%in% ES_list]
  ES_rank <- ES_rank/dim(full_metric)[1]*100
  ES_rank
}
plot_mean_exp <- function(ES_metric,full_metric){
  set.seed(123441)
  random_sample <- sample(dim(full_metric)[1],dim(ES_metric)[1])
  random_metric <- subset(full_metric,features=rownames(full_metric)[random_sample])
  boxplot(rowMeans(ES_metric),rowMeans(as.matrix(random_metric@assays$SCT@data)),rowMeans(as.matrix(full_metric@assays$SCT@data)),notch = T,names=c('Essential genes',"random genes","All genes"))
  t.test(rowMeans(ES_metric),rowMeans(as.matrix(random_metric@assays$SCT@data)))
  
}
