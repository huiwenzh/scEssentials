# correlation of essential genes 
library(stringr)
library(data.table)
library(clusterProfiler)
library(Seurat)
library(scran)
library(SingleCellExperiment)
library(batchelor)
library(org.Mm.eg.db)
library(DescTools)
library(corrr)
library(msigdbr)
library(dplyr)
library(UpSetR)
esg <- readxl::read_excel("deg_annotation_clean.xlsx")
esg.hsa <- esg[esg$organisms=="Homo sapiens",]
esg.hsa <- unique(esg.hsa$...9)

esg.mmu <- esg[esg$organisms=="Mus musculus","symbol"]
esg.mmu <- sapply(esg.mmu,function(x)str_split(x,"/"))
esg.mmu <- unlist(esg.mmu)
esg.mmu <- unique(esg.mmu) #2513
#write.csv(esg.mmu,"essential_genes_mmu_final.csv")
esg.mmu <- bitr(esg.mmu, fromType = "SYMBOL",toType = c("ENTREZID",'ENSEMBL'), OrgDb = org.Mm.eg.db) #2354

MES <- read.table('R:/SNRNA2020-Q2288/scNetworks/across_platforms/GSE75790_ziegenhain_complete_data.txt',row.names = 1 )

pways_sum <- c(table(hallmarks$gs_name))
## Smartseq2##########
SmartSeq2 <- batch_removal(tech = 'SmartSeq2', data = MES)
SmartSeq2_ES <- get_ES_exp(as.matrix(SmartSeq2@assays$SCT@data))
barplot(rowMeans(SmartSeq2_ES)[order(rowMeans(SmartSeq2_ES), decreasing = T)])
plot_mesn_exp(SmartSeq2_ES,SmartSeq2)+title("SmartSeq2")
Smartseq2_rank <- get_ES_rank(rownames(SmartSeq2_ES),SmartSeq2)  
 
boxplot(Smartseq2_rank, names="Essential genes ranking", main="Smartseq2")
SmartSeq2_ES_cor <- ES_cor(SmartSeq2_ES, method = "pearson")
SmartSeq2_ES_cor1 <- ES_cor(SmartSeq2_ES, method = "Fisher")
 
boxplot(SmartSeq2_ES_cor$Freq,SmartSeq2_ES_cor1$Freq, names=c("Pearson","z-transformed pearson"), main="SmartSeq2")

SmartSeq2_ES_cor_rank <- ES_cor_rank(rownames(SmartSeq2_ES),full_metric = SmartSeq2 , method = "pearson")
SmartSeq2_ES_per_pways <- ES_per_pathway(rownames(SmartSeq2_ES),database = "Hallmark")
boxplot(sapply(SmartSeq2_ES_per_pways,length),main="Number of essential genes per pathway")
plot_pathway_random(rownames(SmartSeq2_ES),SmartSeq2,"Hallmark") 
pways_sum <- c(table(hallmarks$gs_name))
boxplot(sapply(SmartSeq2_ES_per_pways,length)/pways_sum*100,main="Percentage of essential genes per pathway")
## Smartseq ##########
SmartSeq  <- batch_removal(tech = 'SmartSeq[AB]', data = MES)
SmartSeq_ES <- get_ES_exp(as.matrix(SmartSeq@assays$SCT@data))
barplot(rowMeans(SmartSeq_ES)[order(rowMeans(SmartSeq_ES), decreasing = T)],main="SmartSeq")
plot_mesn_exp(SmartSeq_ES,SmartSeq )+title("SmartSeq")
Smartseq_rank <- get_ES_rank(rownames(SmartSeq_ES),SmartSeq)
boxplot(Smartseq_rank, names="Essential genes ranking", main="Smartseq")
SmartSeq_ES_cor <- ES_cor(SmartSeq_ES, method = "pearson")
SmartSeq_ES_cor1 <- ES_cor(SmartSeq_ES, method = "Fisher")

boxplot(SmartSeq_ES_cor$Freq,SmartSeq_ES_cor1$Freq, names=c("Pearson","z-transformed pearson"), main="SmartSeq")

#SmartSeq_ES_cor_rank <- ES_cor_rank(rownames(SmartSeq_ES),full_metric = SmartSeq  , method = "pearson")
SmartSeq_ES_per_pways <- ES_per_pathway(rownames(SmartSeq_ES),database = "Hallmark")
boxplot(sapply(SmartSeq_ES_per_pways,length),main="Number of essential genes per pathway")
plot_pathway_random(rownames(SmartSeq_ES),SmartSeq,"Hallmark")+title("Smartseq")
 
boxplot(sapply(SmartSeq_ES_per_pways,length)/pways_sum*100,main="Percentage of essential genes per pathway")
## CELseq2###########
CELseq2  <- batch_removal(tech = 'CELseq2', data = MES)
CELseq2_ES <- get_ES_exp(as.matrix(CELseq2@assays$SCT@data))
barplot(rowMeans(CELseq2_ES)[order(rowMeans(CELseq2_ES), decreasing = T)],main = "CELseq2")
plot_mesn_exp(CELseq2_ES,CELseq2 )+title("CELseq2")
CELseq2_rank <- get_ES_rank(rownames(CELseq2_ES),CELseq2)
boxplot(CELseq2_rank, names="Essential genes ranking", main="CELseq2")
CELseq2_ES_cor <- ES_cor(CELseq2_ES, method = "pearson")
CELseq2_ES_cor1 <- ES_cor(CELseq2_ES, method = "Fisher")

boxplot(CELseq2_ES_cor$Freq,CELseq2_ES_cor1$Freq, names=c("Pearson","z-transformed pearson"), main="CELseq2")
#CELseq2_ES_cor_rank <- ES_cor_rank(rownames(CELseq2_ES),full_metric = CELseq2 , method = "pearson")
CELseq2_ES_per_pways <- ES_per_pathway(rownames(CELseq2_ES),database = "Hallmark")
boxplot(sapply(CELseq2_ES_per_pways,length),main="Number of essential genes per pathway")
plot_pathway_random(rownames(CELseq2_ES),CELseq2,"Hallmark")+title("CELseq2")

boxplot(sapply(CELseq2_ES_per_pways,length)/pways_sum*100,main="Percentage of essential genes per pathway")
 

# DropSeq ###########
DropSeq  <- batch_removal(tech = 'DropSeq', data = MES)
DropSeq_ES <- get_ES_exp(as.matrix(DropSeq@assays$SCT@data))
barplot(rowMeans(DropSeq_ES)[order(rowMeans(DropSeq_ES), decreasing = T)],main = "DropSeq")
plot_mesn_exp(DropSeq_ES,DropSeq )+title("DropSeq")
DropSeq_rank <- get_ES_rank(rownames(DropSeq_ES),DropSeq)
boxplot(DropSeq_rank , names="Essential genes ranking", main="DropSeq")
DropSeq_ES_cor <- ES_cor(DropSeq_ES, method = "pearson")
DropSeq_ES_cor1 <- ES_cor(DropSeq_ES, method = "Fisher")

boxplot(DropSeq_ES_cor$Freq,DropSeq_ES_cor1$Freq, names=c("Pearson","z-transformed pearson"), main="DropSeq")
DropSeq_ES_cor_rank <- ES_cor_rank(rownames(DropSeq_ES),full_metric = DropSeq , method = "pearson")
DropSeq_ES_per_pways <- ES_per_pathway(rownames(DropSeq_ES),database = "Hallmark")
boxplot(sapply(DropSeq_ES_per_pways,length),main="Number of essential genes per pathway")
plot_pathway_random(rownames(DropSeq_ES),DropSeq,"Hallmark")+title("DropSeq")

boxplot(sapply(DropSeq_ES_per_pways,length)/pways_sum*100,main="Percentage of essential genes per pathway")

# MARSseq ###########
MARSseq  <- batch_removal(tech = 'MARSseq', data = MES)
MARSseq_ES <- get_ES_exp(as.matrix(MARSseq@assays$SCT@data))
barplot(rowMeans(MARSseq_ES)[order(rowMeans(MARSseq_ES), decreasing = T)],main = "MARSseq")
plot_mesn_exp(MARSseq_ES,MARSseq )+title("MARSseq")
MARSseq_rank <- get_ES_rank(rownames(MARSseq_ES),MARSseq)
boxplot(MARSseq_rank , names="Essential genes ranking", main="MARSseq")
MARSseq_ES_cor <- ES_cor(MARSseq_ES, method = "pearson")
MARSseq_ES_cor1 <- ES_cor(MARSseq_ES, method = "Fisher")

boxplot(MARSseq_ES_cor$Freq,MARSseq_ES_cor1$Freq, names=c("Pearson","z-transformed pearson"), main="MARSseq")
MARSseq_ES_cor_rank <- ES_cor_rank(rownames(MARSseq_ES),full_metric = MARSseq , method = "pearson")
MARSseq_ES_per_pways <- ES_per_pathway(rownames(MARSseq_ES),database = "Hallmark")
boxplot(sapply(MARSseq_ES_per_pways,length),main="Number of essential genes per pathway")
plot_pathway_random(rownames(MARSseq_ES),MARSseq,"Hallmark")+title("MARSseq")

boxplot(sapply(MARSseq_ES_per_pways,length)/pways_sum*100,main="Percentage of essential genes per pathway")


# SCRBseq ###########
SCRBseq  <- batch_removal(tech = 'SCRBseq', data = MES)
SCRBseq_ES <- get_ES_exp(as.matrix(SCRBseq@assays$SCT@data))
barplot(rowMeans(SCRBseq)[order(rowMeans(SCRBseq), decreasing = T)],main = "SCRBseq")
plot_mesn_exp(SCRBseq_ES,SCRBseq )+title("SCRBseq")
SCRBseq_rank <- get_ES_rank(rownames(SCRBseq_ES),SCRBseq)
boxplot(SCRBseq_rank,names="Essential genes ranking", main="SCRBseq")
SCRBseq_ES_cor <- ES_cor(SCRBseq_ES, method = "pearson")
SCRBseq_ES_cor1 <- ES_cor(SCRBseq_ES, method = "Fisher")

boxplot(SCRBseq_ES_cor$Freq,SCRBseq_ES_cor1$Freq, names=c("Pearson","z-transformed pearson"), main="SCRBseq")
#SCRBseq_ES_cor_rank <- ES_cor_rank(rownames(SCRBseq_ES),full_metric = SCRBseq , method = "pearson")
SCRBseq_ES_per_pways <- ES_per_pathway(rownames(SCRBseq_ES),database = "Hallmark")
boxplot(sapply(SCRBseq_ES_per_pways,length),main="Number of essential genes per pathway")
plot_pathway_random(rownames(SCRBseq_ES),SCRBseq,"Hallmark")+title("SCRBseq")

boxplot(sapply(SCRBseq_ES_per_pways,length)/pways_sum*100,main="Percentage of essential genes per pathway")


# compare and summary######
ES_list <- list(rownames(SmartSeq2_ES),rownames(SmartSeq_ES),
                rownames(CELseq2_ES),rownames(DropSeq_ES),
                rownames(MARSseq_ES),rownames(SCRBseq_ES))
names(ES_list) <- c('SmartSeq2','SmartSeq',"CELseq2",'DropSeq','MARSseq',"SCRBseq")
ES_shared_list <- Reduce( intersect,ES_list)   #1265
upset(fromList(ES_list),nsets = 6, order.by = 'freq')

#ES_ranks heatmap
ES_ranks <- cbind(Smartseq2_rank[names(Smartseq2_rank)%in%ES_shared_list],Smartseq_rank[names(Smartseq_rank)%in%ES_shared_list],CELseq2_rank[names(CELseq2_rank)%in%ES_shared_list],DropSeq_rank[names(DropSeq_rank)%in%ES_shared_list],MARSseq_rank[names(MARSseq_rank)%in%ES_shared_list],SCRBseq_rank[names(SCRBseq_rank)%in%ES_shared_list] )
colnames(ES_ranks) <-  c('SmartSeq2','SmartSeq',"CELseq2",'DropSeq','MARSseq',"SCRBseq")
heatmap(t(ES_ranks),scale ='row')
library(corrplot)
library(RColorBrewer)
corrplot(cor(ES_ranks),addCoef.col = 'black', tl.pos = 'd',tl.col = 'black',cl.pos = 'n',col = colorRampPalette(c("blue","orange"))(20))

barplot(c(dim(SmartSeq2_ES)[1],dim(SmartSeq_ES)[1],dim(CELseq2_ES)[1],dim(DropSeq_ES)[1],dim(MARSseq_ES)[1],dim(SCRBseq_ES)[1]),names.arg = c('SmartSeq2','SmartSeq',"CELseq2",'DropSeq','MARSseq',"SCRBseq"),main = 'Number of essential genes across sequencing method (mouse)')

# remove the "fluctuated genes"
ES_ranks_diff <- apply(ES_ranks,1,function(x)sd(x))
hist(ES_ranks_diff,xlab = '', main='Standard deviation of the distribution of the ranks across sequencing methods (mouse)')
  
ES_ranks_names <- ES_ranks_diff[ES_ranks_diff<sd(ES_ranks_diff)*4]
# update the correlation list 
ES_ranks_final <- ES_ranks[rownames(ES_ranks)%in%names(ES_ranks_names),]
corrplot(cor(ES_ranks_final),addCoef.col = 'black', tl.pos = 'd',tl.col = 'black',cl.pos = 'n',col = colorRampPalette(c("blue","orange"))(20))

# Figures to generate together 
# boxplot to show higher expression than random 
par(mfrow=c(2,3))
plot_mesn_exp(SmartSeq2_ES,SmartSeq2)+title("SmartSeq2")
plot_mesn_exp(SmartSeq_ES,SmartSeq)+title("SmartSeq")
plot_mesn_exp(CELseq2_ES,CELseq2 )+title("CELseq2")
plot_mesn_exp(DropSeq_ES,DropSeq )+title("DropSeq")
plot_mesn_exp(MARSseq_ES,MARSseq )+title("MARSseq")
plot_mesn_exp(SCRBseq_ES,SCRBseq )+title("SCRBseq")

# boxplot to show more involvment 
plot_pathway_random(rownames(MARSseq_ES),MARSseq,"Hallmark")

# save the excel file
write.csv(ES_ranks_names,"essential_genes_mmu_final.csv")
saveRDS(ES_ranks_names,"essential_genes_mmu_final.rds")
### functions ###########
batch_removal <- function(tech, data){
  pattern <- paste0('^',tech)
  data_sub<- data[,grep(pattern,colnames(data))] # both A + B
  # remove ERCC if applicable
  data_sub <- data_sub[!grepl("^gERCC", rownames(data_sub)),]
  
  data_sub <- CreateSeuratObject(counts = data_sub)
  data_sub$batch <- unlist(sapply(data_sub$orig.ident,function(x)str_split(x,"_")))
  
  data_sub <- SCTransform(data_sub, vars.to.regress = "batch",verbose=TRUE)
  data_sub
}

get_ES_exp <- function(metric, type="ENSEMBL"){
  # @metric - normalized and transformed data used for comparison
  # @type - gene symbol type
  #m1 <- metric[rownames(metric)%in%names(ES_ranks_names),]
  index <-  apply(metric,1,function(y)sum(y>0))
  metric <- metric[index>3,]
  m1 <- metric[rownames(metric)%in%esg.mmu[,type],]
  # rownames(m1) <- na.omit(esg.mmu[match(rownames(m1) ,esg.mmu$ENSEMBL),"SYMBOL"])
  m1
}

get_ES_rank <- function(ES_list, full_metric){
  metric_rank <- rank(rowMeans(as.matrix(full_metric@assays$SCT@data)))
  #ES_name <- na.omit(esg.mmu[match(ES_list,esg.mmu$SYMBOL),"ENSEMBL"])
  ES_rank <- metric_rank[names(metric_rank)%in%ES_list]
  ES_rank <- ES_rank/dim(full_metric)[1]*100
}

get_sim_ranks <- function(metric){
  metric_rank <- rank(rowMeans( metric))
  
  all_genes <- rownames(metric)
  set.seed(123456789)
  random_genes <- sample(all_genes[all_genes%notin%ES_shared_list],1000)

  Sim_rank <- metric_rank[names(metric_rank)%in%random_genes]
  Sim_rank <- Sim_rank/dim(metric)[1]*100
  
}

plot_mesn_exp <- function(ES_metric,full_metric){
  set.seed(1234411)
  random_sample <- sample(dim(full_metric)[1],dim(ES_metric)[1])
  random_metric <- subset(full_metric,features=rownames(full_metric)[random_sample])
  boxplot(rowMeans(ES_metric),rowMeans(as.matrix(random_metric@assays$SCT@data)),rowMeans(as.matrix(full_metric@assays$SCT@data)),notch = T,names=c('Essential genes',"random genes","All genes"))
  t.test(rowMeans(ES_metric),rowMeans(as.matrix(random_metric@assays$SCT@data)))
  
}

ES_cor <- function(ES_metric, method ){
  if(method == "pearson" ){
    cor_met <- as.data.frame(as.table(cor(t(ES_metric),method = method)))
    cor_met <- cor_met[cor_met$Freq!=1,]    
  }
  if (method == "Fisher"){
    cor_met <- as.data.frame(as.table(cor(t(ES_metric),method = "pearson")))
    cor_met <- cor_met[cor_met$Freq!=1,]    
    cor_met$Freq <- FisherZ(cor_met$Freq)
  }
  cor_met
  }

ES_cor_rank <- function(ES_list, full_metric, method) {
  if (method == "pearson"){
    cor_met <-  cor(t(as.matrix(full_metric@assays$SCT@data)),method = method) 
    ES_rank <- cor_met[rownames(cor_met)%in%ES_list,]
    ES_rank <- as.data.frame(as.table(cor_met))
    ES_rank <- ES_rank[ES_rank$Freq!=1,]  
    ES_rank$rank <- rank(ES_rank$Freq)
 
  }
  if (method == "Fisher"){
    cor_met <-  cor(t(as.matrix(full_metric@assays$SCT@data)))
    cor_met <- cor_met[cor_met$Freq!=1,]    
    cor_met$Freq <- FisherZ(cor_met$Freq)
    cor_met$rank <- rank(cor_met$Freq)
 
    ES_rank <- cor_met[cor_met$Var1%in%ES_list & cor_met$Var2%in%ES_list , ]
    
  }

  ES_rank
}
 
ES_per_pathway <- function(ES_list, database){
  if (database=="Hallmark"){
    hallmarks <-  msigdbr(species = "Mus musculus", category = "H")%>% dplyr::select(gs_name,ensembl_gene)
    pways_list <- list()
    for (pways in unique(hallmarks$gs_name)){
      pways_list[[pways]] <- hallmarks[hallmarks$gs_name== pways,2]
    }
    overlap <- lapply(pways_list, function(x){x[["ensembl_gene"]][x[["ensembl_gene"]]%in%ES_list]})
    overlap
    }
}

plot_pathway_random <- function(ES_list,full_metric,database){
  set.seed(1234411)
  random_sample <- sample(dim(full_metric)[1],length(ES_list))
  random_metric <- subset(full_metric,features=rownames(full_metric)[random_sample])
  random_metric_per_pways <- ES_per_pathway(rownames(random_metric),"Hallmark")
  ES_per_pways <- ES_per_pathway(ES_list,database = "Hallmark")
 
  boxplot(sapply(ES_per_pways,length),sapply(random_metric_per_pways,length), notch = T,names=c('Essential genes',"random genes" ))
  
}


# SmartSeq2.sce <-  SingleCellExperiment(assays = list(counts = SmartSeq2))
# altExp(SmartSeq2.sce, "spike-in") <- SummarizedExperiment(ERCC)
# 
# SmartSeq2.sce <- computeSumFactors(SmartSeq2.sce)
# SmartSeq2.sce <- logNormCounts(SmartSeq2.sce)
# SmartSeq2.sce$batch <- c(rep("A",80),rep("B",77))
# 
# set.seed(101)
# SmartSeq2.sce.out <- fastMNN(SmartSeq2.sce, batch=SmartSeq2.sce$batch)
# str(reducedDim(f.out2, "corrected"))

`%notin%` <- Negate(`%in%`)
