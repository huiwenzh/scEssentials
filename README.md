# scEssentials
A set of essential genes that is reliable to use in single-cell research, which is also implemented in our coexpression methods [scRoGG](https://github.com/huiwenzh/scRoGG/tree/master/R)  as reference geneset.
## Outline
We use large-scale benchmarking datasets as well as the mouse and human scRNA-seq atlases to refine robust essential genes and characterise their transcriptomic profiles across different cell types, denoted as scEssentials. scEssentials identified for both the mouse and human models are consistently high in expression and exhibit limited variably across more than 60 cell types. We also demonstrate a substantial number of significantly correlated gene pairs within scEssentials, which produce densely connected coexpression networks with functional annotation. Furthermore, we showed high frequencies of scEssentials across 200 pathways. Finally, we develop a score to quantify the relative essentiality of genes within scEssentials, which further validates with significant association with gene mutation frequency and chromatin accessibility. 

This study consists of three components:
  1. We measured the essential genes’ detectability and properties with respect to robustness, corresponding biological functions and coexpression connectivity in more than 100 cell types and 10 sequencing platforms to identify a set of essential genes explicitly for single-cell RNA-seq data (scEssentials) and define their essentiality score. These pre-ranked genes can be applied in any scRNA-seq dataset to either normalise the data or identify the gene expression changes.
  
  2. We evaluated the characteristics of the scEssentials through multi-omics data and investigated the relationship among different molecular levels. High expression levels were identified in all molecular levels and the expression variability was positively correlated with the essentiality score. 
  
  3. We demonstrated the stable expression for the scEssentials genes during ageing to showcase the characteristics of scEssentials in the dynamic and heterogeneous biological processes. 
  
## Data used in the analyses
  
Sequencing methods benchmarking datasets:
  
  1. Mouse: embryonic stem cells (mESCs) that have been sequenced by CELseq2, Dropseq, MARSeq, SCRBseq, Smartseq and Smartseq2. [Ziegenhain et al.](https://pubmed.ncbi.nlm.nih.gov/28212749/) 
  
  2. Human: ‘normal’ B lymphocyte line from breast cancer patients that have been sequenced by 10x Genomics Chromium, Fluidigm C1, Fluidigm C1 HT and Takara Bio ICELL8. [Chen et al.](https://pubmed.ncbi.nlm.nih.gov/33349700/)

Single-cell transcriptomic atlases:

  1. Tabula-Muris: Single-cell transcriptomics of 20 mouse organs. [The Tabula Muris Consortium., Overall coordination., Logistical coordination. et al.](https://www.nature.com/articles/s41586-018-0590-4)
  
  2. Tabula-Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans. [The Tabula Sapiens Consortium](https://pubmed.ncbi.nlm.nih.gov/35549404/)
  
Human gene damage index that measured gene mutation frequency:
  
  1. The human gene damage index as a gene-level approach to prioritizing exome variants. [Itan et al.](https://www.pnas.org/doi/10.1073/pnas.1518646112)

Single-cell mouse chromatin accessibility atlas:
  
  1.  A Single-Cell Atlas of In Vivo Mammalian Chromatin Accessibility. [Cusanovich et al.] (https://pubmed.ncbi.nlm.nih.gov/30078704/)
 
Single-cell mouse ageing atlas:

  1. A single-cell transcriptomic atlas characterizes ageing tissues in the mouse. [The Tabula Muris Consortium] (https://www.nature.com/articles/s41586-020-2496-1)
