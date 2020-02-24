#To find genes with similar expression patterns to the JUND cluster 
#(genes that covaried with median JUND cluster signature determined by Spearman correlation
#coeffcient > 0.5) in the 10-cell samples reprofiled by RNAseq

library(genefilter)
library(gplots)
library(EnvStats)
library(RColorBrewer)
library(psych)
library(pheatmap)
library(rstudioapi)

current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path))

rsemTPM_mcf10A_1 <- read.csv('/Users/ejp4gf/Google Drive/Janes Lab Drive/RNAseq/rsemTPM_mcf10A_1.csv',header=TRUE,sep=',')
rsemTPM_mcf10A_1_filt <- rsemTPM_mcf10A_1[(rowSums(rsemTPM_mcf10A_1[,20:24]>=5) >= 3),] #if genes are expressed at greater than 5 TPM in 3 or more pooled control samples, keep them

#getting rid of these two genes that had means across samples of zero (almost all values of samples were 0 but passed prefiltering of pooled samples)
rsemTPM_mcf10A_1_filt=rsemTPM_mcf10A_1_filt[-2079,] 
rsemTPM_mcf10A_1_filt=rsemTPM_mcf10A_1_filt[-3468,] 

#Reorganize TPM matrix
rsemTPM_mcf10A_1_filt[,10]=NULL #get rid of No RT sample
pooledsamples_TPM=rsemTPM_mcf10A_1_filt[,19:23]
rsemTPM_mcf10A_1_filt[,19:23]=NULL
means=rowMeans(rsemTPM_mcf10A_1_filt[,10:27])

#Row normalize TPM
rsem_TPM_rownorm=rsemTPM_mcf10A_1_filt
rsem_TPM_rownorm[,10:27]=rsemTPM_mcf10A_1_filt[,10:27]/means

#JUND cluster signature analysis---------
JUNDcluster= c('Coprs',' Slirp','Ndufa1','Sec61g','Nop10','Mrpl33','Tbca','Rps27l','Eef2','Rpl13','Tmem258','Rps21','Rpl38','Fam89b','Cdkn1a','Rps29','Rnh1','S100a6','Cox8a','Tceb2','Ilf2','Adirf','Api5','Atp5h','Krt5','Eif3m','Hspe1','Atp5e','Nme2','Marcks','Jund','Mus81','Rps6','Chchd1','Oaf','Aldh3b1','Keap1')
JUNDcluster=toupper(JUNDcluster)

indicesofJUND=which(rsem_TPM_rownorm[,3] %in% JUNDcluster)
JUNDcluster_rownorm=rsem_TPM_rownorm[c(indicesofJUND),]

#median JUND signature
median_JUNDsignature=apply(JUNDcluster_rownorm[,10:27],2,median)
median_JUNDsignature=data.frame(median_JUNDsignature)

genes_rownorm=cbind(name=rsem_TPM_rownorm$symbol,rsem_TPM_rownorm[,10:27]) #create matrix with just gene names and row norm expression values
data_mat=data.matrix(t(rsem_TPM_rownorm[,10:27])) #must transpose in order to use cor function later which computes correlation of columns of matrix
colnames(data_mat)=rsem_TPM_rownorm[,3]
df=data.frame(data_mat)

#Spearman analysis
correl_Spearman=corr.test(median_JUNDsignature,data_mat,method="spearman",adjust="BH") #same r values as using cor, using Benjamini-Hoch correction for p-values

correl_indices_spearman_10per=which(correl_Spearman$p<.1) #FDR cutoff 10%
genes_Spearman=rsem_TPM_rownorm[,3][correl_indices_spearman_10per] #566 genes out of 7228 TPM filtered

correl_Spearman=corr.test(median_JUNDsignature,data_mat, method='spearman')

correl_indices_Spearman=which(correl_Spearman$r>.5)

genes_Spearman=rsem_TPM_rownorm[,3][correl_indices_Spearman] #633 genes

expanded_cluster_Spearman=rsemTPM_mcf10A_1_filt[correl_indices_Spearman,]

#JUND cluster gene analysis 
JUND_Spearman=which(genes_Spearman %in% JUNDcluster)
genes_Spearman[JUND_Spearman] #25 genes of JUND cluster
correl_Spearman_JUND=correl_Spearman[which(rsem_TPM_rownorm[,3] %in% JUNDcluster)]

#this is just the matrix filtered based on expression in pooled controls, no normalization
expanded_Spearman_datamat=data.matrix(expanded_cluster_Spearman[,10:27]) 

#zscore
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

#Z scored un-rownormalized filtered matrix, then used pheatmap to cluster 
expanded_Spearman_datamat.norm <- t(apply(expanded_Spearman_datamat, 1, cal_z_score))

pheatmap(expanded_Spearman_datamat.norm, clustering_distance_rows = 'euclidean',clustering_distance_cols = 'euclidean',
         clustering_method="ward.D",labels_row=genes_Spearman,show_rownames = T, fontsize_row=2,col=rev(brewer.pal(11 ,"RdBu")),border_color=NA)