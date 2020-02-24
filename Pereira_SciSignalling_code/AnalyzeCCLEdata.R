#CCLE 
#https://ocg.cancer.gov/ctd2-data-project/translational-genomics-research-institute-tgen-quantified
#Translational Genomics Research Institute (TGen): Quantified Cancer Cell Line Encyclopedia (CCLE) RNA-seq Data
#Many applications analyze quantified transcript-level abundances to make inferences.  
#Having completed this computation across the large sample set, the CTD2 Center at the Translational Genomics Research 
#Institute presents the quantified data in a straightforward, consolidated form for these types of analyses.

#Experimental Approaches-
#After downloading RNA-seq data for 935 cell lines from the Cancer Cell Line Encyclopedia (link is external) (CCLE), 
#transcript-level abundance was quantified using Salmon1. All data were aligned using Salmon 0.4.2 using Homo Sapiens 
#GRCh37.74 (link is external) for reference. Raw BAM files used to generate this data is avaliable at GDC. 
#The resulting 935 quantification files, named by sample ID, have 4 columns for ensemble gene ID, length, number of reads,
#and transcripts per million (TPM). Other Salmon arguments were "--libType IU" (inward, unstranded).

#Download CCLE_id_mapping.txt, ccle_gene_quantification.zip from link above, also need ENSG_ID2Name.txt file

setwd("/Users/ejp4gf/Google Drive/Janes Lab Drive/CCLE")

#I only want to analyze breast cancer cell lines so use CCLE_id_mapping.txt file to identify 
cellines<-read.delim('CCLE_id_mapping.txt',sep='',header=T,stringsAsFactors = F)
breast_cellines<-cellines[grep("BREAST",cellines$ccl_name_primary_site),]

ensgtoid <- read.delim('ENSG_ID2Name.txt',sep='\t',header=F,stringsAsFactors = F)

row.names(ena_data) <- ena_data$Gene_ID
ena_data <- ena_data[,-1] #can remove bc they are row names now
row.names(ensgtoid) <- ensgtoid$V1
ena_data <- merge(as.data.frame(ena_data), as.data.frame(ensgtoid), by='row.names')


#9-26-19 analysis- using rsem_genes_tpm data from https://portals.broadinstitute.org/ccle/data
setwd("/Users/ejp4gf/Google Drive/Janes Lab Drive/CCLE")
CCLE_TPM<-read.delim('CCLE_RNAseq_rsem_genes_tpm_20180929.txt',header=TRUE,sep="\t")
rownames(CCLE_TPM) <- gsub("\\..*","",CCLE_TPM$gene_id)
CCLE_TPM <- CCLE_TPM[,-c(1,2)] #can remove bc they are row names now

ensgtoid <- read.delim('ENSG_ID2Name.txt',sep='\t',header=F,stringsAsFactors = F)
row.names(ensgtoid) <- ensgtoid$V1
CCLE_TPM <- merge(as.data.frame(CCLE_TPM), as.data.frame(ensgtoid), by='row.names')

CCLE_TPM <- CCLE_TPM[c(1021,1022,1:1020)] #CHANGE BASED ON YOUR INDICES
#tidy up for now
colnames(CCLE_TPM)[1:2] <- c("ensmbl","symbol")
CCLE_TPM <- CCLE_TPM[,-3]

breast_cellines<-CCLE_TPM[,grep("BREAST",colnames(CCLE_TPM))]
breast_CCLE_TPM<-cbind(CCLE_TPM$ensmbl,CCLE_TPM$symbol,breast_cellines)

TNBC_rnaseq=breast_CCLE_TPM[,c(1,2,4,7,10,13,16,17,18,23,24,25,28,30,31,35,40,42,45,47)]