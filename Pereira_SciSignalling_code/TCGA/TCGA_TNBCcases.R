#Process TCGA triple-negative breast cancer data 

#install.packages("gplots")
library(gplots)
library(readtext)
library(stringr)
library(data.table)
library(dplyr)
library(rstudioapi)

current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path))

#Load FPKM data from Xena and annotate with gene ID
FPKM<-read.table("TCGA-BRCA.htseq_fpkm.tsv",sep = '\t', header = TRUE,stringsAsFactors = F)
ensgtoid <- read.delim('ENSG_ID2Name.txt',sep='\t',header=F,stringsAsFactors = F)
row.names(FPKM) <- sub("\\..*","",FPKM$Ensembl_ID)
#row.names(FPKM)<-FPKM$Row.names
FPKM <- FPKM[,-1] #can remove bc they are row names now
row.names(ensgtoid) <- ensgtoid$V1
FPKM_gene<- merge(as.data.frame(FPKM), as.data.frame(ensgtoid), by='row.names')

FPKM_gene <- FPKM_gene[c(1219,1220,1:1218)] #CHANGE BASED ON YOUR INDICES
#tidy up for now
colnames(FPKM_gene)[1:2] <- c("ensmbl","symbol")
FPKM_gene <- FPKM_gene[,-3]
write.csv(FPKM_gene,'TCGA-BRCA.htseq-fpkm_geneid.csv',row.names = F) 
#COMMENT OUT ABOVE AND START FROM TCGA-BRCA.htseq-fpkm_geneid.csv IN FUTURE

#Get TNBC samples using clinical data (PAM50 subtype), use BRCA_clinicalmatrix that Sham sent me that has the subtype information
clinical_all<- read.delim('/BRCA_clinicalmatrix',header=TRUE, sep="\t",stringsAsFactors=FALSE)

#Different ways I tried to get triple-negative cases, first tried to pull out all by Basal PAM50
#but then realized not all were actually triple negative by ER,PR,HER2

# for (i in (3:1219)) {
#   patientID <- as.character(colnames(FPKM_gene[i]))
#   FPKM_gene[57289,i] <- t(as.data.frame(clinical_all[grep((str_sub(patientID,-11,-5)),(clinical_all[,1])),20]))
# }
# 
# FPKM_gene[57289,][FPKM_gene[57289,]==""] <- "NA"
# 
# FPKM_basal <- cbind(FPKM_gene[,c(1,which(FPKM_gene[57289,]== "Basal"))])
# 
# #MATCH BASED ON PATIENT IDs FROM TRIPLE NEGATIVE MATRIX KEVIN PUT ON SERVER 
# #(RSEM-normalized RNAseq data for 122 breast cancers classified as “triple-negative” by their pathology report)
# TNBCseq=read.csv(file="TCGA_TNBC_RNAseq.txt", header = TRUE, sep="\t") #Read in TGCA data
# patientID_TNBC<-as.character(colnames(TNBCseq[2:123]))
# 
# FPKM_TNBC<-FPKM_gene[,1:2]
# for (k in (3:1219)){
#   patientID <-as.character(colnames(FPKM_gene[k]))
#   if (str_sub(patientID,1,15) %in% patientID_TNBC) {
#     FPKM_TNBC[,k]<-FPKM_gene[,k]
# }
# }

FPKM_gene<-read.csv('TCGA-BRCA.htseq-fpkm_geneid.csv',header=T) 
patientID_FPKM<-as.character(colnames(FPKM_gene[3:1219]))
#Just talked to Sham, filtering based on these patient IDs from the breast cancers classified as "TN" by
#path does not always mean they get basal by PAM50, she said she contacted TCGA to ask and they said
#this is because PAM50 was made for microarray data and not sequencing data.....
#FPKM_TNBC <- FPKM_gene[,match(str_sub(patientID_FPKM,1,15),patientID_TNBC,nomatch=0)]

#select columns that are negative for ER, HER2, PR
triple_neg_samples<- clinical_all %>% filter(ER_Status_nature2012=="Negative" & HER2_Final_Status_nature2012=="Negative" & PR_Status_nature2012=="Negative")
#replace all dashes with periods in Sample IDs in clinical table to match IDs in FPKM table
triple_neg_samples$sampleID<-gsub("[-]",".",triple_neg_samples$sampleID)
#add 2 because the first two columns of FPKM_gene are ensemble ID and geneID
FPKM_TNBC_path<- FPKM_gene[,which(!is.na(match(str_sub(patientID_FPKM,1,15),triple_neg_samples$sampleID)),TRUE)+2]
FPKM_TNBC_path$ensembl<-FPKM_gene$ensmbl
FPKM_TNBC_path$symbol<-FPKM_gene$symbol
FPKM_TNBC_path <- FPKM_TNBC_path[c(123,124,1:122)] #FINAL MATRIX WITH ensembl ID, geneID, and TNBC FPKM

# START HERE IF INTERESTED IN TNBC cases!!!!!!!!!!!!! Write FPKM_TNBC_path b/c loading TCGA-BRCA.htseq-fpkm_geneid.csv takes a long time
#write.csv(FPKM_TNBC_path,"FPKM_TNBC_path.csv")
FPKM_TNBC_path<-read.csv("FPKM_TNBC_path.csv")

#TPM normalize
FPKM_TNBC_path_TPM<-(sweep(FPKM_TNBC_path[3:124],2,colSums(FPKM_TNBC_path[3:124]),`/`))*1e6
FPKM_TNBC_path_TPM$Gene<-FPKM_TNBC_path$symbol

genes<-c("KEAP1","NFE2L2","MAFF","MAFG","MAFK","ATM","CHEK2","TP53","MDM2","PPM1D","CDKN1A","TXN","SOD1","PRDX1","HMOX1")

TCGA_model<-FPKM_TNBC_path_TPM[match(genes, FPKM_TNBC_path_TPM$Gene),]
row.names(TCGA_model)<-genes

#Write patient IDs and info of TNBC tumors
patientIDs_TNBC<-colnames(FPKM_TNBC_path)[3:124]
write.csv(patientIDs_TNBC,"TCGA_TNBC_patientIDs_1-1-20.csv")
write.csv(triple_neg_samples,"TCGA_TNBC_clinicaldata_1-1-20.csv")

#Reorder clinical table so it's in same order as in FPKM_TNBC_path ("Tumor 1-122")
triple_neg_samples_ordered<-triple_neg_samples[match(str_sub(patientIDs_TNBC,1,15),triple_neg_samples$sampleID),]
write.csv(triple_neg_samples_ordered,"TCGA_TNBC_clinicaldata_1-1-20.csv")