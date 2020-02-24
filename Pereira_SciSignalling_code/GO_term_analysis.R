#GO term analysis
#Generated GO terms for expanded JUND cluster (based on Spearman correlation >0.5 with median JUND signature)
#Analyzing the NEW GO terms that came up from analysis of expanded cluster- refine this list to get rid of 
#children GO terms from original microarray GO term analysis and get rid of GO terms with <35 genes annotated 
#to them

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GO.db", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GOfuncR", version = "3.8")

#Install Homo-sapiens database to get genes within GO terms 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Homo.sapiens", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GeneAnswers", version = "3.8")

library(GO.db)
library(readxl)
library(GOfuncR)
library(prob)
library(Homo.sapiens)
library(dplyr)
library(GeneAnswers)
library(convertToGeneAnswers)

setwd("/Users/ejp4gf/Google Drive/Janes Lab Drive/Bioinformatics")

#Read in GO terms generated from microarray cluster
GOannot_JUND<-readxl::read_excel("GOterms_microarray_vs_exp_intersectionanalysis.xlsx", sheet = "microarrayGOterms",col_names="GOterm")
GOannot_JUND_str<-apply(GOannot_JUND,1,toString) #convert to characters so I can use strsplit

#Read in the new GO terms generated from expanded cluster
GOannot_expanded<-readxl::read_excel("GOterms_microarray_vs_exp_intersectionanalysis.xlsx", sheet = "newGOterms")
GOannot_expanded_str<-apply(GOannot_expanded[,1],1,toString) #convert to characters so I can use strsplit

#Just get GOID of each GO annotation
GOIDs_JUND<-unlist(lapply(GOannot_JUND_str,function(x){unlist(strsplit(x,'[()]'))[2]}))
  
GOIDs_expanded<-unlist(lapply(GOannot_expanded_str,function(x){unlist(strsplit(x,'[()]'))[2]}))
GOIDs_expanded_FDR<-data.frame(GOIDs_expanded,GOannot_expanded$FDR_adj_pval)

#Get parent nodes of expanded cluster new GO terms
parents<-get_parent_nodes(GOIDs_expanded)

#Get children nodes of microarray cluster GO terms (not sure which comparison will be easier)
children<-get_child_nodes(GOIDs_JUND)

#Find intersection between children of JUND cluster GO terms and new GO terms from expanded cluster
int<-intersect(children$child_go_id,GOIDs_expanded)

#Remove GO terms that are children of terms from JUND cluster (106 terms remaining of 172)
GOIDs_expanded_refined<-GOIDs_expanded[-which(GOIDs_expanded %in% children$child_go_id)]

#Get genes that are in each remaining GO term
genes<-get_anno_genes(GOIDs_expanded_refined, database = 'Homo.sapiens')

genesperGO<-genes %>% group_by(go_id) %>% summarise(genes=paste(gene,collapse=", "),num=length(gene))

#Filter out GO terms with <35 genes annotated to them (101 terms remaining of 101)
lessthan35<-genesperGO[c(which(genesperGO$num<=35)),]
GOIDs_expanded_refined<-GOIDs_expanded_refined[-which(GOIDs_expanded_refined %in% lessthan35$go_id)]

GOIDs_expanded_refined_df<-data.frame(GOIDs_expanded_refined)
GOIDs_expanded_refined<-get_names(GOIDs_expanded_refined)[,1:2]

#Append FDR-adjusted pvalues to refined GO term database
FDRs<-GOIDs_expanded_FDR$GOannot_expanded.FDR_adj_pval
GOIDs_expanded_refined$FDRpval<-FDRs[which(as.character(GOIDs_expanded_FDR$GOIDs_expanded) %in% GOIDs_expanded_refined$go_id)]


#Get rid of one UNCLASSIFIED GO term, left with 100 GO terms
GOIDs_expanded_refined<-GOIDs_expanded_refined[!is.na(GOIDs_expanded_refined[,2]),]


