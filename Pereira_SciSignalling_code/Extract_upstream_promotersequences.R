#From fasta file of promoter sequences (downloaded from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/),
#extract promoter sequences of expanded JUND cluster genes (all of the transcript variants)
#Also extract 250bp after TSS into each gene (from hg19.2bit downloaded from genome.ucsc)

install.packages("seqinr")
library("seqinr")

setwd("/Users/ejp4gf/Google Drive/Janes Lab Drive/Bioinformatics")
#upstreamseq_5000=read.fasta(file="upstream5000.fa",as.string=TRUE)
upstreamseq_1000=read.fasta(file="upstream1000.fa",as.string=TRUE)

#Get extended JUND cluster list, choose smallest (oldest) accession number
JUNDextendedcluster_NM=read.delim(file="JUNDexpandedcluster_NM.txt",header=TRUE)
colnames(JUNDextendedcluster_NM)=c("GeneName","NM_")
JUNDextendedcluster_NM$GeneName=as.factor(JUNDextendedcluster_NM$GeneName)

#Before aggregating by gene name, get rid of NM_, then aggregate and compare digits --> use smallest NM number (oldest RefSeq)
JUNDextendedcluster_NM$NMnum=gsub("^.{0,3}","",JUNDextendedcluster_NM$NM_)
agg_list=aggregate(. ~JUNDextendedcluster_NM$GeneName, data=JUNDextendedcluster_NM,toString)
agg_list$GeneName=NULL
agg_list$NM_=NULL
colnames(agg_list)=c("GeneName","NM_num")

#Find minimum
agg_list$min=sapply(1:nrow(agg_list),function(x){min(as.numeric(strsplit(agg_list$NM_num,',')[[x]]))})

#pad with front zeros to get back to correct ID (converting to numeric in previous step removed leading zeros)
agg_list$min_zero=ifelse(nchar(agg_list$min)>6,sprintf("%09d",agg_list$min),sprintf("%06d",agg_list$min))

#add NM_ prefix back  
agg_list$NM_min=paste0("NM_",agg_list$min_zero)   

#create dataframe with just gene name and lowest/oldest NM
JUND_NM=data.frame(agg_list$GeneName,agg_list$NM_min)
colnames(JUND_NM)=c("GeneName","NM_min")

#Get promoter sequence for all accession numbers of JUND cluster

pat=("^((?:[^_]*_){2}).*")
just_NM=sub(pat,'\\1',names(upstreamseq_1000))
just_NM=gsub(".{1}$","",just_NM)

JUNDcluster_fasta=upstreamseq_1000[which(just_NM %in% JUND_NM$NM_min)]

JUNDcluster_NM=sub(pat,"\\1",names(JUNDcluster_fasta))
JUNDcluster_NM=gsub(".{1}$","",JUNDcluster_NM)
JUNDcluster_promoters=getSequence(JUNDcluster_fasta,as.string=TRUE)
JUNDcluster_promoters_reform=lapply(JUNDcluster_promoters,unlist) 


JUNDextendedcluster_NMstring=paste(JUND_NM$NM_min)
gene_names=paste(JUND_NM$GeneName[match(JUNDcluster_NM, JUNDextendedcluster_NMstring)])

JUNDcluster_namesprom=cbind(gene_names,JUNDcluster_NM,JUNDcluster_promoters_reform)
JUNDcluster_df=data.frame(JUNDcluster_namesprom) #dataframe containing gene names, NM RefSeq ID and promoter sequence
colnames(JUNDcluster_df)=c("Gene symbol","RefSeq NM","Promoter sequence 1000 bases")

JUNDcluster_df_char=apply(JUNDcluster_df,2,as.character)

write.table(JUNDcluster_df_char,file="JUNDextendedcluster_promoterseq1000.txt",sep="\t",col.name=TRUE)

write.fasta(sequences=JUNDcluster_promoters_reform,names=gene_names,file.out="JUNDextendedcluster_promoterseq5000.fasta",open="w") #create file using gene names 

write.fasta(sequences=JUNDcluster_promoters_reform,names=JUNDcluster_NM,file.out="JUNDextendedcluster_promoterseq5000_NM.fasta",open="w") #create file using NM identifiers

#Get rid of entries with duplicated NM numbers (for entry into MEME all have to have unique identifiers)
# dup=duplicated(JUNDcluster_NM)
# rows_to_remove=which(dup==TRUE)
# 
# JUNDcluster_promoters_reform_dupsremoved=JUNDcluster_promoters_reform[-c(rows_to_remove)]
# JUNDcluster_NM_dupsremoved=JUNDcluster_NM[-c(rows_to_remove)]
# 
# write.fasta(sequences=JUNDcluster_promoters_reform_dupsremoved,names=JUNDcluster_NM_dupsremoved,file.out="JUNDextendedcluster_promoterseq5000_NM_dupsremoved.fasta",open="w") #create file using NM identifiers
# 
