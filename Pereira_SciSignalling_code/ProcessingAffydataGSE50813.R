#This code is for processing and analyzing microarray data from GSE50813, but at the end
#has code to display data in different useful ways (plotting information for each gene using
#small multiples, among other ways)

#Change directory to where files from GSE are downloaded
setwd("/Users/ejp4gf/Google\ Drive/Janes\ Lab\ Drive/SV40TAg\ mouse\ model\ of\ TNBC")
#setwd("/Users/elizabethpereira/Google\ Drive/Janes\ Lab\ Drive/SV40TAg\ mouse\ model\ of\ TNBC")
source("http://bioconductor.org/biocLite.R")
biocLite('affy')
library(affy)
biocLite("limma")
library(limma)
biocLite("annotate")
library(annotate)
install.packages("ggplot2")
library(ggplot2)
library(reshape2)
library(dplyr)
install.packages("dplyr")
library(dplyr)


affy.data=ReadAffy() #Read raw affymetrix data, creates an AffyBatch object
nvals=rma(affy.data) #RMA normalize to get values in GSE50813 family matrix

ned=exprs(nvals)   #Extract expression valus from Affy Object 
nsamp=sampleNames(nvals)  #Extract sample names
nprobes=featureNames(nvals) #Extract probeIDs

annot=read.csv(file="GPL1261-tbl-1.txt", header = FALSE, sep="\t") #Read in tab-delimited annotation file with ProbeIDs, GeneNames, and EntrezIDs

info=annot[,c(1,11,12)]
colnames(info)=c("Probe ID", "Gene Name", "Entrez ID")
expwithGeneNames=data.frame(info[,2],ned)
colnames(expwithGeneNames)[1]="GeneName"

#Write processed data matrix for Kevin 12-9-19
write.csv(expwithGeneNames,"C3TagProcessedMicroarray.csv")

#JUND cluster analysis-----------

JUNDcluster= c('Nfe2l2','Coprs', 'Nqo1',' Slirp','Ndufa1','Sec61g','Nop10','Mrpl33','Tbca','Rps27l','Eef2','Rpl13','Tmem258','Rps21','Rpl38','Fam89b','Cdkn1a','Rps29','Rnh1','S100a6','Cox8a','Tceb2','Ilf2','Adirf','Api5','Atp5h','Krt5','Eif3m','Hspe1','Atp5e','Nme2','Marcks','Jund','Mus81','Rps6','Chchd1','Oaf','Aldh3b1', 'Trp53','Keap1','Kras')

#Put expression data of JUND cluster genes in a matrix
JUNDclustermat<-NULL
for(i in 1:length(JUNDcluster)){
	JUNDclustermat=rbind(JUNDclustermat,expwithGeneNames[which(expwithGeneNames$GeneName ==JUNDcluster[i]),])
	}
#add genes who have other identifiers in GeneName column
JUNDclustermat=rbind(JUNDclustermat,expwithGeneNames[which(grepl('Sec61g',expwithGeneNames$GeneName)),])
JUNDclustermat=rbind(JUNDclustermat,expwithGeneNames[which(grepl('Gm10071 /// Rpl13',expwithGeneNames$GeneName)),])
JUNDclustermat=rbind(JUNDclustermat,expwithGeneNames[which(grepl('Rps21',expwithGeneNames$GeneName)),])
JUNDclustermat=rbind(JUNDclustermat,expwithGeneNames[which(grepl('Fam89b',expwithGeneNames$GeneName)),])
JUNDclustermat=rbind(JUNDclustermat,expwithGeneNames[which(grepl('Adirf',expwithGeneNames$GeneName)),])
	
GeneName=JUNDclustermat[,1]
write.table(JUNDclustermat, file='JUNDclusterexp.txt', sep="\t")

#If multiple probes corresponded to the same gene, then the expression values of these probe sets were averaged
JUNDclusteravg <- aggregate(. ~ GeneName, data = JUNDclustermat, mean)
write.table(JUNDclusteravg, file='JUNDclusteravgexp.txt', sep="\t")

#Data processing
JUNDclusteravg_df=data.frame(JUNDclusteravg[,2:25])
M=matrix_means_num=mapply(JUNDclusteravg,FUN=as.numeric)

n2=c("Control_8week","Control_8week","Control_8week","Control_8week","Control_8week", "Tumor_8week","Tumor_8week","Tumor_8week","Tumor_8week","Tumor_8week","Tumor_12week","Tumor_12week","Tumor_12week","Tumor_12week","Tumor_12week", "Tumor_16week","Tumor_16week","Tumor_16week","Tumor_16week","Tumor_16week","Tumor_20week","Tumor_20week","Tumor_20week","Tumor_20week")
colnames(JUNDclusteravg_df)=n2

control_med=NULL
for(i in 1:39){
	control_med[i]=median(M[i,2:6])
}

control_mean=NULL
JUNDclusteravg_df_control=NULL
for(i in 1:39){
	control_mean[i]=mean(M[i,2:6])
	JUNDclusteravg_df_control=JUNDclusteravg_df[,1:5]
}

week8_tumor_mean=NULL
JUNDclusteravg_df_week8=NULL
for(i in 1:39){
	week8_tumor_mean[i]=mean(M[i,7:11])
	JUNDclusteravg_df_week8=JUNDclusteravg_df[,6:10]
}

week12_tumor_mean=NULL
JUNDclusteravg_df_week12=NULL
for(i in 1:39){
	week12_tumor_mean[i]=mean(M[i,12:16])
	JUNDclusteravg_df_week12=JUNDclusteravg_df[,11:15]
}

week16_tumor_mean=NULL
JUNDclusteravg_df_week16=NULL
for(i in 1:39){
	week16_tumor_mean[i]=mean(M[i,17:21])
	JUNDclusteravg_df_week16=JUNDclusteravg_df[,16:20]
}

week20_tumor_mean=NULL
JUNDclusteravg_df_week20=NULL
for(i in 1:39){
	week20_tumor_mean[i]=mean(M[i,22:25])
	JUNDclusteravg_df_week20=JUNDclusteravg_df[,21:24]
}

matrix_means=cbind(control_mean,week8_tumor_mean,week12_tumor_mean,week16_tumor_mean,week20_tumor_mean)
matrix_means_num=as.matrix(matrix_means)
matrix_means_num=mapply(matrix_means_num,FUN=as.numeric)
matrix_means_num=matrix(data=matrix_means_num,ncol=5,nrow=39)
centered_matrix_means=sweep(matrix_means_num,1,control_med,"-")  #Center around median relative abundance in 8 week control animals 
centered_matrix_means_df=data.frame(centered_matrix_means)
centered_matrix_means_df_wgenes=data.frame(cbind(JUNDclusteravg[,1],data.frame(centered_matrix_means)))
n=c("Control_8week", "Tumor_8week","Tumor_12week", "Tumor_16week","Tumor_20week")
colnames(centered_matrix_means_df)=n

#Plot all centered means
df=melt(centered_matrix_means_df)
df$rowid=JUNDclusteravg[,1]
pdf(file="JUNDclustermeanstogether.pdf")
ggplot(df,aes(variable,value,group=factor(rowid)))+geom_point(aes(color=factor(rowid)))+theme(text=element_text(size=8),axis.text.x=element_text(angle=45,hjust=1))
dev.off()

#Plot individually
pdf(file="JUNDclustermeans.pdf")
ggplot(df,aes(variable,value,group=factor(rowid)))+geom_point(aes(color=factor(variable)),shape=15,size=2)+facet_wrap(~ rowid)+theme(text=element_text(size=8),axis.text.x=element_text(angle=45,hjust=1),legend.position="none")
dev.off()

#Plot all data points
JUNDclustermatrix=as.matrix(JUNDclusteravg_df)
JUNDclustercentered_matrix=sweep(JUNDclustermatrix,1,control_med,"-")  #Center around median relative abundance in 8 week control animals 
JUNDclustercentered_matrix_df=data.frame(JUNDclustercentered_matrix)
df2=melt(JUNDclustercentered_matrix_df)
df2$rowid=JUNDclusteravg[,1]

n3=c("Control_8week"=1:195,"Tumor_8week"=196:390,"Tumor_12week"=391:585,"Tumor_16week"=586:780,"Tumor_20week"=781:936)
days=factor(n3)
levels(days)=list(Control_8week=1:195,Tumor_8week=196:390,Tumor_12week=391:585,Tumor_16week=586:780,Tumor_20week=781:936)
groups=rep(c("ControlWeek8","TumorWeek8","TumorWeek12","TumorWeek16","TumorWeek20"),each=195)
groups=groups[1:936]
df2$groups=groups
df2$groups=as.character(df2$groups)
df2$groups=factor(df2$groups, levels=unique(df2$groups))
gd=df2 %>%                            #Calculate means to then plot
	group_by(df2$rowid,groups) %>%
	summarise(x=mean(df2$value))

pdf(file="JUNDclustertimepoints.pdf")
ggplot()+geom_point(data=df2,aes(groups,value,group=factor(df2$rowid)))+facet_wrap(~ df2$rowid)+ theme(text=element_text(size=8),axis.text.x=element_text(angle=45,hjust=1),legend.position="none")
dev.off()

ggplot(df2,aes(groups,value,group=factor(df2$rowid)))+facet_wrap(~ df2$rowid)+theme(text=element_text(size=8),axis.text.x=element_text(angle=45,hjust=1),legend.position="none")+
geom_point(data=gd)

+geom_point(data=df,aes(variable,value,group=factor(rowid),shape=15,size=2))
 
#Pull out BCAR3 data for Amy Bouton--------------
BCAR3=expwithGeneNames[which(expwithGeneNames$GeneName =='Bcar3'),]
write.table(BCAR3, file='BCAR3exp.txt', sep="\t")
#If multiple probes corresponded to the same gene, then the expression values of these probe sets were averaged
BCAR3avg <- aggregate(. ~ GeneName, data = BCAR3, mean)
write.table(BCAR3avg, file='BCAR3avgexp.txt', sep="\t")

#Looking at ER, PR, HER2 status
ESR1=expwithGeneNames[which(expwithGeneNames$GeneName =='Esr1'),]
PGR=expwithGeneNames[which(expwithGeneNames$GeneName =='Pgr'),]
ERBB2=expwithGeneNames[which(expwithGeneNames$GeneName =='Erbb2'),]
#write.table(BCAR3, file='BCAR3exp.txt', sep="\t")
#If multiple probes corresponded to the same gene, then the expression values of these probe sets were averaged
ESR1avg <- aggregate(. ~ GeneName, data = ESR1, mean)
PGRavg <- aggregate(. ~ GeneName, data = PGR, mean)
ERBB2avg <- aggregate(. ~ GeneName, data = ERBB2, mean)