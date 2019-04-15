
# install.packages("lme4")
# 
library(BiocManager)
# BiocManager::install("tximport")
# 
# install.packages("devtools")
devtools::install_github("stephenturner/annotables")

# Load modules
library(Tmisc)
library(annotables)
library(tidyverse)
library(tximport)

# Collapse the annotation table to the ensembl gene id (mouse)
annotables::grcm38 %>% distinct(ensgene)
anno_collapsed_to_ensgene <- annotables::grcm38 %>%
  group_by(ensgene) %>%
  summarize_all(funs(. %>% unique %>% paste(collapse="; ")))

sampnames <- list.files(pattern="isoforms.results") %>% str_replace(".isoforms.results","")
quantfiles <- list.files(pattern="isoforms.results") %>% set_names(sampnames)

# # Join mouse genes and ERCCs
tx2gene <- annotables::grcm38_tx2gene
ercctx2gene <- read_tsv("/Users/ss4cf/Desktop/Lab/Data/BCC\ Samples/BCC_RNAseq/RSEM\ data/ercc_tx2gene.tsv")
tx2gene <- rbind(tx2gene, ercctx2gene)

# # edit above for human
# annotables::grch38 %>% distinct(ensgene)
# anno_collapsed_to_ensgene <- annotables::grch38 %>% 
#   group_by(ensgene) %>% 
#   summarize_all(funs(. %>% unique %>% paste(collapse="; ")))

# file names for data to import
setwd('/Users/ss4cf/Desktop/Lab/Dylan\ Schaff/RNAseq/KP1_impos')
sampnames <- list.files(pattern="isoforms.results") %>% str_replace(".isoforms.results","")
quantfiles <- list.files(pattern="isoforms.results") %>% set_names(sampnames)

# # Join human genes and ERCCs
# tx2gene <- annotables::grch38_tx2gene
# ercctx2gene <- read_tsv("~/Desktop/Lab/Data/BCC Samples/BCC_RNAseq/ercc_tx2gene.txt")
# tx2gene <- rbind(tx2gene, ercctx2gene)


# Import transcript data

txi <- tximport(quantfiles,
                type = "rsem",
                tx2gene = tx2gene,
                importer = read_tsv,
                txIn = TRUE,
                txOut = FALSE,
                ignoreTxVersion = TRUE, 
                countsFromAbundance="no") 

# Annotate genes and save .csv file
rsem.counts <- txi$counts %>% as.data.frame()
rsem.counts$ensgene <- rownames(rsem.counts)
rsem.counts <- left_join(rsem.counts, anno_collapsed_to_ensgene) %>%
  select(ensgene, entrez, symbol, chr, start, end, strand, biotype, description, everything()) %>%
  write_csv("janes_dataset8_kp1_impos_rsem_counts_both.csv")

#compare both
temp <- read.csv('/Users/ss4cf/Desktop/Lab/Dylan\ Schaff/RNAseq/OLDjanes_dataset7_rsem_counts.csv',header=TRUE,sep=',',stringsAsFactors = F)

rsemindex <- 10:57
mets <- process_rsem(rsem.counts,rsemindex)
spheres <- process_rsem(temp,rsemindex)

nGenes_mets <- count_genes(mets,rsemindex,10)
nGenes_spheres <- count_genes(spheres,rsemindex,10)
