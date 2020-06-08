# Set working directory
setwd("/Volumes/GoogleDrive/My Drive/Janes Lab/Lab meetings/2020-06-05 code")

# Load dataset (+ do some other things that aren't important here)
source("./import_opcBulk.R")

# Get only the samples (columns) of interest
rsem_counts <- bulk$rsem[,c(1:9,9+which(bulk$info$day == 150 & bulk$info$genotype == "WT"),9+which(bulk$info$day == 150 & bulk$info$genotype == "CKO"))]

# Visualize the "read-depth" (~= "sample loading")
readDepth <- colSums(rsem_counts[,10:ncol(rsem_counts)])
barplot(height = as.numeric(readDepth),ylab = "Number of reads",xlab = "13 samples")

# Transcript-per-million (TPM) normalization*
# rsem_counts[,10] <- rsem_counts[,10] / sum(rsem_counts[,10]) * 1e6
# rsem_counts[,11] <- rsem_counts[,11] / sum(rsem_counts[,11]) * 1e6
# rsem_counts[,12] <- rsem_counts[,12] / sum(rsem_counts[,12]) * 1e6

# ... or
source("./normalizeTPM.R")
tpm <- normalizeTPM(rsem = rsem_counts,index_counts = 10:ncol(rsem_counts))
readDepth2 <- colSums(tpm[,10:ncol(tpm)])
barplot(height = as.numeric(readDepth2),ylab = "Number of reads",xlab = "13 samples")

# Correlation of two wild-type samples
plot(x = tpm[,10],y = tpm[,11],xlab = "Sample #1 TPM",ylab = "Sample #2 TPM")

# Log2 transformation
log2_tpm <- tpm
log2_tpm[,10:ncol(log2_tpm)] <- log2(log2_tpm[,10:ncol(log2_tpm)] + 1)
# Plot log-transformed correlations
plot(x = log2_tpm[,10],y = log2_tpm[,11],xlab = "Sample #1 Log2(TPM+1)",ylab = "Sample #2 Log2(TPM+1)",pch = 16,cex = 0.5,col = "#00000044")



# DESeq2 implementaion -----

{
  rsem_counts <- bulk$rsem[,c(1:9,9+which(bulk$info$day == 150 & bulk$info$genotype == "WT"),9+which(bulk$info$day == 150 & bulk$info$genotype == "CKO"))]
  
  # Format for DESeq2
  row.names(rsem_counts) <- rsem_counts$symbol
  rsem_counts <- rsem_counts[,10:ncol(rsem_counts)]
  
  rsem_counts <- as.matrix(rsem_counts)
  mode(rsem_counts) <- "integer"
  
  # Filter out undetected genes
  rsem_counts <- rsem_counts[rowSums(rsem_counts) > 0,]
  
  # Create a dataframe with column information
  dataset_info <- data.frame(row.names = colnames(rsem_counts),
                             genotype = unlist(lapply(strsplit(colnames(rsem_counts),split = "_"),function(x) x[3])),
                             sex = c(as.character(bulk$info$sex[bulk$info$day == 150 & bulk$info$genotype == "WT"]),as.character(bulk$info$sex[bulk$info$day == 150 & bulk$info$genotype == "CKO"])))
  
  # Set the wild-type condition as the "baseline"
  dataset_info$genotype <- relevel(x = dataset_info$genotype,ref = "WT")
  
  # BiocManager::install("DESeq2")
  library(DESeq2)
  
  # Format dataset as a DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = rsem_counts,colData = dataset_info,design = ~genotype)
  
  # Run DESeq2
  dds <- DESeq(object = dds)
  
  # Extract results
  res <- results(object = dds)
  head(res)
  
  # Reorder from smallest -> largest p-value
  resOrdered <- res[order(res$padj),]
  print(resOrdered)
  
  # Remove "NA" rows
  resOrdered <- resOrdered[!is.na(resOrdered$padj),]
  
  # Get significantly DE genes
  resOrdered <- resOrdered[resOrdered$padj < 0.05,]
  nrow(resOrdered)
  
  # Save only gene names
  DEgenes_DESeq2 <- row.names(resOrdered)
  
  # Useful built-in plots
  DESeq2::plotMA(object = dds)
  DESeq2::plotCounts(dds = dds,gene = "Epha6",intgroup = "genotype")
  
  # Volcano plot
  plot(x = res$log2FoldChange,y = -log10(res$padj),
       xlim = c(-10,10),ylim = c(0,100),pch = 16,cex = 0.5,col = "#00000066",
       xlab = "Log2FoldChange",ylab = "-Log10(p-value)")
  
}

# edgeR implementation -----

{
  rsem_counts <- bulk$rsem[,c(1:9,9+which(bulk$info$day == 150 & bulk$info$genotype == "WT"),9+which(bulk$info$day == 150 & bulk$info$genotype == "CKO"))]
  
  # Format for edgeR
  row.names(rsem_counts) <- rsem_counts$symbol
  rsem_counts <- rsem_counts[,10:ncol(rsem_counts)]
  
  rsem_counts <- as.matrix(rsem_counts)
  mode(rsem_counts) <- "integer"
  
  # Filter out undetected genes
  rsem_counts <- rsem_counts[rowSums(rsem_counts) > 0,]
  
  # Create a dataframe with column information
  dataset_info <- data.frame(row.names = colnames(rsem_counts),
                             genotype = unlist(lapply(strsplit(colnames(rsem_counts),split = "_"),function(x) x[3])),
                             sex = c(as.character(bulk$info$sex[bulk$info$day == 150 & bulk$info$genotype == "WT"]),as.character(bulk$info$sex[bulk$info$day == 150 & bulk$info$genotype == "CKO"])))
  
  # Set the wild-type condition as the "baseline"
  dataset_info$genotype <- relevel(x = dataset_info$genotype,ref = "WT")
  
  # BiocManager::install("edgeR")
  library(edgeR)
  
  # Create edgeR object
  dgList <- DGEList(counts = rsem_counts,genes = row.names(rsem_counts))
  
  # edgeR's normalization method based on "effective library size"
  dgList <- calcNormFactors(dgList, method="TMM")
  
  # Create a study "design matrix"
  designMat <- model.matrix(~dataset_info$genotype)
  
  # Estimate dispersion for negative binomial model
  # Why negative binomial? https://bioramble.wordpress.com/2016/01/30/why-sequencing-data-is-modeled-as-negative-binomial/
  dgList <- estimateDisp(dgList, design = designMat)
  
  # Fit the generalized linear model to the data
  fit <- glmFit(y = dgList, design = designMat)
  
  # Likelihood ratio test of the model fit
  lrt <- glmLRT(glmfit = fit)
  topTags(object = lrt)
  head(lrt$table)
  
  # Determine significationly DE genes with multiple hypothesis test correction
  deGenes <- decideTestsDGE(lrt, p=0.05)
  
  # Get only gene names
  DEgenes_edgeR <- rownames(lrt)[as.logical(deGenes)]
  
}

# Comparison of DESeq2 and edgeR results -----

library(VennDiagram)

v <- venn.diagram(x = list("DESeq2" = DEgenes_DESeq2,"edgeR" = DEgenes_edgeR),
                  filename = "./venn_DESeq2_edgeR.png",height = 1600,width = 1600,units = "px",resolution = 300)


# Format data frame for overlap test
df <- data.frame(x = c(rep("DESeq2",length(DEgenes_DESeq2)),
                       rep("edgeR",length(DEgenes_edgeR))),
                 y = c(DEgenes_DESeq2,DEgenes_edgeR))
source("./overlap.test.R")
# Test significance of overlap
overlap.test(df = df,nGenes = 10000,nSim = 1000)

# (Unranked) gene set enrichment analysis -----
library(hypeR)

# Display available gene sets
msigdb_info()

# Download all mouse gene sets
# msigdb_download_all(species = "Mus musculus",output_dir = ".")

# Import Hallmark gene sets ("R Data Set file")
HALLMARK <- readRDS(file = "./H.v7.1.1.rds")

# Run enrichment
hyp_obj <- hypeR(signature = DEgenes_DESeq2,gsets = HALLMARK,fdr_cutoff = 0.05)

# View results
View(hyp_obj$as.data.frame())

# Visualize as heatmap -----
log2_tpm_subset <- log2_tpm[log2_tpm$symbol %in% DEgenes_DESeq2, 10:ncol(log2_tpm)]

library(pheatmap)
library(RColorBrewer) # https://colorbrewer2.org/

set.seed(0) # random seed, for repeatable "random" results
pheatmap(mat = as.matrix(log2_tpm_subset[sample(x = nrow(log2_tpm_subset),size = 500),]),
         show_rownames = F,
         # annotation_col = dataset_info,
         # color = rev(brewer.pal(n = 9,name = "RdBu"))
         )

# Center and scale on a gene-by-gene basis
log2_tpm_subset_scale <- t(scale(x = t(log2_tpm_subset)))

set.seed(0)
pheatmap(mat = as.matrix(log2_tpm_subset_scale[sample(x = nrow(log2_tpm_subset_scale),size = 500),]),
         show_rownames = F,
         annotation_col = dataset_info,
         color = rev(brewer.pal(n = 9,name = "RdBu")))

# Scale only for display
set.seed(0)
pheatmap(mat = as.matrix(log2_tpm_subset[sample(x = nrow(log2_tpm_subset),size = 500),]),
         show_rownames = F,
         annotation_col = dataset_info,
         color = rev(brewer.pal(n = 9,name = "RdBu")),
         scale = "row")





# Principle Components Analysis
log2_tpm_subset <- log2_tpm[log2_tpm$symbol %in% DEgenes_DESeq2,10:ncol(log2_tpm)]

# devtools::install_github("vqv/ggbiplot")
library(ggbiplot)

pc <- prcomp(x = t(log2_tpm_subset))
ggbiplot(pc,var.axes = F,groups = dataset_info$genotype)

# UMAP
library(uwot)

umap_res <- umap(X = t(log2_tpm_subset),n_neighbors = 4)
plot(x = umap_res[,1],y = umap_res[,2])

point_shape <- c(1,16)[as.numeric(dataset_info$genotype)]

plot(x = umap_res[,1],y = umap_res[,2],
     pch = point_shape)
