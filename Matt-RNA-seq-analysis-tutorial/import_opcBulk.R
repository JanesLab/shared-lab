# Import bulk samples

import_opcBulk <- function() {
  
  source("./normalizeTPM.R")
  
  # File paths
  bulk.path <- "./rsem_opcBulk.csv"
  bulk.info.path <- "./info_opcBulk.csv"
  
  # Import
  bulk <- read.csv(bulk.path,stringsAsFactors = F)
  bulk.info <- read.csv(bulk.info.path,stringsAsFactors = F)
  
  # Keep only tdT+ samples
  bulk <- bulk[,c(1:9,9+which(bulk.info$celltype == "tdTpositive"))]
  bulk.info <- bulk.info[bulk.info$celltype == "tdTpositive",]
  
  # Remove outliers
  outliers <- c("opcBulk_90dpi_CKO_20647_tdTpositive_GAGCT","opcBulk_90dpi_WT_21167_tdTpositive")
  bulk <- bulk[,!(names(bulk) %in% outliers)]
  bulk.info <- bulk.info[!(bulk.info$name %in% outliers),]
  
  # Rename samples if they're from the same mouse as an outlier
  for (iOutlier in outliers) {
    if (length(strsplit(x = iOutlier,split = "_")[[1]]) == 6) {
      iRootName <- paste(strsplit(x = iOutlier,split = "_")[[1]][-6],collapse = "_")
      names(bulk)[which(grepl(iRootName,names(bulk)))] <- iRootName
      bulk.info$name[which(grepl(iRootName,bulk.info$name))] <- iRootName
    }
  }
  
  bulk_collapse <- bulk
  bulk_collapse.info <- bulk.info
  
  # Convert all columns to characters
  bulk_collapse.info[,names(bulk_collapse.info)] <- lapply(bulk_collapse.info[,names(bulk_collapse.info)],as.character)
  
  # Collapse replicates from same mouse
  bulk_collapse_unique_replicateIndex.info <- unique(bulk_collapse.info[,c("day","genotype","sex","brainID")])
  indsToRemove <- c()
  for (i in 1:nrow(bulk_collapse_unique_replicateIndex.info)) {
    inds <- which(bulk_collapse.info$day == bulk_collapse_unique_replicateIndex.info$day[i] &
                    bulk_collapse.info$genotype == bulk_collapse_unique_replicateIndex.info$genotype[i] &
                    bulk_collapse.info$sex == bulk_collapse_unique_replicateIndex.info$sex[i] &
                    bulk_collapse.info$brainID == bulk_collapse_unique_replicateIndex.info$brainID[i])
    if (length(inds) > 1) {
      bulk_collapse[,9+inds[1]] <- rowSums(bulk_collapse[,9+inds])
      bulk_collapse.info[inds[1],"nCells"] <- paste(bulk_collapse.info[inds,"nCells"],collapse="|")
      bulk_collapse.info[inds[1],"runID"] <- paste(bulk_collapse.info[inds,"runID"],collapse="|")
      bulk_collapse.info[inds[1],"barcode"] <- paste(bulk_collapse.info[inds,"barcode"],collapse="|")
      bulk_collapse.info[inds[1],"replicate"] <- paste(bulk_collapse.info[inds,"replicate"],collapse="|")
      indsToRemove <- c(indsToRemove,inds[2:length(inds)])
      
      iNewName <- paste(strsplit(x = names(bulk_collapse)[9+inds[1]],split = "_")[[1]][-6],collapse = "_")
      names(bulk_collapse)[9+inds[1]] <- iNewName
      bulk_collapse.info$name[inds[1]] <- iNewName
    }
  }
  bulk_collapse <- bulk_collapse[,-(9+indsToRemove)]
  bulk_collapse.info <- bulk_collapse.info[-indsToRemove,]
  
  # Convert all columns to factors
  bulk_collapse.info[,names(bulk_collapse.info)] <- lapply(bulk_collapse.info[,names(bulk_collapse.info)],factor)
  
  bulk_collapse_tpm <- normalizeTPM(rsem = bulk_collapse,index_counts = 10:ncol(bulk_collapse))
  
  bulk_collapse_tpm_log2 <- cbind(bulk_collapse_tpm[,1:9],log2(bulk_collapse_tpm[,10:ncol(bulk_collapse_tpm)] + 1))
  
  return(list(rsem = bulk_collapse,
              info = bulk_collapse.info,
              tpm = bulk_collapse_tpm,
              log2 = bulk_collapse_tpm_log2))
}

bulk <- import_opcBulk()