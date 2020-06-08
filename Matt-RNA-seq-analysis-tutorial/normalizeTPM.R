normalizeTPM <- function(rsem, index_counts) {
  # rsem: data frame with rsem values
  # index_counts: vector of indices with rsem counts
  
  # Remove MT genes and ERCCs
  # rsem_mt <- rsem[rsem$chr == "MT",]
  # rsem_ercc <- rsem[rsem$chr == "ERCC",]
  # rsem <-  rsem[rsem$chr != "MT" & rsem$chr != "ERCC",]
  
  # Calculate total reads for each sample
  totalreads <- colSums(rsem[,index_counts])
  
  # Normalize to TPM
  rsem[,index_counts] <- sweep(rsem[,index_counts]*1e6,2,totalreads,"/")
  
  # Scale MT genes and ERCCs by the same factor
  # if (nrow(rsem_mt) > 0) {
  #   rsem_mt[,index_counts] <- sweep(rsem_mt[,index_counts]*1e6,2,totalreads,"/")
  # }
  # if (nrow(rsem_ercc) > 0) {
  #   rsem_ercc[,index_counts] <- sweep(rsem_ercc[,index_counts]*1e6,2,totalreads,"/")
  # }
  
  # Re-insert 'normalized' MT genes and ERCCs
  # rsem <- rbind(rsem,rsem_mt,rsem_ercc)
  
  return(rsem)
}