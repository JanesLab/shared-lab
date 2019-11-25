removeDuplicates_improved <- function(rsem, index_counts) {
  # rsem: data frame with rsem values
  # index_counts: vector of indices with rsem counts

  duplicatedSymbols <- unique(rsem$symbol[duplicated(rsem$symbol)])
  indsToDelete <- numeric(0)
  indsToKeep <- numeric(0)
  # Pre-allocate matrix for colSum-ed rows
  editVals <- matrix(data = rep(0,length(index_counts)*length(duplicatedSymbols)),nrow = length(duplicatedSymbols),ncol = length(index_counts))
  for (iSymbol in 1:length(duplicatedSymbols)) {
    duplicatedSymbolRows <- which(rsem$symbol == duplicatedSymbols[iSymbol])
    
    # Save row that has the shortest chromosome name / lowest Ensembl ID
    indPriority <- order(nchar(rsem$chr[duplicatedSymbolRows]),rsem$ensgene[duplicatedSymbolRows])
    
    indsToKeep <- c(indsToKeep,duplicatedSymbolRows[indPriority[1]])
    indsToDelete <- c(indsToDelete,duplicatedSymbolRows[indPriority[2:length(indPriority)]])
    editVals[iSymbol,] <- colSums(rsem[duplicatedSymbolRows,index_counts])
  }
  rsem[indsToKeep,index_counts] <- editVals
  rsem <- rsem[-indsToDelete,]
  
  return(rsem)
}