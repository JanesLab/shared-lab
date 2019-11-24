removeDuplicates_improved <- function(rsem, index_counts) {
  # rsem: data frame with rsem values
  # index_counts: vector of indices with rsem counts

  duplicatedSymbols <- unique(rsem$symbol[duplicated(rsem$symbol)])
  indsToDelete <- numeric(0)
  for (iSymbol in duplicatedSymbols) {
    duplicatedSymbolRows <- which(rsem$symbol == iSymbol)
    
    # Save row that has the shortest chromosome name / lowest Ensembl ID
    indPriority <- order(nchar(rsem$chr[duplicatedSymbolRows]),rsem$ensgene[duplicatedSymbolRows])
    
    indToKeep <- duplicatedSymbolRows[indPriority[1]]
    indsToDelete <- c(indsToDelete,duplicatedSymbolRows[indPriority[2:length(indPriority)]])
    rsem[indToKeep,index_counts] <- colSums(rsem[duplicatedSymbolRows,index_counts])
  }
  rsem <- rsem[-indsToDelete,]
  
  return(rsem)
}