fixERCCs <- function(rsem) {
  # rsem: data frame with rsem values
  
  rsem$chr[substr(rsem$ensgene,1,5) == "ERCC-"] <- "ERCC"
  rsem$symbol[substr(rsem$ensgene,1,5) == "ERCC-"] <- rsem$ensgene[substr(rsem$ensgene,1,5) == "ERCC-"]
  
  return(rsem)
}