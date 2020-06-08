# Generalized statistical test for gene list overlaps
# Does not test "exclusive"/"exact" overlaps (e.g. genes overlapping in groups 1 and 2 that do NOT show up in group 3)

# Format of df
# df
#     x             y
# 1   geneList1     gene1
# 2   geneList1     gene2
# 3   geneList1     gene3
# 4   geneList1     gene4
# 5   ...           ...
# 6   geneList2     gene1
# 7   geneList2     gene2
# 8   geneList2     gene3
# 9   geneList2     gene4

calculateIntersections <- function(df) {
  intersections <- rep(0,2^length(levels(df$x))-length(levels(df$x))-1)
  # calculate intersections
  iIntersection <- 1
  for (iCombn in 2:length(levels(df$x))) {
    comparisons <- combn(levels(df$x),iCombn)
    for (iComparison in 1:ncol(comparisons)) {
      appearances <- table(df$y[df$x %in% comparisons[,iComparison]])
      intersectGenes <- sum(appearances >= nrow(comparisons))
      intersections[iIntersection] <- intersectGenes
      iIntersection <- iIntersection + 1
    }
  }
  return(intersections)
}

overlap.test <- function(df,nGenes,nSim) {
  # Calculate true intersections
  intersections <- calculateIntersections(df)
  
  # Create vector to store monte carlo simulations
  sigList <- matrix(0,nrow = length(intersections),ncol = nSim)
  
  # Make a "gene list" for the simulations
  Genes <- 1:nGenes
  
  # Calculate number of genes to draw for each group
  groups <- as.numeric(table(df$x))
  
  # Run simulations
  for (iSim in 1:nSim) {
    idf <- data.frame(x = numeric(),y = numeric())
    # For each group, randomly sample genes
    for (iGroup in 1:length(groups)) {
      iGenes <- sample(Genes,groups[iGroup])
      idf <- rbind(idf,data.frame(x = rep(iGroup,length(iGenes)),
                                  y = iGenes))
    }
    idf$x <- factor(idf$x)
    
    # Calculate simulated intersection
    sigList[,iSim] <- calculateIntersections(idf)
  }
  
  # Calculate p values
  p <- rep(1,length(intersections))
  for (iComparison in 1:length(intersections)) {
    p[iComparison] <- sum(sigList[iComparison,] >= intersections[iComparison]) / nSim
  }
  
  # Store results in data frame
  iIntersection <- 1
  res <- data.frame(comparison = rep("",length(intersections)),intersection = intersections,p.value = p,stringsAsFactors = F)
  for (iCombn in 2:length(levels(df$x))) {
    comparisons <- combn(levels(df$x),iCombn)
    strs <- apply(comparisons,2,function(x) paste(x,collapse=" | "))
    res$comparison[iIntersection:(iIntersection+length(strs)-1)] <- strs
    iIntersection <- iIntersection + length(strs)
  }
  return(res)
}