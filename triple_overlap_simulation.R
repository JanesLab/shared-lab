# Overlap of three groups

# Total number of genes
nTotalGenes <- 31581

# Number of significant genes in each group
nGroup1 <- 587
nGroup2 <- 126
nGroup3 <- 6660

# Number of simulations to run
nSimulations <- 10000

# -----

# Use a vector of numbers/indices instead of actual gene names
geneList <- 1:nTotalGenes

# Create vector that will hold the number of overlapping genes for each simulation
sigAll <- rep(0,nSimulations)
sigGroups_1_2 <- rep(0,nSimulations) # Includes the triple overlapping genes
sigGroups_1_3 <- rep(0,nSimulations) # Includes the triple overlapping genes
sigGroups_2_3 <- rep(0,nSimulations) # Includes the triple overlapping genes

# Run simulations
for (iSim in 1:nSimulations) {
  # Randomly sample from geneList
  geneList1 <- sample(geneList,nGroup1)
  geneList2 <- sample(geneList,nGroup2)
  geneList3 <- sample(geneList,nGroup3)
  
  # Calculate overlaps
  sigAll[iSim] <- length(intersect(geneList1,intersect(geneList2,geneList3))) # can only intersect 2 groups at a time
  sigGroups_1_2[iSim] <- length(intersect(geneList1,geneList2))
  sigGroups_1_3[iSim] <- length(intersect(geneList1,geneList3))
  sigGroups_2_3[iSim] <- length(intersect(geneList2,geneList3))
}

# Visualize expected distribution of triple overlap
library(ggplot2)
ggplot(data = data.frame(x = sigAll),aes(x=x)) + geom_histogram() + xlab("Number of triple-overlapping genes")

# Calculate p value for triple overlap
observedOverlap <- 9
sum(sigAll >= observedOverlap) / nSimulations
