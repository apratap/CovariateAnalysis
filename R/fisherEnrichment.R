# Function to perform Fishers enrichment analysis
fisherEnrichment <- function(genesInSignificantSet, # A character vector of differentially expressed or some significant genes to test
                             genesInGeneSet, # A character vector of genes in gene set like GO annotations, pathways etc...
                             genesInBackground # Background genes that are 
){
  genesInSignificantSet = intersect(genesInSignificantSet, genesInBackground) # back ground filtering
  genesInNonSignificantSet = base::setdiff(genesInBackground, genesInSignificantSet)
  genesInGeneSet = intersect(genesInGeneSet, genesInBackground) # back ground filtering
  genesOutGeneSet = base::setdiff(genesInBackground,genesInGeneSet)
  
  pval = fisher.test(
    matrix(c(length(intersect(genesInGeneSet, genesInSignificantSet)),             
             length(intersect(genesInGeneSet, genesInNonSignificantSet)),
             length(intersect(genesOutGeneSet, genesInSignificantSet)),
             length(intersect(genesOutGeneSet, genesInNonSignificantSet))), 
           nrow=2, ncol=2),
    alternative="greater")
  OR = (length(intersect(genesInGeneSet, genesInSignificantSet)) * length(intersect(genesOutGeneSet, genesInNonSignificantSet))) / (length(intersect(genesInGeneSet, genesInNonSignificantSet)) * length(intersect(genesOutGeneSet, genesInSignificantSet)))
  return(data.frame(pval = pval$p.value,
                    ngenes = length(genesInGeneSet),
                    noverlap = length(intersect(genesInGeneSet, genesInSignificantSet)),
                    Odds.Ratio = OR,
                    Genes = paste(intersect(genesInGeneSet, genesInSignificantSet), collapse = '|')
  )
  )
}
