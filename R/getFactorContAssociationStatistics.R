# Find Inter Class Correlation between factor and continuous covariates
# Inspired from http://stats.stackexchange.com/questions/108007/correlations-with-categorical-variables
getFactorContAssociationStatistics <- function(factorContNames,COVARIATES, na.action='remove',
                                               alpha = 0.05){

  factorContNames = factorContNames %>% unlist %>% as.character()

  if (na.action == "remove")
    COVARIATES = na.omit(COVARIATES[,factorContNames])

  class.names = sapply(COVARIATES, is.factor)
  stats = ICC(colnames(COVARIATES)[class.names], colnames(COVARIATES)[!class.names], COVARIATES, alpha = alpha)

  return(data.frame(Estimate = stats$ICC, Pval = stats$p.val))
}
