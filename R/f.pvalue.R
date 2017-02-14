f.pvalue <- function (dat, mod, mod0, block = NULL) {
  n <- dim(dat)[2]
  m <- dim(dat)[1]
  
  df1 <- dim(mod)[2]
  df0 <- dim(mod0)[2]
  
  p <- rep(0, m)
  Id <- diag(n)
  if (!is.null(block)){
    voom.expr <- limma::voom(dat, design=mod, plot=F)
    correlation = CovariateAnalysis::parallelDuplicateCorrelation(voom.expr, block = block)
    voom.expr <- limma::voom(dat, design=mod, block = block, correlation = correlation$cor, plot=F)
    voom.fit <- limma::lmFit(voom.expr)
    resid <- limma::residuals.MArrayLM(voom.fit, voom.expr$E)
  } else {
    voom.expr <- limma::voom(dat, design=mod, plot=F)
    voom.fit <- limma::lmFit(voom.expr)
    resid <- limma::residuals.MArrayLM(voom.fit, voom.expr$E)
  }
  rss1 <- rowSums(resid * resid)
  rm(resid)
  if (!is.null(block)){
    voom.expr <- limma::voom(dat, design=mod, plot=F)
    correlation = CovariateAnalysis::parallelDuplicateCorrelation(voom.expr, block = block)
    voom.expr <- limma::voom(dat, design=mod, block = block, correlation = correlation$cor, plot=F)
    voom.fit <- limma::lmFit(voom.expr)
    resid0 <- limma::residuals.MArrayLM(voom.fit, voom.expr$E)
  } else {
    voom.expr <- limma::voom(dat, design=mod, plot=F)
    voom.fit <- limma::lmFit(voom.expr)
    resid0 <- limma::residuals.MArrayLM(voom.fit, voom.expr$E)
  }
  rss0 <- rowSums(resid0 * resid0)
  rm(resid0)
 
  fstats <- ((rss0 - rss1)/(df1 - df0))/(rss1/(n - df1))
  p <- 1 - pf(fstats, df1 = (df1 - df0), df2 = (n - df1))
  return(p)
}