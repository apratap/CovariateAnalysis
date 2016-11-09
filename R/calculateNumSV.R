# Function to calculate number of surrogate variables using "be" method
## Modified from jleek/sva repository num.sv.R
## Has option for parallel for loop by foreach package and mixed effect model using lme4 package
calculateNumSV <- function (dat, mod, block = NULL, method = c("be"), vfilter = NULL, B = 20, seed = NULL) {
  
  if (!is.null(vfilter)) {
    if (vfilter < 100 | vfilter > dim(dat)[1]) {
      stop(paste("The number of genes used in the analysis must be between 100 and", 
                 dim(dat)[1], "\n"))
    }
    tmpv = rowVars(dat)
    ind = which(rank(-tmpv) < vfilter)
    dat = dat[ind, ]
  }
  
  method <- match.arg(method)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (method == "be" && is.null(block)) {
    voom.expr <- limma::voom(dat, design=mod, plot=F)
    voom.fit <- limma::lmFit(voom.expr)
    res <- limma::residuals.MArrayLM(voom.fit, voom.expr$E)
    uu <- svd(res)
    ndf <- unique(voom.fit$df.residual)
    dstat <- uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
    dstat0 <- matrix(0, nrow = B, ncol = ndf)
    dstat0 <- foreach(i=1:B, .combine = cbind) %dopar% {
      res0 <- t(apply(res, 1, sample, replace = FALSE))
      voom.fit0 <- limma::lmFit(res0, design = mod)
      res0 <- limma::residuals.MArrayLM(voom.fit0, res0)
      uu0 <- svd(res0)
      tmp <- uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
    }
    dstat0 <- t(dstat0)
    psv <- rep(1, ncol(res))
    for (i in 1:ndf) {
      psv[i] <- mean(dstat0[, i] >= dstat[i])
    }
    for (i in 2:ndf) {
      psv[i] <- max(psv[(i - 1)], psv[i])
    }
    nsv <- sum(psv <= 0.1)
    return(nsv)
  }
  else if (method == "be" && !is.null(block)){
    voom.expr <- limma::voom(dat, design=mod, plot=F)
    correlation = CovariateAnalysis::parallelDuplicateCorrelation(voom.expr, block = block)
    voom.expr <- limma::voom(dat, design=mod, block = block, correlation = correlation$cor, plot=F)
    voom.fit <- limma::lmFit(voom.expr)
    res <- limma::residuals.MArrayLM(voom.fit, voom.expr$E)
    uu <- svd(res)
    ndf <- unique(voom.fit$df.residual)
    dstat <- uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
    dstat0 <- foreach(i=1:B, .combine = cbind) %dopar% {
      res0 <- t(apply(res, 1, sample, replace = FALSE))
      voom.fit0 <- limma::lmFit(res0, design=mod, block = block, correlation = correlation$cor)
      res0 <- limma::residuals.MArrayLM(voom.fit0, res0)
      uu0 <- svd(res0)
      tmp <- uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
    }
    dstat0 <- t(dstat0)
    psv <- rep(1, ncol(res))
    for (i in 1:ndf) {
      psv[i] <- mean(dstat0[, i] >= dstat[i])
    }
    for (i in 2:ndf) {
      psv[i] <- max(psv[(i - 1)], psv[i])
    }
    nsv <- sum(psv <= 0.1)
    return(nsv)
  }
}